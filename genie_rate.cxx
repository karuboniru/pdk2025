#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TSpline.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <generator>
#include <iostream>
#include <print>
#include <ranges>
#include <sstream>

#include <boost/program_options.hpp>

struct config {
  std::string flux_file;
  std::string genie_xsec_file;
  double emin{}, emax{};
};

config parse_command_line(int argc, char **argv) {
  namespace po = boost::program_options;
  config cfg;

  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
    ("help,h", "Display this help message")
    ("flux-file,f", po::value<std::string>(&cfg.flux_file)->required(), "Input flux file")
    ("genie-xsec-file,g", po::value<std::string>(&cfg.genie_xsec_file)->required(), "Input GENIE cross section ROOT file")
    ("emin,i", po::value<double>(&cfg.emin)->default_value(0.1), "Minimum energy for integration (GeV)")
    ("emax,a", po::value<double>(&cfg.emax)->default_value(25.0), "Maximum energy for integration (GeV)")
    ;
  // clang-format on

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "Error: " << e.what() << "\n";
    std::cerr << desc << "\n";
    exit(EXIT_FAILURE);
  }

  if (vm.contains("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }

  return cfg;
}

constexpr auto neutrino_name =
    std::to_array({"nu_mu", "nu_mu_bar", "nu_e", "nu_e_bar"});
constexpr auto target_name = std::to_array({"H1", "O16"});
constexpr auto target_A = std::to_array({1, 16});
constexpr auto target_fraction =
    std::to_array({2. / 18., 16. / 18.}); // in nucleon number (H2O)
constexpr auto channel_name = std::to_array({"tot_cc", "tot_nc"});

constexpr double seconds_in_year = 3600 * 365 * 24;
constexpr double avogadro_number = 6.02214076e23;
constexpr double n_nucleon_per_Mton = 1e12 / 18.01 * avogadro_number * 18.0;
constexpr double var_4pi = 4 * M_PI;
constexpr double exposure_factor =
    seconds_in_year * n_nucleon_per_Mton * var_4pi;
constexpr double cross_section_unit = 1e-42; // code to m2
constexpr double factor = cross_section_unit * exposure_factor;

std::generator<std::pair<double, std::array<double, 4>>>
read_flux(std::string filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::println(std::cerr, "Error: Could not open file {}", filename);
    co_return;
  }
  for (std::string line; std::getline(file, line);) {
    if (line.empty() || line[0] == '#' || line[1] == 'E' || line[0] == 'a')
      continue; // Skip empty lines and comments
    std::pair<double, std::array<double, 4>> ret{};
    std::istringstream iss(line);
    if (iss >> ret.first >> ret.second[0] >> ret.second[1] >> ret.second[2] >>
        ret.second[3]) {
      co_yield ret;
    } else {
      std::println(std::cerr, "Error: Invalid line format in file {} :\n\t {}",
                   filename, line);
      co_return;
    }
  }
}

TSpline3 make_flux_spline(const std::vector<double> &log_energy,
                          const std::vector<double> &value) {
  TGraph graph(log_energy.size());
  for (size_t i = 0; i < log_energy.size(); ++i) {
    graph.SetPoint(i, log_energy[i], value[i]);
  }
  TSpline3 spline("flux_spline", &graph);
  return spline;
}

TSpline3 from_tgraph(const TGraph *graph) {
  TSpline3 spline("", graph);
  return spline;
}

TF1 make_flux_function(const std::vector<double> &log_energy,
                       const std::vector<double> &value, double emin,
                       double emax) {
  return TF1{"",
             [spline = make_flux_spline(log_energy, value)](const double *x,
                                                            const double *) {
               return spline.Eval(std::log10(x[0]));
             },
             emin, emax, 0};
}

int main(int argc, char **argv) {
  auto &&[flux_file, genie_xsec_file, emin, emax] =
      parse_command_line(argc, argv);
  std::vector<double> energy;
  std::array<std::vector<double>, 4> flux_values;
  for (const auto &[E, flux] :
       read_flux(flux_file) | std::views::take_while([&](const auto &d) {
         return d.first < emax * 2;
       })) {
    energy.push_back(std::log10(E));
    for (size_t i = 0; i < 4; ++i) {
      flux_values[i].push_back(flux[i]);
    }
  }

  auto flux_funcs = flux_values |
                    std::views::transform(
                        [&energy, &emin, &emax](const auto &values) mutable {
                          return make_flux_function(energy, values, emin, emax);
                        }) |
                    std::ranges::to<std::vector>();

  auto genie_cross_section_file = std::unique_ptr<TFile, void (*)(TFile *)>{
      TFile::Open(genie_xsec_file.c_str(), "READ"), [](TFile *f) {
        if (f) {
          f->Close();
          delete f;
        }
      }};

  if (!genie_cross_section_file || genie_cross_section_file->IsZombie()) {
    std::println(std::cerr, "Error: Could not open GENIE cross section file");
    return 1;
  }

  auto get_obj = [&](const std::string &name) {
    auto obj = genie_cross_section_file->Get<TGraph>(name.c_str());
    if (!obj) {
      throw std::runtime_error("Could not find object: " + name);
    }
    return obj;
  };

  auto xsec_func_ccnc =
      std::views::cartesian_product(neutrino_name, channel_name) |
      std::views::transform([&](auto &&names) {
        const auto &[nu_name, chan_name] = names;
        return TF1{
            "",
            [xsec_data =
                 std::views::zip(
                     target_name | std::views::transform([&](auto &&t_name) {
                       return std::format("{}_{}/{}", nu_name, t_name,
                                          chan_name);
                     }) | std::views::transform([&](auto &&full_name) {
                       return from_tgraph(get_obj(full_name));
                     }),
                     target_fraction, target_A) |
                 std::ranges::to<std::vector>()](const double *x,
                                                 const double *) {
              auto E = x[0];
              return std::ranges::fold_left_first(
                         xsec_data |
                             std::views::transform([E](const auto &data_frac) {
                               const auto &[graph, frac, A] = data_frac;
                               double xs = graph.Eval(E);
                               return xs * frac / A;
                             }),
                         std::plus<>{})
                  .value();
            },
            emin, emax, 0};
      });

  auto result =
      std::views::zip(std::views::cartesian_product(flux_funcs, channel_name),
                      xsec_func_ccnc) |
      std::views::transform([&](const auto &pair) {
        const auto &[flux, xsec] = pair;
        auto &&[flux_func, chan_name] = flux;
        TF1 rate_func{"",
                      [flux_func, xsec](const double *x, const double *) {
                        return flux_func.Eval(x[0]) * xsec.Eval(x[0]);
                      },
                      emin, emax, 0};
        rate_func.SetNpx(10000000);
        auto res = rate_func.Integral(emin, emax, 1e-10);
        auto fluxint = flux_func.Integral(emin, emax);
        auto xsec_avg = res / fluxint;
        return std::make_pair(res * factor, xsec_avg);
      }) |
      std::ranges::to<std::vector>();
  std::println("Flux integrated exposure per year per Mton events in range "
               "[{}, {}]: ",
               emin, emax);
  for (const auto &pair : std::views::zip(
           std::views::cartesian_product(neutrino_name, channel_name),
           result)) {
    const auto &[names, value] = pair;
    const auto &[nu_name, chan_name] = names;
    std::println("\t{:<9} {:<6}: {:.6e} with xsec = {} ", nu_name, chan_name,
                 value.first, value.second);
  }

  std::println("total: {:.6e}",
               std::ranges::fold_left_first(result | std::views::elements<0>,
                                            std::plus<>{})
                   .value());

  std::ranges::for_each(
      std::views::zip(neutrino_name,
                      flux_funcs | std::views::transform([&](auto &f) {
                        f.SetNpx(100000);
                        return f.Integral(emin, emax);
                      })),
      [](const auto &pair) {
        const auto &[chan_name, flux_int] = pair;
        std::println("Flux integral for channel {} : {}", chan_name, flux_int);
      });
}
