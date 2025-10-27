#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <TROOT.h>
#include <array>
#include <format>
#include <print>
#include <ranges>
#include <regex>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPie.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "EvtTracker2event.h"
#include "cmdline.h"
#include "commondefine.h"
#include "event.h"
#include "smear.h"

ROOT::RDF::RResultPtr<TH1D> make_plot(auto df, ROOT::RDF::TH1DModel model,
                                      const std::string &varname,
                                      const std::string &prefix = "") {
  model.fName = prefix + varname;
  model.fTitle = prefix + varname;
  return df.Histo1D(model, varname);
}

std::string normalize_name(std::string orig) {
  std::replace(orig.begin(), orig.end(), ' ', '_');
  orig = std::regex_replace(orig, std::regex("\\-"), "_minus");
  orig = std::regex_replace(orig, std::regex("\\+"), "_positive");
  return "h_" + orig;
}

auto iter_pair(auto &iter_range) {
  auto en = iter_range | std::views::enumerate;
  return std::views::cartesian_product(en, en) |
         std::views::filter([](const auto &pair) {
           auto [first, second] = pair;
           return std::get<0>(first) < std::get<0>(second);
         }) |
         std::views::transform([](const auto &pair) {
           auto [first, second] = pair;
           return std::make_pair(std::get<1>(first), std::get<1>(second));
         });
}

void make_pie_plot(auto &data, const std::string &filename) {
  auto col =
      std::to_array({kP6Blue, kP6Yellow, kP6Red, kP6Grape, kP6Gray, kP6Violet});
  auto pie = std::make_unique<TPie>("final state", "final state", data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    pie->SetEntryVal(i, data[i].second);
    pie->SetEntryFillColor(i, col[i % col.size()]);
    pie->SetEntryFillStyle(i, 1000 + (i / col.size()));
    pie->SetEntryLabel(i, data[i].first.c_str());
  }
  auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 700, 700);
  canvas->cd();
  // pie->SetRadius(0.25);
  pie->SetCircle(0.5, 0.5 - .1, 0.35);
  auto leg = pie->MakeLegend(.6, .6, .9, .9);
  pie->SetLabelFormat("%perc");
  pie->SetLabelsOffset(-.2);
  pie->Draw("");
  leg->Draw("SAME");
  canvas->SaveAs(filename.c_str());
}

pair_momentum_t operator+(const pair_momentum_t &a, const pair_momentum_t &b) {
  return {a.first + b.first, a.second + b.second};
}

int main(int argc, char **argv) {
  constexpr double to_deg = 180. / M_PI;
  initializeGaussianSmearStrategy();
  ROOT::EnableImplicitMT(4);
  TH1::AddDirectory(false);
  auto [input_files, output_path] = parse_command_line(argc, argv);

  auto tracker_df = TrackerPrepare(ROOT::RDataFrame{"outtree", input_files});
  ROOT::RDF::Experimental::AddProgressBar(tracker_df);
  auto df_all =
      tracker_df
          .Define("raw_proton_momentum",
                  [](const NeutrinoEvent &event) {
                    auto proton = event.in_range(2212);
                    return proton.begin()->second.P();
                  },
                  {"EventRecord"})
          .Define("raw_mass_proton",
                  [](const NeutrinoEvent &event) {
                    auto proton = event.in_range(2212);
                    return proton.begin()->second.M();
                  },
                  {"EventRecord"})
          .Define("final_state_system",
                  [](const NeutrinoEvent &event) {
                    const auto &post = event.get_post();
                    ROOT::Math::PxPyPzEVector total_momentum;
                    for (const auto &p :
                         post | std::views::filter([](const auto &entry) {
                           return entry.first != 2212 && entry.first != 2112;
                         }) | std::views::values) {
                      total_momentum += p;
                    }
                    return total_momentum;
                  },
                  {"EventRecord"})
          .Define("raw_pi0_before_fsi",
                  [](const NeutrinoEvent &event) {
                    auto pi0 = event.out_range(111);
                    return pi0.begin()->second.P();
                  },
                  {"EventRecord"})
          .Define("raw_final_state_mass",
                  [](const ROOT::Math::PxPyPzEVector &final_state_system) {
                    return final_state_system.M();
                  },
                  {"final_state_system"})
          .Define("raw_final_state_p",
                  [](const ROOT::Math::PxPyPzEVector &final_state_system) {
                    return final_state_system.P();
                  },
                  {"final_state_system"})
          .Define(
              "channel_name",
              [](NeutrinoEvent &e) { return e.get_channelname_no_nucleon(); },
              {"EventRecord"});
  auto count_per_channel = df_all.Aggregate(
      [](std::map<std::string, std::size_t> &data, const std::string &col) {
        data[col]++;
      },
      [](std::vector<std::map<std::string, std::size_t>> &to_merge) {
        for (auto &target = to_merge[0];
             const auto &item : to_merge | std::views::drop(1)) {
          for (const auto &[key, value] : item) {
            target[key] += value;
          }
        }
      },
      "channel_name", std::map<std::string, std::size_t>{});
  auto make_hist = []() { return TH1D("", "", 400, 0, 1.2); };
  auto pi0p_per_channel =
      df_all
          .Define("tmp_data_",
                  [](const std::string &channel, double p) {
                    return std::make_pair(channel, p);
                  },
                  {"channel_name", "raw_pi0_before_fsi"})
          .Aggregate(
              [&](std::map<std::string, TH1D> &data,
                  const std::pair<std::string, double> &col) {
                auto iter = data.find(col.first);
                if (iter == data.end()) {
                  iter = data.emplace(col.first, make_hist()).first;
                }
                iter->second.Fill(col.second);
              },
              [&](std::vector<std::map<std::string, TH1D>> &to_merge) {
                auto &target = to_merge[0];
                for (auto &item : to_merge | std::views::drop(1)) {
                  for (auto &[key, hist] : item) {
                    auto iter = target.find(key);
                    if (iter == target.end()) {
                      iter = target.emplace(key, make_hist()).first;
                    }
                    iter->second.Add(&hist);
                  }
                }
              },
              "tmp_data_", std::map<std::string, TH1D>{});

  auto event_count = df_all.Count();

  auto df_epi_final_state =
      df_all
          .Filter(
              [](const NeutrinoEvent &event) {
                return event.count_det(-11) == 1 && event.count_det(22) == 2 &&
                       event.count_det(211) == 0 && event.count_det(-211) == 0;
              },
              {"EventRecord"}, "epi final state selection")
          .Define("electron",
                  [](const NeutrinoEvent &event) {
                    auto electron = event.det_range(-11);
                    return electron.begin()->second;
                  },
                  {"EventRecord"})
          .Define("photons",
                  [](const NeutrinoEvent &event)
                      -> std::pair<pair_momentum_t, pair_momentum_t> {
                    auto photons = event.det_range(22);
                    if (photons.size() < 2) {
                      throw std::runtime_error(
                          "Not enough photons in the event");
                    }
                    auto photon1 = photons.begin()->second;
                    auto photon2 = (++photons.begin())->second;

                    if (photon1.second.P() < photon2.second.P()) {
                      // return std::make_pair(photon1, photon2);
                      std::swap(photon1, photon2);
                    }
                    return std::make_pair(photon1, photon2);
                  },
                  {"EventRecord"})
          .Define(
              "lead_photon",
              [](const std::pair<pair_momentum_t, pair_momentum_t> &photons) {
                return photons.first;
              },
              {"photons"})
          .Define(
              "sublead_photon",
              [](const std::pair<pair_momentum_t, pair_momentum_t> &photons) {
                return photons.second;
              },
              {"photons"})
          .Define("pi0_system",
                  [](const pair_momentum_t &leading_photon,
                     const pair_momentum_t &subleading_photon) {
                    return leading_photon + subleading_photon;
                  },
                  {"lead_photon", "sublead_photon"})
          .Define("epi_system",
                  [](const pair_momentum_t &electron,
                     const pair_momentum_t &pi0_system) {
                    return electron + pi0_system;
                  },
                  {"electron", "pi0_system"});

  auto &&[df_epi_with_vars, to_snapshot, mass_list, p_list] =
      DefineForEPi(ROOT::RDF::RNode{df_epi_final_state});

  std::ranges::copy(std::to_array({"raw_proton_momentum", "raw_mass_proton",
                                   "raw_pi0_before_fsi"}),
                    std::back_inserter(to_snapshot));
  std::ranges::copy(
      std::to_array({"raw_pi0_before_fsi", "raw_proton_momentum"}),
      std::back_inserter(p_list));
  std::ranges::copy(std::to_array({"raw_mass_proton"}),
                    std::back_inserter(mass_list));

  // auto [filtered_signal, upper, lower] =
  auto signals = FilterSignalKinematics(df_epi_with_vars);
  auto &&[filtered_signal, upper, lower] = signals;
  auto cut_efficiency = filtered_signal.Report();
  auto cut_efficiency_lower = lower.Report();
  auto cut_efficiency_upper = upper.Report();

  auto all_nofsi = df_epi_with_vars.Filter(
      [](const NeutrinoEvent &event) { return event.is_transparent(); },
      {"EventRecord"}, "transparent FSI  (effectively no FSI)");

  ROOT::RDF::TH1DModel inv_mass_model{"inv_mass_epip_system", "inv mass", 400,
                                      0.0, 1.2};
  ROOT::RDF::TH1DModel momentum_model{"inv_mass_epip_system", "momentum", 400,
                                      0.0, 1.2};
  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms{};

  for (const auto &varname :
       std::to_array({"raw_mass_proton", "raw_final_state_mass"})) {
    histograms.emplace_back(make_plot(df_all, inv_mass_model, varname));
  }

  for (const auto &varname :
       std::to_array({"raw_final_state_p", "raw_proton_momentum",
                      "raw_pi0_before_fsi"})) {
    histograms.emplace_back(make_plot(df_all, momentum_model, varname));
    histograms.emplace_back(
        make_plot(all_nofsi, momentum_model, varname, "noint_"));
  }

  for (auto &p_var : p_list) {
    histograms.emplace_back(
        make_plot(df_epi_with_vars, momentum_model, p_var, "epi_"));
    histograms.emplace_back(
        make_plot(all_nofsi, momentum_model, p_var, "noint_"));
  }

  for (auto &m_var : mass_list) {
    histograms.emplace_back(
        make_plot(df_epi_with_vars, inv_mass_model, m_var, "epi_"));
    histograms.emplace_back(
        make_plot(all_nofsi, inv_mass_model, m_var, "noint_"));
  }

  df_epi_with_vars.Snapshot("outtree", output_path + ".tree.root", to_snapshot);

  TFile output_file{output_path.c_str(), "RECREATE"};
  for (auto &hist : histograms) {
    hist->SetDirectory(&output_file);
    hist->Write();
    hist->SetDirectory(nullptr);
  }
  auto dir = output_file.mkdir("per_channel");
  dir->cd();
  for (auto &&[channel, hist] : pi0p_per_channel) {
    hist.Write(normalize_name(channel).c_str());
  }
  output_file.Close();

  std::println("Wrote output to {}", output_path);

  auto vec_count_per_channel =
      count_per_channel |
      std::ranges::to<std::vector<std::pair<std::string, size_t>>>();
  std::ranges::sort(
      vec_count_per_channel,
      [](const auto &a, const auto &b) -> bool { return a.second > b.second; });
  std::println("Channel counts: {}", vec_count_per_channel);

  auto entries_to_plot = vec_count_per_channel | std::views::take(4) |
                         std::ranges::to<std::vector>();
  entries_to_plot.emplace_back(
      "Other",
      event_count.GetValue() -
          std::accumulate(
              entries_to_plot.begin(), entries_to_plot.end(), 0,
              [](size_t sum, const auto &pair) { return sum + pair.second; }));
  std::println("Entries to plot: {}", entries_to_plot);
  make_pie_plot(entries_to_plot, output_path + ".pie.eps");

  // and the cut efficiency
  std::println("Cut efficiency report:");
  cut_efficiency->Print();
  std::println("Cut efficiency report(lower region):");
  cut_efficiency_lower->Print();
  std::println("Cut efficiency report(upper region):");
  cut_efficiency_upper->Print();
}
