#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <TDatabasePDG.h>
#include <TMemFile.h>
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
#include "decaygen.h"
#include "event.h"
#include "kf.h"
#include "smear.h"
#include "common_tools.hxx"

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

int main(int argc, char **argv) {
  constexpr double to_deg = 180. / M_PI;
  initializeGaussianSmearStrategy();
  ROOT::EnableImplicitMT(guess_nproc_from_env());
  // trigger initialization of everything
  EvtGenInterface::get_instance();
  TDatabasePDG::Instance();
  TH1::AddDirectory(false);
  auto [input_files, output_path, genie_mode] = parse_command_line(argc, argv);

  auto tracker_df =
      genie_mode
          ? TrackerPrepareGENIE(ROOT::RDataFrame{"gRooTracker", input_files})
          : TrackerPrepare(ROOT::RDataFrame{"outtree", input_files});
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

  auto angle_2gamma = [](const momentum_t &a, const momentum_t &b) {
    return std::acos(a.Vect().Unit().Dot(b.Vect().Unit())) * to_deg;
  };

  auto df_epi_final_state =
      df_all
          .Filter(
              [](const size_t nrings) { return nrings == 2 || nrings == 3; },
              {"nrings"}, "2 or 3 rings in detector")
          .Filter(
              [](const size_t nrings, const size_t nshower_rings) {
                return is_mupi ? (nshower_rings + 1 == nrings)
                               : (nshower_rings == nrings);
              },
              {"nrings", "nshower_rings"},
              std::format("{} non-shower-like rings", is_mupi ? "one" : "no"))
          .Filter(
              [](const size_t nmichel_electrons) {
                return nmichel_electrons == (is_mupi ? 1 : 0);
              },
              {"nmichel_electrons"},
              std::format("{} michel electrons", is_mupi ? "one" : "no"))
          .Define("rec",
                  [](const NeutrinoEvent &event) {
                    auto rec = event.Rec_lpi_event(is_mupi);
                    return rec;
                  },
                  {"EventRecord"})
          // .Define("rec_KF",
          //         [](RecResult rec) -> std::optional<RecResult> {
          // if (rec.subleading_gamma.has_value()) {
          //   auto kf_result =
          //       kf_pi0({rec.leading_gamma.m_pair.second,
          //               rec.subleading_gamma->m_pair.second});
          //   if (kf_result.has_value()) {
          //     auto new_pi0 =
          //         kf_result.value()[0] + kf_result.value()[1];
          //     rec.rec_pi0->second = new_pi0;
          //     rec.leading_gamma.m_pair.second =
          //     kf_result.value()[0];
          //     rec.subleading_gamma->m_pair.second =
          //         kf_result.value()[1];
          //     return rec;
          //   }
          //   // std::println("KF failed");
          // }
          // return std::nullopt;
          // },
          // {"rec_raw"})
          // .Define("rec",
          //         [](const std::optional<RecResult> &rec_opt,
          //            const RecResult &rec_raw) {
          //           return rec_opt.value_or(rec_raw);
          //         },
          //         {"rec_KF", "rec_raw"})
          .Define("electron",
                  [](const RecResult &rec) { return rec.lepton.m_pair; },
                  {"rec"})
          .Define("lead_photon",
                  [](const RecResult &rec) { return rec.leading_gamma.m_pair; },
                  {"rec"})
          .Define("sublead_photon",
                  [](const RecResult &rec) {
                    return rec.subleading_gamma.value_or(RingInfo{}).m_pair;
                  },
                  {"rec"})
          .Define("pi0_system",
                  [](const RecResult &rec) {
                    return rec.rec_pi0.value_or(pair_momentum_t{});
                  },
                  {"rec"})
          .Define("epi_system",
                  [](const RecResult &rec) {
                    return rec.lepton.m_pair +
                           rec.rec_pi0.value_or(pair_momentum_t{});
                  },
                  {"rec"});
  //         .Define("epi_system_nokf",
  //                 [](const RecResult &rec) {
  //                   return rec.lepton.m_pair +
  //                          rec.rec_pi0.value_or(pair_momentum_t{});
  //                 },
  //                 {"rec_raw"})
  //         .Define("epi_system_idealpi0",
  //                 [](const RecResult &rec) {
  //                   return rec.lepton.m_pair +
  //                          std::make_pair(rec.rec_pi0->first,
  //                                         rec.rec_pi0->first);
  //                 },
  //                 {"rec_raw"});

  // auto angle_before_fit =
  //     df_epi_final_state
  //         .Filter([](const size_t nrings) { return nrings == 3; },
  //         {"nrings"},
  //                 "3 rings in detector (pre. cond. for KF)")
  //         .Define("angle_2gamma_true",
  //                 [&](const RecResult &rec) {
  //                   if (rec.subleading_gamma.has_value()) {
  //                     return angle_2gamma(rec.leading_gamma.m_pair.first,
  //                                         rec.subleading_gamma->m_pair.first);
  //                   }
  //                   return -1.0;
  //                 },
  //                 {"rec_raw"})
  //         .Define("angle_2gamma_smeared",
  //                 [&](const RecResult &rec) {
  //                   if (rec.subleading_gamma.has_value()) {
  //                     return angle_2gamma(rec.leading_gamma.m_pair.second,
  //                                         rec.subleading_gamma->m_pair.second);
  //                   }
  //                   return -1.0;
  //                 },
  //                 {"rec_raw"});

  // auto angle_after_fit =
  //     angle_before_fit
  //         .Filter(
  //             [](const std::optional<RecResult> &rec_opt) {
  //               return rec_opt.has_value();
  //             },
  //             {"rec_KF"}, "KF succeeded")
  //         .Define("angle_2gamma_smeared_KF",
  //                 [&](const std::optional<RecResult> &rec_opt) {
  //                   const auto &rec = rec_opt.value();
  //                   return angle_2gamma(rec.leading_gamma.m_pair.second,
  //                                       rec.subleading_gamma->m_pair.second);
  //                 },
  //                 {"rec_KF"});

  // auto kf_rate = angle_after_fit.Report();

  ROOT::RDF::TH1DModel angle_model{"angle_2gamma", "angle between 2 gammas",
                                   180, 0.0, 180.};

  auto &&[df_epi_with_vars, to_snapshot, mass_list, p_list] =
      DefineForEPi(ROOT::RDF::RNode{df_epi_final_state});

  std::ranges::copy(std::to_array({"raw_proton_momentum", "raw_mass_proton",
                                   "raw_pi0_before_fsi", "nrings"}),
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

  histograms.emplace_back(
      make_plot(df_all, {"nrings", "nrings", 20, -0.5, 19.5}, "nrings"));

  // histograms.emplace_back(
  //     make_plot(angle_before_fit, angle_model, "angle_2gamma_true", ""));
  // histograms.emplace_back(
  //     make_plot(angle_before_fit, angle_model, "angle_2gamma_smeared", ""));
  // histograms.emplace_back(
  //     make_plot(angle_after_fit, angle_model, "angle_2gamma_smeared_KF",
  //     ""));

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

  auto df_epi_with_vars_3ring =
      df_epi_with_vars.Filter([](const size_t nrings) { return nrings == 3; },
                              {"nrings"}, "3 rings in detector");

  auto smared_pi0_system_m_stddev =
      df_epi_with_vars_3ring.StdDev("smared_pi0_system_m");
  auto smeared_pi0_system_m_mean =
      df_epi_with_vars_3ring.Mean("smared_pi0_system_m");
  auto smeared_epi_system_m_stddev =
      df_epi_with_vars_3ring.StdDev("smared_epi_system_m");
  auto smeared_epi_system_m_mean =
  df_epi_with_vars_3ring.Mean("smared_epi_system_m");

  // count of epi_p < 0.125
  auto low_epi_p_count =
      df_epi_with_vars
          .Filter([](const double epi_p) { return epi_p < 0.125; },
                  {"smared_epi_system_p"}, "epi_p < 0.125")
          .Count();
  // count of epi_p >= 0.125 && epi_p < 0.250
  auto high_epi_p_count =
      df_epi_with_vars
          .Filter(
              [](const double epi_p) {
                return epi_p >= 0.125 && epi_p < 0.250;
              },
              {"smared_epi_system_p"}, "0.125 <= epi_p < 0.250")
          .Count();

  signals[0].Snapshot("outtree", output_path + ".tree.root", to_snapshot);

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

  std::println("KF success rate report:");
  // kf_rate->Print();

  std::println("Smeared pi0 system mass: mean = {}, stddev = {}",
               smeared_pi0_system_m_mean.GetValue(),
               smared_pi0_system_m_stddev.GetValue());
  std::println("Smeared epi system mass: mean = {}, stddev = {}",
               smeared_epi_system_m_mean.GetValue(),
               smeared_epi_system_m_stddev.GetValue());
  std::println("high to low epi_p count ratio: {}",
               static_cast<double>(high_epi_p_count.GetValue()) /
                   low_epi_p_count.GetValue());
}
