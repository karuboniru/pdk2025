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
  return df.Histo1D(model, varname, "weight");
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
  ROOT::EnableImplicitMT(4);
  TH1::AddDirectory(false);
  auto [input_files, output_path] = parse_command_line(argc, argv);
  auto nfile = input_files.size();

  auto tracker_df =
      TrackerPrepareNeutrino(ROOT::RDataFrame{"out_tree", input_files});
  ROOT::RDF::Experimental::AddProgressBar(tracker_df);
  auto df_all = tracker_df
                    .Define("channel_name",
                            [](NeutrinoEvent &e) {
                              return e.get_channelname_no_nucleon();
                            },
                            {"EventRecord"})
                    .Define("NeutrinoEnergy",
                            [](const ROOT::RVec<double> &StdHepP4_) {
                              return StdHepP4_[3];
                            },
                            {"StdHepP4"});
  auto count_per_channel =
      df_all
          .Define("channel_name_weight",
                  [](const std::string &name,
                     double weight) -> std::pair<std::string, double> {
                    return std::make_pair(name, weight);
                  },
                  {"channel_name", "weight"})
          .Aggregate(
              [](std::map<std::string, double> &data,
                 const std::pair<std::string, double> &col) {
                auto &[name, weight] = col;
                data[name] += weight;
              },
              [](std::vector<std::map<std::string, double>> &to_merge) {
                for (auto &target = to_merge[0];
                     const auto &item : to_merge | std::views::drop(1)) {
                  for (const auto &[key, value] : item) {
                    target[key] += value;
                  }
                }
              },
              "channel_name_weight", std::map<std::string, double>{});

  // auto event_count = df_all.Count();
  auto weight_sum = df_all.Sum("weight");

  auto df_epi_final_state =
      FilterTrackedRDF{df_all}
          .SetWeightColumnName("weight")
          .FilterTracked(
              [](const size_t nrings) { return nrings == 2 || nrings == 3; },
              {"nrings"}, "2 or 3 rings in detector")
          .FilterTracked(
              [](const size_t nrings, const size_t nshower_rings) {
                return nshower_rings == nrings;
              },
              {"nrings", "nshower_rings"}, "all shower-like rings")
          .FilterTracked(
              [](const size_t nmichel_electrons) {
                return nmichel_electrons == 0;
              },
              {"nmichel_electrons"}, "no michel electrons")
          .Define(
              "rec",
              [](const NeutrinoEvent &event) { return event.Rec_lpi_event(); },
              {"EventRecord"})
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
                  [](const pair_momentum_t &electron,
                     const pair_momentum_t &pi0_system) {
                    return electron + pi0_system;
                  },
                  {"electron", "pi0_system"});
  auto &&[df_epi_with_vars, to_snapshot, mass_list, p_list] =
      DefineForEPi(df_epi_final_state);
  to_snapshot.push_back("weight");
  to_snapshot.push_back("NeutrinoEnergy");
  to_snapshot.push_back("channel");

  auto signals = FilterSignalKinematics(df_epi_with_vars);
  auto weight_sum_signal = signals |
                           std::views::transform([](ROOT::RDF::RNode &node) {
                             return node.Sum("weight");
                           }) |
                           std::ranges::to<std::vector>();

  ROOT::RDF::TH1DModel inv_mass_model{"inv_mass_epip_system", "inv mass", 400,
                                      0.0, 1.0};
  ROOT::RDF::TH1DModel momentum_model{"inv_mass_epip_system", "momentum", 400,
                                      0.0, 1.0};
  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms{};

  histograms.reserve(p_list.size() + mass_list.size());
  for (auto &p_var : p_list) {
    histograms.emplace_back(
        make_plot(df_epi_with_vars, momentum_model, p_var, "epi_"));
  }

  for (auto &m_var : mass_list) {
    histograms.emplace_back(
        make_plot(df_epi_with_vars, inv_mass_model, m_var, "epi_"));
  }

  signals[0].Snapshot("outtree", output_path + ".tree.root", to_snapshot);

  TFile output_file{output_path.c_str(), "RECREATE"};
  for (auto &hist : histograms) {
    hist->SetDirectory(&output_file);
    hist->Write();
    hist->SetDirectory(nullptr);
  }
  output_file.Close();

  std::println("Wrote output to {}", output_path);

  auto vec_count_per_channel =
      count_per_channel |
      std::ranges::to<std::vector<std::pair<std::string, double>>>();
  std::ranges::sort(
      vec_count_per_channel,
      [](const auto &a, const auto &b) -> bool { return a.second > b.second; });
  std::println("Channel counts: {}", vec_count_per_channel);

  auto entries_to_plot = vec_count_per_channel | std::views::take(4) |
                         std::ranges::to<std::vector>();
  entries_to_plot.emplace_back(
      "Other",
      weight_sum.GetValue() - std::accumulate(entries_to_plot.begin(),
                                              entries_to_plot.end(), 0.,
                                              [](size_t sum, const auto &pair) {
                                                return sum + pair.second;
                                              }));
  std::println("Entries to plot: {}", entries_to_plot);
  make_pie_plot(entries_to_plot, output_path + ".pie.eps");

  auto weight_ratio_signal = weight_sum_signal |
                             std::views::transform([&weight_sum](auto w) {
                               return w.GetValue() / weight_sum.GetValue();
                             }) |
                             std::ranges::to<std::vector>();
  std::println("Signal weight {}, {}, {}", weight_sum_signal[0].GetValue(),
               weight_sum_signal[1].GetValue(),
               weight_sum_signal[2].GetValue());
  std::println("Signal weight ratios: {}", weight_ratio_signal);
  std::println("total weight/nfile = {} / {} = {}", weight_sum.GetValue(),
               nfile, weight_sum.GetValue() / nfile);

  for (auto &signal : signals) {
    signal.Report();
    std::println("");
  }
}
