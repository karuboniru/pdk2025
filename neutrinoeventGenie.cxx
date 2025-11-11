#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
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

int main(int argc, char **argv) {
  constexpr double to_deg = 180. / M_PI;
  initializeGaussianSmearStrategy();
  // ROOT::EnableImplicitMT(3);
  TH1::AddDirectory(false);
  auto [input_files, output_path, genie_mode] = parse_command_line(argc, argv);

  auto tracker_df =
      TrackerPrepareGENIE(ROOT::RDataFrame{"gRooTracker", input_files});
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

  auto event_count = df_all.Count();

  auto df_epi_final_state =
      df_all
          .Filter(
              [](const size_t nrings) { return nrings == 2 || nrings == 3; },
              {"nrings"}, "2 or 3 rings in detector")
          .Filter(
              [](const size_t nrings, const size_t nshower_rings) {
                return nshower_rings == nrings;
              },
              {"nrings", "nshower_rings"}, "all shower-like rings")
          .Filter(
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
      DefineForEPi(ROOT::RDF::RNode{df_epi_final_state});

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

  signals[0].Snapshot("outtree", output_path + ".tree.root", to_snapshot);

  TFile output_file{output_path.c_str(), "RECREATE"};
  for (auto &hist : histograms) {
    hist->SetDirectory(&output_file);
    hist->Write();
    hist->SetDirectory(nullptr);
  }

  std::println("Wrote output to {}", output_path);

  // and the cut efficiency
  std::println("Cut efficiency report:");
  cut_efficiency->Print();
  std::println("Cut efficiency report(lower region):");
  cut_efficiency_lower->Print();
  std::println("Cut efficiency report(upper region):");
  cut_efficiency_upper->Print();
}
