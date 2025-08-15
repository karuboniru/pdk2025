#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <TROOT.h>
#include <array>
#include <format>
#include <iostream>
#include <print>
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
#include "event.h"
#include "smear.h"

ROOT::RDF::RResultPtr<TH1D> make_plot(auto df, ROOT::RDF::TH1DModel model,
                                      const std::string &varname,
                                      const std::string &prefix = "") {
  model.fName = prefix + varname;
  model.fTitle = prefix + varname;
  return df.Histo1D(model, varname);
}

void make_pie_plot(auto &data, std::string filename) {
  // constexpr std::array<int, 10> col{kRed,   kBlue, kViolet, kYellow, kOrange,
  //                                   kGreen, kGray, kTeal,   kPink};
  auto col =
      std::to_array({kP6Blue, kP6Yellow, kP6Red, kP6Grape, kP6Gray, kP6Violet});
  auto pie = std::make_unique<TPie>("final state", "final state", data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    pie->SetEntryVal(i, data[i].second);
    pie->SetEntryFillColor(i, col[i % col.size()]);
    pie->SetEntryFillStyle(i, 1000 + i / col.size());
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
  initializeGaussianSmearStrategy();
  ROOT::EnableImplicitMT(4);
  auto [input_files, output_path] = parse_command_line(argc, argv);

  auto tracker_df = TrackerPrepare(ROOT::RDataFrame{"outtree", input_files});

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
                    // if (pi0.size() == 0) {
                    //   return ROOT::Math::PxPyPzEVector{0, 0, 0, 0};
                    // }
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
      {"channel_name"}, std::map<std::string, std::size_t>{});
  auto event_count = df_all.Count();

  auto all_with_vars =
      df_all
          .Filter(
              [](const NeutrinoEvent &event) {
                return event.count_det(-11) == 1 && event.count_det(22) == 2 &&
                       event.count_det(211) == 0 && event.count_det(-211) == 0;
              },
              {"EventRecord"})
          .Define("raw_electron_momentum",
                  [](const NeutrinoEvent &event) {
                    auto electron = event.post_range(-11);
                    return electron.begin()->second.P();
                  },
                  {"EventRecord"})
          .Define("smear_electron_momentum",
                  [](const NeutrinoEvent &event) {
                    auto electron = event.det_range(-11);
                    return electron.begin()->second.P();
                  },
                  {"EventRecord"})
          .Define("angle_2gamma_no_smear",
                  [](const NeutrinoEvent &event) {
                    const auto &detector = event.get_before_smear();
                    if (detector.count(22) < 2) {
                      return -1.0; // Not enough photons
                    }
                    auto it = detector.equal_range(22);
                    auto p4gamma1 = it.first->second;
                    auto p4gamma2 = (++it.first)->second;
                    auto costheta =
                        p4gamma1.Vect().Unit().Dot(p4gamma2.Vect().Unit());
                    return std::acos(costheta) * 180.0 /
                           M_PI; // Convert to degrees
                  },
                  {"EventRecord"})
          .Define("angle_2gamma_after_smear",
                  [](const NeutrinoEvent &event) {
                    const auto &detector = event.get_det();
                    if (detector.count(22) < 2) {
                      return -1.0; // Not enough photons
                    }
                    auto it = detector.equal_range(22);
                    auto p4gamma1 = it.first->second;
                    auto p4gamma2 = (++it.first)->second;
                    auto costheta =
                        p4gamma1.Vect().Unit().Dot(p4gamma2.Vect().Unit());
                    return std::acos(costheta) * 180.0 /
                           M_PI; // Convert to degrees
                  },
                  {"EventRecord"})
          .Define("raw_pi0_mass",
                  [](const NeutrinoEvent &event) {
                    auto &p4pi0 = event.get_leading(111);
                    return p4pi0.M();
                  },
                  {"EventRecord"})
          .Define("angle_epi0_truth",
                  [](const NeutrinoEvent &event) {
                    const auto &p4e = event.get_leading(-11);
                    const auto &p4pi0 = event.get_leading(111);
                    auto costheta = p4e.Vect().Unit().Dot(p4pi0.Vect().Unit());
                    return std::acos(costheta) * 180.0 / M_PI;
                  },
                  {"EventRecord"})
          .Define("raw_pi0_p",
                  [](const NeutrinoEvent &event) {
                    auto &p4pi0 = event.get_leading(111);
                    return p4pi0.P();
                  },
                  {"EventRecord"})
          .Define("rec_pi0",
                  [](const NeutrinoEvent &event) {
                    const auto &detector = event.get_det();
                    ROOT::Math::PxPyPzEVector total_momentum;
                    for (const auto &p :
                         detector | std::views::filter([](const auto &entry) {
                           auto pdg = entry.first;
                           return pdg == 22;
                         }) | std::views::values) {
                      total_momentum += p;
                    }
                    return total_momentum;
                  },
                  {"EventRecord"})
          .Define(
              "rec_pi0_M",
              [](const ROOT::Math::PxPyPzEVector &p4pi0) { return p4pi0.M(); },
              {"rec_pi0"})
          .Define(
              "rec_pi0_p",
              [](const ROOT::Math::PxPyPzEVector &p4pi0) { return p4pi0.P(); },
              {"rec_pi0"})
          .Define("angle_epi0_rec",
                  [](const NeutrinoEvent &event,
                     const ROOT::Math::PxPyPzEVector &p4pi0) {
                    const auto &p4e = event.get_leading_det(-11);
                    auto costheta = p4e.Vect().Unit().Dot(p4pi0.Vect().Unit());
                    return std::acos(costheta) * 180.0 / M_PI;
                  },
                  {"EventRecord", "rec_pi0"})
          .Define("rec_epi_system",
                  [](const NeutrinoEvent &event,
                     const ROOT::Math::PxPyPzEVector &p4pi0) {
                    const auto &p4e = event.get_leading_det(-11);
                    return (p4e + p4pi0);
                  },
                  {"EventRecord", "rec_pi0"})
          .Define("rec_mass_epi_system",
                  [](const ROOT::Math::PxPyPzEVector &epi_system) {
                    return epi_system.M();
                  },
                  {"rec_epi_system"})
          .Define("rec_p_epi_system",
                  [](const ROOT::Math::PxPyPzEVector &epi_system) {
                    return epi_system.P();
                  },
                  {"rec_epi_system"});

  ROOT::RDF::TH1DModel inv_mass_model{"inv_mass_epip_system", "inv mass", 100,
                                      0.3, 1.0};
  ROOT::RDF::TH1DModel momentum_model{"inv_mass_epip_system", "momentum", 100,
                                      0.0, 1.0};
  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms{};

  for (const auto &varname :
       std::to_array({"raw_mass_proton", "raw_final_state_mass"})) {
    histograms.emplace_back(make_plot(df_all, inv_mass_model, varname));
    histograms.emplace_back(
        make_plot(all_with_vars, inv_mass_model, varname, "epi_"));
  }

  for (const auto &varname :
       std::to_array({"raw_final_state_p", "raw_proton_momentum",
                      "raw_pi0_before_fsi"})) {
    histograms.emplace_back(make_plot(df_all, momentum_model, varname));
    histograms.emplace_back(
        make_plot(all_with_vars, momentum_model, varname, "epi_"));
  }

  histograms.emplace_back(make_plot(
      all_with_vars, {"inv_mass_epip_system", "momentum", 100, 0.0, 0},
      "rec_pi0_M", "epi_"));
  histograms.emplace_back(
      make_plot(all_with_vars, inv_mass_model, "rec_mass_epi_system", "epi_"));

  for (const auto &varname : {"rec_p_epi_system", "rec_pi0_p", "raw_pi0_p"}) {
    histograms.emplace_back(
        make_plot(all_with_vars, momentum_model, varname, "epi_"));
  }
  all_with_vars.Snapshot(
      "outtree", output_path + ".tree.root",
      {"raw_proton_momentum", "raw_mass_proton", "rec_pi0_M", "rec_pi0_p",
       "raw_pi0_p", "rec_mass_epi_system", "rec_p_epi_system", "raw_pi0_mass",
       "raw_final_state_mass", "angle_epi0_rec", "angle_epi0_truth",
       "angle_2gamma_after_smear", "angle_2gamma_no_smear",
       "raw_electron_momentum", "smear_electron_momentum"});

  TFile output_file{output_path.c_str(), "RECREATE"};
  for (auto &hist : histograms) {
    hist->SetDirectory(&output_file);
    hist->Write();
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
}
