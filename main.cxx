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
  ROOT::EnableImplicitMT();
  auto [input_files, output_path] = parse_command_line(argc, argv);

  auto tracker_df = TrackerPrepare(ROOT::RDataFrame{"outtree", input_files});

  auto df_all =
      tracker_df
          .Define("raw_mass_proton",
                  [](const NeutrinoEvent &event) {
                    auto proton = event.in_range(2212);
                    return proton.begin()->second.M();
                  },
                  {"EventRecord"})
          .Define("final_state_mass",
                  [](const NeutrinoEvent &event) {
                    const auto &post = event.get_post();
                    TLorentzVector total_momentum;
                    for (const auto &p :
                         post | std::views::filter([](const auto &entry) {
                           return entry.first != 2212 && entry.first != 2112;
                         }) | std::views::values) {
                      total_momentum += p;
                    }
                    return total_momentum.M();
                  },
                  {"EventRecord"})
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
                return event.count_post(-11) != 0 && event.count_post(111) != 0;
              },
              {"EventRecord"})
          .Define("final_state_mass_epip_system",
                  [](const NeutrinoEvent &event) {
                    const auto &post = event.get_post();
                    TLorentzVector total_momentum;
                    for (const auto &p :
                         post | std::views::filter([](const auto &entry) {
                           auto pdg = entry.first;
                           return pdg == -11 || pdg == 111;
                         }) | std::views::values) {
                    }
                    return total_momentum.M();
                  },
                  {"EventRecord"})
          .Define("final_state_mass_leading_epip_system",
                  [](const NeutrinoEvent &event) {
                    auto total_momentum =
                        event.get_leading(-11) + event.get_leading(111);
                    return total_momentum.M();
                  },
                  {"EventRecord"});

  ROOT::RDF::TH1DModel inv_mass_model{
      "inv_mass_epip_system", "Invariant mass of e+pi- system", 100, 0.3, 1.0};

  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms{};
  histograms.emplace_back(make_plot(df_all, inv_mass_model, "raw_mass_proton"));
  histograms.emplace_back(
      make_plot(df_all, inv_mass_model, "final_state_mass"));

  auto final_hists = std::to_array({"final_state_mass", "raw_mass_proton",
                                    "final_state_mass_epip_system",
                                    "final_state_mass_leading_epip_system"});
  for (const auto &varname : final_hists) {
    histograms.emplace_back(
        make_plot(all_with_vars, inv_mass_model, varname, "epi_"));
  }

  // df_all.Snapshot("outtree", output_path + ".tree.root",
  //                 {"raw_mass_proton", "final_state_mass", "channel_name"});

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
  make_pie_plot(entries_to_plot, output_path + ".pie.pdf");
}
