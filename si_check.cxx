#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TMemFile.h>
#include <TROOT.h>
#include <format>
#include <memory>
#include <print>
#include <ranges>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPie.h>
#include <TSystem.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "EvtTracker2event.h"
#include "ROOT/RResultPtr.hxx"
#include "cmdline.h"
#include "common_tools.hxx"
#include "commondefine.h"
#include "event.h"
#include "smear.h"

double cos_theta_between_vectors(const ROOT::Math::PxPyPzEVector &v1,
                                 const ROOT::Math::PxPyPzEVector &v2) {
  double dot_product = v1.Vect().Dot(v2.Vect());
  double magnitude_product = v1.P() * v2.P();
  if (magnitude_product == 0) {
    return 0.0; // Avoid division by zero
  }
  double cos_theta = dot_product / magnitude_product;
  // Clamp the value to the valid range for acos to avoid NaN due to numerical
  // errors
  cos_theta = std::clamp(cos_theta, -1.0, 1.0);
  return cos_theta;
}

ROOT::RDF::RNode get_initial_frame(bool genie_mode,
                                   const std::vector<std::string> &filenames) {
  if (genie_mode) {
    return TrackerPrepareGENIE(ROOT::RDataFrame{"gRooTracker", filenames});
  }
  try {
    return TrackerPrepare(ROOT::RDataFrame{"outtree", filenames});
  } catch (...) {
    // some historical reasons...
    return TrackerPrepareNeutrino(ROOT::RDataFrame{"out_tree", filenames});
  }
}

std::tuple<ROOT::RDF::RNode, std::unique_ptr<TChain>, std::unique_ptr<TChain>>
get_initial_frame_corr(bool genie_mode,
                       const std::vector<std::string> &filenames,
                       const std::vector<std::string> &corr_files) {
  if (corr_files.empty()) {
    return {get_initial_frame(genie_mode, filenames), nullptr, nullptr};
  }

  auto main_tree = genie_mode ? "gRooTracker" : "out_tree";
  auto corr_tree = "out_tree";
  auto chain_main = std::make_unique<TChain>(main_tree);
  for (const auto &file : filenames) {
    chain_main->Add(file.c_str());
  }
  auto chain_corr = std::make_unique<TChain>(corr_tree);
  for (const auto &file : corr_files) {
    chain_corr->Add(file.c_str());
  }
  chain_main->AddFriend(chain_corr.get(), "corr");
  auto df = ROOT::RDataFrame(*chain_main);
  if (genie_mode) {
    return {TrackerPrepareGENIE(df), std::move(chain_main),
            std::move(chain_corr)};
  }
  return {TrackerPrepareNeutrino(df), std::move(chain_main),
          std::move(chain_corr)};
}

std::vector<double> log_bin_edges(int count, double max) {
  std::vector<double> edges(count + 1);
  double log_min = std::log10(1e-3);
  double log_max = std::log10(max);
  double step = (log_max - log_min) / count;
  for (int i = 0; i <= count; ++i) {
    edges[i] = std::pow(10, log_min + i * step);
  }
  return edges;
}

std::vector<ROOT::RDF::RResultPtr<TH1>>
plot_from_df_impl(ROOT::RDF::RNode node, const std::string &name,
                  const std::string &weight_column,
                  const std::string &suffix_plot = "norm") {
  auto edges_5 = log_bin_edges(50, 5.0);
  auto edges_15 = log_bin_edges(50, 15.0);
  return {
      node.Histo1D({std::format("{}_{}_{}", name, "Q2", suffix_plot).c_str(),
                    ";Q^{2};a.u.", 50, 0.0, 15.0},
                   "Q2", weight_column),
      node.Histo1D({std::format("{}_{}_{}", name, "W", suffix_plot).c_str(),
                    ";W;a.u.", 30, 0.0, 5.0},
                   "W", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}_log", name, "Q2", suffix_plot).c_str(),
           ";Q^{2};a.u.", 50, edges_15.data()},
          "Q2", weight_column),
      node.Histo1D({std::format("{}_{}_{}_log", name, "W", suffix_plot).c_str(),
                    ";W;a.u.", 50, edges_5.data()},
                   "W", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}", name, "init_p", suffix_plot).c_str(),
           ";p_{N};a.u.", 50, 0.0, 1.0},
          "initial_proton_p", weight_column),
      node.Histo1D({std::format("{}_{}_{}", name, "npM", suffix_plot).c_str(),
                    ";M_{#nu N};a.u.", 100, 0.0, 6.0},
                   "np_system_m", weight_column),
      node.Histo1D({std::format("{}_{}_{}", name, "npP", suffix_plot).c_str(),
                    ";M_{#nu N};a.u.", 100, 0.0, 8.0},
                   "np_system_p", weight_column),
      node.Histo2D(
          {std::format("{}_{}_{}", name, "hist2d_initp_npP", suffix_plot)
               .c_str(),
           ";p_{N};M_{#nu N};a.u.", 40, 0.0, 1.0, 40, 0.0, 8.0},
          "initial_proton_p", "np_system_p", weight_column),

      //
      // before SI
      node.Histo1D(
          {std::format("{}_{}_{}", name, "total_p_before_SI", suffix_plot)
               .c_str(),
           ";p_{tot}^{before SI};a.u.", 100, 0.0, 8.0},
          "total_p_before_SI", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}", name, "total_m_before_SI", suffix_plot)
               .c_str(),
           ";M_{tot}^{before SI};a.u.", 100, 0.0, 6.0},
          "total_m_before_SI", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}", name, "W_before_SI", suffix_plot).c_str(),
           ";W^{before SI};a.u.", 100, 0.0, 6.0},
          "W_before_SI", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}", name, "total_p_after_SI", suffix_plot)
               .c_str(),
           ";p_{tot}^{after SI};a.u.", 100, 0.0, 8.0},
          "total_p_after_SI", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}", name, "total_m_after_SI", suffix_plot)
               .c_str(),
           ";M_{tot}^{after SI};a.u.", 100, 0.0, 6.0},
          "total_m_after_SI", weight_column),
      node.Histo1D(
          {std::format("{}_{}_{}", name, "W_after_SI", suffix_plot).c_str(),
           ";W^{after SI};a.u.", 100, 0.0, 6.0},
          "W_after_SI", weight_column),

  };
}

std::vector<ROOT::RDF::RResultPtr<TH1>>
plot_from_df(const ROOT::RDF::RNode &node, const std::string &name) {
  std::vector<ROOT::RDF::RResultPtr<TH1>> ret;
  for (auto &&obj : plot_from_df_impl(node, name, "weight")) {
    ret.emplace_back(obj);
  }
  return ret;
}

auto calc_tot_p4(const ROOT::RVec<int> &StdHepStatus,
                 const ROOT::RVecD &StdHepP4) {
  ROOT::Math::PxPyPzEVector total_momentum{};
  for (auto &&[status, p4] :
       std::views::zip(StdHepStatus,
                       std::ranges::iota_view(StdHepP4.data()) |
                           std::views::chunk(4) |
                           std::views::transform([](auto &&chunk) {
                             return ROOT::Math::PxPyPzEVector{
                                 *chunk[0], *chunk[1], *chunk[2], *chunk[3]};
                           })) |
           std::views::filter(
               [](auto &&entry) { return std::get<0>(entry) == 1; })) {
    total_momentum += p4;
  }
  return total_momentum;
}

auto calc_tot_p4_had_system(const ROOT::RVec<int> &StdHepStatus,
                            const ROOT::RVec<int> &StdHepPdg,
                            const ROOT::RVecD &StdHepP4) {
  ROOT::Math::PxPyPzEVector total_momentum{};
  for (auto &&[status, pdg, p4] :
       std::views::zip(StdHepStatus, StdHepPdg,
                       std::ranges::iota_view(StdHepP4.data()) |
                           std::views::chunk(4) |
                           std::views::transform([](auto &&chunk) {
                             return ROOT::Math::PxPyPzEVector{
                                 *chunk[0], *chunk[1], *chunk[2], *chunk[3]};
                           })) |
           std::views::filter([](auto &&entry) {
             auto &[status, pdg, p4] = entry;
             if (status != 1)
               return false;
             int abs_pdg = std::abs(pdg);
             return abs_pdg != 11 && abs_pdg != 13 && abs_pdg != 15;
           })) {
    total_momentum += p4;
  }
  return total_momentum;
}

int main(int argc, char **argv) {
  gSystem->ResetSignal(kSigBus);
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  initializeGaussianSmearStrategy();
  ROOT::EnableImplicitMT(guess_nproc_from_env());
  TH1::AddDirectory(false);
  // auto [input_files, input_corr, output_path, genie_mode] =
  //     parse_command_line(argc, argv)
  auto &&cfg = parse_command_line(argc, argv);
  std::println("{}", cfg);
  auto &&[input_files, input_corr, output_path, genie_mode] = cfg;

  // auto tracker_df = get_initial_frame(genie_mode, input_files);
  auto &&[tracker_df, chain_main, chain_corr] =
      get_initial_frame_corr(genie_mode, input_files, input_corr);
  try {
    tracker_df = tracker_df.Define("weight", []() { return 1.0; }, {});
  } catch (...) {
    // weight already exists
  }

  if (chain_corr) {
    tracker_df =
        tracker_df
            .Define("initial_proton_p4",
                    [](const ROOT::RVecD &arr) {
                      return ROOT::Math::PxPyPzEVector(arr[4], arr[5], arr[6],
                                                       arr[7]);
                    },
                    {"corr.StdHepP4"})
            .Define("initial_proton_p",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.P(); },
                    {"initial_proton_p4"})
            .Define("initial_proton_m",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.M(); },
                    {"initial_proton_p4"})
            .Define("init_neutrino",
                    [](const ROOT::RVecD &arr) {
                      return ROOT::Math::PxPyPzEVector(arr[0], arr[1], arr[2],
                                                       arr[3]);
                    },
                    {"corr.StdHepP4"})
            .Define("np_system_p4",
                    [](const ROOT::Math::PxPyPzEVector &init_nu,
                       const ROOT::Math::PxPyPzEVector &init_p) {
                      return init_nu + init_p;
                    },
                    {"init_neutrino", "initial_proton_p4"})
            .Define("np_system_p",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.P(); },
                    {"np_system_p4"})
            .Define("np_system_m",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.M(); },
                    {"np_system_p4"})
            .Define("final_state_lepton_p4",
                    [](const ROOT::RVecD &arr) {
                      return ROOT::Math::PxPyPzEVector(arr[8], arr[9], arr[10],
                                                       arr[11]);
                    },
                    {"corr.StdHepP4"})
            .Define("Q2",
                    [](const ROOT::Math::PxPyPzEVector &init_nu,
                       const ROOT::Math::PxPyPzEVector &final_lepton) {
                      return -1.0 * (init_nu - final_lepton).M2();
                    },
                    {"init_neutrino", "final_state_lepton_p4"})
            .Define(
                "W",
                [](const ROOT::Math::PxPyPzEVector &init_nu,
                   const ROOT::Math::PxPyPzEVector &init_p,
                   const ROOT::Math::PxPyPzEVector &final_lepton) {
                  return (init_nu - final_lepton + init_p).M();
                },
                {"init_neutrino", "initial_proton_p4", "final_state_lepton_p4"})
            .Define(
                "scaled_w",
                [](double weight, int n_capture) { return weight * n_capture; },
                {"weight", "n_capture"})
            .Define("total_p4_before_SI", calc_tot_p4,
                    {"corr.StdHepStatus", "corr.StdHepP4"})
            .Define("total_p_before_SI",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.P(); },
                    {"total_p4_before_SI"})
            .Define("total_m_before_SI",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.M(); },
                    {"total_p4_before_SI"})
            .Define("total_p4_had_system_before_SI", calc_tot_p4_had_system,
                    {"corr.StdHepStatus", "corr.StdHepPdg", "corr.StdHepP4"})
            .Define("W_before_SI",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.M(); },
                    {"total_p4_before_SI"})
            .Define("total_p4_after_SI", calc_tot_p4,
                    {"StdHepStatus", "StdHepP4"})
            .Define("total_p_after_SI",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.P(); },
                    {"total_p4_after_SI"})
            .Define("total_m_after_SI",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.M(); },
                    {"total_p4_after_SI"})
            .Define("total_p4_had_system_after_SI", calc_tot_p4_had_system,
                    {"StdHepStatus", "StdHepPdg", "StdHepP4"})
            .Define("W_after_SI",
                    [](const ROOT::Math::PxPyPzEVector &p4) { return p4.M(); },
                    {"total_p4_after_SI"});
  }

  std::vector<ROOT::RDF::RResultPtr<TH1>> histograms;

  auto add_plots = [&](const std::vector<ROOT::RDF::RResultPtr<TH1>> &plots) {
    for (const auto &plot : plots) {
      histograms.emplace_back(plot);
    }
  };
  auto total_weight = tracker_df.Sum("weight");
  // histograms.emplace_back(plot_from_df(tracker_df, "all"));
  add_plots(plot_from_df(tracker_df, "all"));

  TFile output_file(output_path.c_str(), "RECREATE");
  for (auto &&hist : histograms) {
    auto hist_cloned = dynamic_cast<TH1 *>(hist->Clone());
    hist_cloned->SetDirectory(&output_file);
  }
  output_file.Write();
  output_file.Close();

  return 0;
}
