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

std::vector<ROOT::RDF::RResultPtr<TH1>> plot_from_df(ROOT::RDF::RNode node,
                                                     const std::string &name) {
  return {node.Histo1D({name.c_str(), ";p_{N};a.u.", 50, 0.0, 1.0},
                       "initial_proton_p", "weight"),
          node.Histo1D({(name + "M").c_str(), ";M_{#nu N};a.u.", 100, 0.0, 6.0},
                       "np_system_m", "weight"),
          node.Histo1D({(name + "P").c_str(), ";M_{#nu N};a.u.", 100, 0.0, 8.0},
                       "np_system_p", "weight"),
          node.Histo2D({std::format("{}_{}", name, "hist2d_initp_npP").c_str(),
                        ";p_{N};M_{#nu N};a.u.", 40, 0.0, 1.0, 40, 0.0, 8.0},
                       "initial_proton_p", "np_system_p", "weight")};
}

int main(int argc, char **argv) {
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
                    {"np_system_p4"});
    ;
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
  auto df_sliced_1 =
      tracker_df
          .Define("raw_proton_momentum",
                  [](const NeutrinoEvent &event) {
                    auto proton = event.in_range(2212);
                    if (proton.size() == 0) {
                      return 0.0;
                    }
                    return proton.begin()->second.P();
                  },
                  {"EventRecord"})
          .Define("raw_mass_proton",
                  [](const NeutrinoEvent &event) {
                    auto proton = event.in_range(2212);
                    if (proton.size() == 0) {
                      return 0.0;
                    }
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
                    if (pi0.size() == 0) {
                      return 0.0;
                    }
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
          .Filter(
              [](const size_t nrings) { return nrings == 2 || nrings == 3; },
              {"nrings"}, "2 or 3 rings in detector");
  add_plots(plot_from_df(df_sliced_1, "ring_cut"));
  auto df_sliced2 = df_sliced_1.Filter(
      [](const size_t nrings, const size_t nshower_rings) {
        return is_mupi ? (nshower_rings + 1 == nrings)
                       : (nshower_rings == nrings);
      },
      {"nrings", "nshower_rings"},
      std::format("{} non-shower-like rings", is_mupi ? "one" : "no"));
  add_plots(plot_from_df(df_sliced2, "shower_like_cut"));
  auto df_sliced_3 = df_sliced2.Filter(
      [](const size_t nmichel_electrons) {
        return nmichel_electrons == (is_mupi ? 1 : 0);
      },
      {"nmichel_electrons"},
      std::format("{} michel electrons", is_mupi ? "one" : "no"));
  add_plots(plot_from_df(df_sliced_3, "michel_electron_cut"));
  auto df_sliced_4 =
      df_sliced_3
          .Define("rec",
                  [](const NeutrinoEvent &event) {
                    auto rec = event.Rec_lpi_event(is_mupi);
                    return rec;
                  },
                  {"EventRecord"})
          .Define("electron",
                  [](const RecResult &rec) { return rec.lepton.m_pair; },
                  {"rec"})
          .Define("pi0_system",
                  [](const RecResult &rec) {
                    return rec.rec_pi0.value_or(rec.leading_gamma.m_pair);
                  },
                  {"rec"})
          .Define(
              "E_lepton",
              [](const momentum_pair &momentum) { return momentum.second.E(); },
              {"electron"})
          .Define(
              "E_pi0",
              [](const momentum_pair &momentum) { return momentum.second.E(); },
              {"pi0_system"})
          .Define(
              "cos_theta_lepton_pi0",
              [](const momentum_pair &lepton, const momentum_pair &pi0_system) {
                return cos_theta_between_vectors(lepton.second,
                                                 pi0_system.second);
              },
              {"electron", "pi0_system"})
          .Define(
              "is_transparent",
              [](const NeutrinoEvent &event) { return event.is_transparent(); },
              {"EventRecord"})
          .Define(
              "pi0_mass",
              [](const momentum_pair &momentum) { return momentum.second.M(); },
              {"pi0_system"})
          .Filter(
              [](const double pi0_mass, const size_t nrings) {
                return nrings == 2 || (pi0_mass > 0.085 && pi0_mass < 0.185);
              },
              {"pi0_mass", "nrings"}, "Pi0 mass between 85 MeV and 185 MeV");
  add_plots(plot_from_df(df_sliced_4, "pi0_mass_cut"));
  auto df_sliced_5 =
      df_sliced_4
          .Define("rec_proton_system",
                  [](const RecResult &rec) {
                    return rec.lepton.m_pair.second +
                           rec.leading_gamma.m_pair.second +
                           rec.subleading_gamma.value_or({}).m_pair.second;
                  },
                  {"rec"})
          .Define("smear_proton_momentum",
                  [](const ROOT::Math::PxPyPzEVector &proton_system) {
                    return proton_system.P();
                  },
                  {"rec_proton_system"})
          .Define("smear_proton_mass",
                  [](const ROOT::Math::PxPyPzEVector &proton_system) {
                    return proton_system.M();
                  },
                  {"rec_proton_system"})
          .Filter([](double m) { return m > 0.8 && m < 1.05; },
                  {"smear_proton_mass"},
                  "Reconstructed proton mass between 800 MeV and 1050 MeV");
  add_plots(plot_from_df(df_sliced_5, "proton_mass_cut"));
  auto df_sliced_6 = df_sliced_5.Filter([](double p) { return p < 0.25; },
                                        {"smear_proton_momentum"},
                                        "Proton momentum less than 250 MeV");

  add_plots(plot_from_df(df_sliced_6, "proton_momentum_cut"));

  auto cut_report = df_sliced_6.Report();

  TFile output_file(output_path.c_str(), "RECREATE");
  for (auto &&hist : histograms) {
    auto hist_cloned = dynamic_cast<TH1 *>(hist->Clone());
    hist_cloned->SetDirectory(&output_file);
  }
  output_file.Write();
  output_file.Close();

  cut_report->Print();

  return 0;
}
