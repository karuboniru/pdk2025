#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <TDatabasePDG.h>
#include <TMemFile.h>
#include <TROOT.h>
#include <format>
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
#include "cmdline.h"
#include "common_tools.hxx"
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

int main(int argc, char **argv) {
  initializeGaussianSmearStrategy();

  ROOT::EnableImplicitMT(guess_nproc_from_env());
  TH1::AddDirectory(false);
  auto [input_files, output_path, genie_mode] = parse_command_line(argc, argv);

  auto tracker_df = get_initial_frame(genie_mode, input_files);
  ROOT::RDF::Experimental::AddProgressBar(tracker_df);
  try {
    tracker_df = tracker_df.Define("weight", []() { return 1.0; }, {});
  } catch (...) {
    // weight already exists
  }

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
              {"electron", "pi0_system"});

  auto hist3D =
      df_all.Histo3D({"E_lepton_vs_E_pi0_vs_cos_theta_lepton_pi0",
                      "Lepton Energy vs Pi0 Energy vs Cosine of Angle between "
                      "Lepton and "
                      "Pi0;E_{lepton} [GeV];E_{#pi^{0}} "
                      "[GeV];cos(#theta_{lepton, #pi^{0}})",
                      50, 0., 1.5, 50, 0., 1.5, 50, -1., 1.},
                     "E_lepton", "E_pi0", "cos_theta_lepton_pi0", "weight");

  hist3D->SaveAs(output_path.c_str());
  return 0;
}
