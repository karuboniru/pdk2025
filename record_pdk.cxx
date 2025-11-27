#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TDatabasePDG.h>
#include <TMemFile.h>
#include <TROOT.h>
#include <array>
#include <optional>
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
#include "cmdline.h"
#include "data.h"
#include "event.h"
#include "kf.h"
#include "smear.h"

int main(int argc, char **argv) {
  constexpr double to_deg = 180. / M_PI;
  initializeGaussianSmearStrategy();
  // ROOT::TTreeProcessorMT::SetTasksPerWorkerHint(-1);
  ROOT::EnableImplicitMT();
  // trigger initialization of everything
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
          .Define(
              "channel_name",
              [](NeutrinoEvent &e) { return e.get_channelname_no_nucleon(); },
              {"EventRecord"})
          .Define(
              "nrings_cut",
              [](const size_t nrings) { return nrings == 2 || nrings == 3; },
              {"nrings"})
          .Define("shower_ring_cut",
                  [](const size_t nrings, const size_t nshower_rings) {
                    return is_mupi ? (nshower_rings + 1 == nrings)
                                   : (nshower_rings == nrings);
                  },
                  {"nrings", "nshower_rings"})
          .Define("nmichel_electrons_cut",
                  [](const size_t nmichel_electrons) {
                    return nmichel_electrons == (is_mupi ? 1 : 0);
                  },
                  {"nmichel_electrons"})
          .Define("topo_cuts",
                  [](bool a, bool b, bool c) { return a && b && c; },
                  {"nrings_cut", "shower_ring_cut", "nmichel_electrons_cut"})
          .Define("rec_raw",
                  [](const NeutrinoEvent &event,
                     bool topo) -> std::optional<RecResult> {
                    if (!topo) {
                      return std::nullopt;
                    }
                    return event.Rec_lpi_event(is_mupi);
                  },
                  {"EventRecord", "topo_cuts"})
          .Define("truth",
                  [](const std::optional<RecResult> &rec_opt) -> EventRec {
                    if (!rec_opt.has_value()) {
                      return EventRec{};
                    }
                    const auto &rec = rec_opt.value();
                    auto lepton = rec.lepton.m_pair.first;
                    auto gamma_1 = rec.leading_gamma.m_pair.first;
                    auto gamma_2 = rec.subleading_gamma.transform(
                        [](const RingInfo &ring) { return ring.m_pair.first; });
                    return EventRec{lepton, gamma_1, gamma_2};
                  },
                  {"rec_raw"})
          .Define(
              "smeared",
              [](const std::optional<RecResult> &rec_opt) -> EventRec {
                if (!rec_opt.has_value()) {
                  return EventRec{};
                }
                const auto &rec = rec_opt.value();
                auto lepton = rec.lepton.m_pair.second;
                auto gamma_1 = rec.leading_gamma.m_pair.second;
                auto gamma_2 = rec.subleading_gamma.transform(
                    [](const RingInfo &ring) { return ring.m_pair.second; });
                return EventRec{lepton, gamma_1, gamma_2};
              },
              {"rec_raw"})
          .Define("kf_gammas_with_chi2",
                  [](const EventRec &smeared_opt)
                      -> std::optional<
                          std::tuple<std::array<momentum_t, 2>, double>> {
                    if (!smeared_opt.is_valid || !smeared_opt.has_gamma2) {
                      return std::nullopt;
                    }
                    decltype(kf_pi0(
                        {smeared_opt.gamma1, smeared_opt.gamma2})) best_fit =
                        std::nullopt;
                    for (int i = 0; i < 6; ++i) {
                      auto fit_result =
                          kf_pi0({smeared_opt.gamma1, smeared_opt.gamma2});
                      if (fit_result.has_value()) {
                        // choose the best fit based on smallest chi2
                        if (!best_fit.has_value() ||
                            std::get<1>(fit_result.value()) <
                                std::get<1>(best_fit.value())) {
                          best_fit = fit_result;
                        }
                      }
                    }
                    return best_fit;
                  },
                  {"smeared"})
          .Define("kf_gammas",
                  [](const std::optional<
                      std::tuple<std::array<momentum_t, 2>, double>>
                         &kf_gammas_with_chi2_opt)
                      -> std::optional<std::array<momentum_t, 2>> {
                    if (kf_gammas_with_chi2_opt.has_value()) {
                      return std::get<0>(kf_gammas_with_chi2_opt.value());
                    }
                    return std::nullopt;
                  },
                  {"kf_gammas_with_chi2"})
          .Define("kf_chi2",
                  [](const std::optional<
                      std::tuple<std::array<momentum_t, 2>, double>>
                         &kf_gammas_with_chi2_opt) -> double {
                    if (kf_gammas_with_chi2_opt.has_value()) {
                      return std::get<1>(kf_gammas_with_chi2_opt.value());
                    }
                    return -1.0;
                  },
                  {"kf_gammas_with_chi2"})
          .Define("kf_3d",
                  [](const EventRec &smeared_opt) -> gamma_dof {
                    if (!smeared_opt.is_valid || !smeared_opt.has_gamma2) {
                      return gamma_dof{};
                    }
                    decltype(kf_pi0_3D(
                        {smeared_opt.gamma1, smeared_opt.gamma2})) best_fit =
                        std::nullopt;
                    for (int i = 0; i < 6; ++i) {
                      auto fit_result =
                          kf_pi0_3D({smeared_opt.gamma1, smeared_opt.gamma2});
                      if (fit_result.has_value()) {
                        // choose the best fit based on smallest chi2
                        if (!best_fit.has_value() ||
                            std::get<1>(fit_result.value()) <
                                std::get<1>(best_fit.value())) {
                          best_fit = fit_result;
                        }
                      }
                    }
                    return best_fit
                        .transform([](const auto &t) { return std::get<0>(t); })
                        .value_or(gamma_dof{});
                  },
                  {"smeared"})
          .Define(
              "kf",
              [](const EventRec &smeared_opt,
                 const std::optional<std::array<momentum_t, 2>> &kf_gammas_opt)
                  -> EventRec {
                // return EventRec{};
                if (!smeared_opt.is_valid) {
                  return EventRec{};
                }
                if (kf_gammas_opt.has_value()) {
                  auto kf_gammas = kf_gammas_opt.value();
                  return EventRec{smeared_opt.lepton, kf_gammas[0],
                                  kf_gammas[1]};
                }
                return EventRec{};
              },
              {"smeared", "kf_gammas"})
          .Define("weight", []() { return 1.0; }, {});
  auto kf_report = df_all
                       .Filter([](size_t nrings) { return nrings == 3; },
                               {"nrings"}, "nrings_cut")
                       .Filter([](const EventRec &rec) { return rec.is_valid; },
                               {"kf"}, "valid KF")
                       .Report();
  df_all.Snapshot("outtree", output_path,
                  {// reconstructed variables
                   "truth", "smeared", "kf",
                   // raw variables
                   "raw_proton_momentum", "raw_mass_proton",
                   "final_state_system", "raw_pi0_before_fsi", "channel_name",
                   // cut and topology info
                   "nrings", "nshower_rings", "nmichel_electrons", "nrings_cut",
                   "shower_ring_cut", "nmichel_electrons_cut",
                   // dummy weight
                   "weight", "kf_chi2", "kf_3d"});
  kf_report->Print();
}
