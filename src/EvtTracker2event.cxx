#include "EvtTracker2event.h"
#include "RtypesCore.h"
#include "event.h"
#include "local_rand.h"
#include <Math/LorentzVector.h>
#include <cmdline.h>
#include <cstdlib>
#include <print>
#include <ranges>
#include <span>

double p_no_neutron_tag(size_t n_neutron) {
  return std::pow(1 - n_tagging_eff, n_neutron);
}

size_t n_michel_electrons_tagged(size_t nmichel_electrons_raw) {
  // assuming tagging efficiency of 88%, binomial distribution
  auto &rng = get_thread_local_random();
  return rng.Binomial(nmichel_electrons_raw, 0.88);
}

auto general_define(auto &&df) {
  return df
      .Define("nrings",
              [](const NeutrinoEvent &event) {
                return event.count_rings_in_detector();
              },
              {"EventRecord"})
      .Define("nshower_rings",
              [](const NeutrinoEvent &event) {
                return event.count_shower_rings_in_detector();
              },
              {"EventRecord"})
      .Define("nmichel_electrons_raw",
              [](const NeutrinoEvent &event) {
                return event.get_n_michel_electrons();
              },
              {"EventRecord"})
      .Define("nmichel_electrons", n_michel_electrons_tagged,
              {"nmichel_electrons_raw"})
      .Define("n_neutron", [](NeutrinoEvent &e) { return e.count_post(2112); },
              {"EventRecord"})
      .Define("p_no_neutron_tag", p_no_neutron_tag, {"n_neutron"});
}

ROOT::RDF::RNode TrackerPrepare(ROOT::RDF::RNode df) {
  return general_define(df.Define(
      "EventRecord",
      [](int StdHepN, const ROOT::RVec<int> &StdHepPdg,
         const ROOT::RVec<int> &StdHepStatus,
         const ROOT::RVec<double> &StdHepP4_) {
        NeutrinoEvent e{};
        for (const auto &[pdg, status, p4] :
             std::views::zip(StdHepPdg, StdHepStatus,
                             std::ranges::iota_view(StdHepP4_.data()) |
                                 std::views::chunk(4) |
                                 std::views::transform([](auto &&chunk) {
                                   return ROOT::Math::PxPyPzEVector{
                                       *chunk[1] * unit, *chunk[2] * unit,
                                       *chunk[3] * unit, *chunk[0] * unit};
                                 }))) {
          switch (status) {
          case 0:
            e.add_in(pdg, p4);
            break;
          case 1:
            if (nofsi)
              break;
            e.add_post(pdg, p4);
            break;
          case 2:
            e.add_out(pdg, p4);
            if (nofsi) {
              e.add_post(pdg, p4);
            }
            break;
          default:
            throw;
          }
        }
        e.finalize_and_decay_in_detector();
        return e;
      },
      {"nparticles", "pdg", "status", "P"}));
}

ROOT::RDF::RNode TrackerPrepareGENIE(ROOT::RDF::RNode df) {
  return general_define(df.Define(
      "EventRecord",
      [](const ROOT::RVec<int> &StdHepPdg, const ROOT::RVec<int> &StdHepStatus,
         const ROOT::RVec<double> &StdHepP4_) {
        NeutrinoEvent e{};
        for (auto &&[pdg, status, p4] :
             std::views::zip(StdHepPdg, StdHepStatus,
                             std::ranges::iota_view(StdHepP4_.data()) |
                                 std::views::chunk(4) |
                                 std::views::transform([](auto &&chunk) {
                                   return ROOT::Math::PxPyPzEVector{
                                       *chunk[0] * unit, *chunk[1] * unit,
                                       *chunk[2] * unit, *chunk[3] * unit};
                                 }))) {
          switch (status) {
          case 3:
            e.add_in(pdg, p4);
            break;
          case 1:
            if (!nofsi) {
              e.add_post(pdg, p4);
            }
            if (auto abs_pdg = std::abs(pdg);
                // if  StdHepPdg[0] == 2212 then before FSI is after FSI
                // and leptons don't undergo FSI
                StdHepPdg[0] == 2212 || abs_pdg == 11 || abs_pdg == 13) {
              e.add_out(pdg, p4);
              if (nofsi) {
                e.add_post(pdg, p4);
              }
            }
            break;
          case 14:
            e.add_out(pdg, p4);
            if (nofsi) {
              e.add_post(pdg, p4);
            }
            break;
          default:
            break;
          }
        }

        e.finalize_and_decay_in_detector();
        return e;
      },
      {"StdHepPdg", "StdHepStatus", "StdHepP4"}));
}

ROOT::RDF::RNode TrackerPrepareNeutrino(ROOT::RDF::RNode df) {
  return general_define(df.Define(
      "EventRecord",
      [](int StdHepN, const ROOT::RVec<int> &StdHepPdg,
         const ROOT::RVec<int> &StdHepStatus,
         const ROOT::RVec<double> &StdHepP4_) {
        NeutrinoEvent e{};
        for (auto &&[pdg, status, p4] : std::views::zip(
                 StdHepPdg, StdHepStatus,
                 std::ranges::iota_view(StdHepP4_.data()) |
                     std::views::chunk(4) |
                     std::views::transform(
                         [](auto &&chunk) -> ROOT::Math::PxPyPzEVector {
                           return {*chunk[0] * unit, *chunk[1] * unit,
                                   *chunk[2] * unit, *chunk[3] * unit};
                         }))) {
          switch (status) {
          case 0:
            e.add_in(pdg, p4);
            break;
          case 1:
            if (nofsi)
              break;
            e.add_post(pdg, p4);
            break;
          case 2:
            e.add_out(pdg, p4);
            if (nofsi) {
              e.add_post(pdg, p4);
            }
            break;
          default:
            throw;
          }
        }
        e.finalize_and_decay_in_detector();
        return e;
      },
      {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4"}));
}
