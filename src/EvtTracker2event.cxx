#include "EvtTracker2event.h"
#include "event.h"
#include <Math/LorentzVector.h>

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
      .Define("nmichel_electrons",
              [](const NeutrinoEvent &event) {
                return event.get_n_michel_electrons();
              },
              {"EventRecord"});
  ;
}

ROOT::RDF::RNode TrackerPrepare(ROOT::RDF::RNode df) {
  return general_define(
      df.Define("EventRecord",
                [](int StdHepN, const ROOT::RVec<int> &StdHepPdg,
                   const ROOT::RVec<int> &StdHepStatus,
                   const ROOT::RVec<double> &StdHepP4_) {
                  NeutrinoEvent e{};
                  for (int i = 0; i < StdHepN; i++) {
                    ROOT::Math::PxPyPzEVector p4{
                        StdHepP4_[(i * 4) + 1], StdHepP4_[(i * 4) + 2],
                        StdHepP4_[(i * 4) + 3], StdHepP4_[(i * 4UL)]};
                    switch (StdHepStatus[i]) {
                    case 0:
                      e.add_in(StdHepPdg[i], p4);
                      break;
                    case 1:
                      e.add_post(StdHepPdg[i], p4);
                      break;
                    case 2:
                      e.add_out(StdHepPdg[i], p4);
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

ROOT::RDF::RNode TrackerPrepareNeutrino(ROOT::RDF::RNode df) {
  return general_define(
      df.Define("EventRecord",
                [](int StdHepN, const ROOT::RVec<int> &StdHepPdg,
                   const ROOT::RVec<int> &StdHepStatus,
                   const ROOT::RVec<double> &StdHepP4_) {
                  NeutrinoEvent e{};
                  for (int i = 0; i < StdHepN; i++) {
                    ROOT::Math::PxPyPzEVector p4{
                        StdHepP4_[(i * 4UL)], StdHepP4_[(i * 4) + 1],
                        StdHepP4_[(i * 4) + 2], StdHepP4_[(i * 4) + 3]};
                    switch (StdHepStatus[i]) {
                    case 0:
                      e.add_in(StdHepPdg[i], p4);
                      break;
                    case 1:
                      e.add_post(StdHepPdg[i], p4);
                      break;
                    case 2:
                      e.add_out(StdHepPdg[i], p4);
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
