#include "EvtTracker2event.h"
#include "event.h"
#include <Math/LorentzVector.h>
ROOT::RDF::RNode TrackerPrepare(ROOT::RDF::RNode df) {
  return df.Define("EventRecord",
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
                   {"nparticles", "pdg", "status", "P"});
}

ROOT::RDF::RNode TrackerPrepareNeutrino(ROOT::RDF::RNode df) {
  return df.Define("EventRecord",
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
                   {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4"});
}
