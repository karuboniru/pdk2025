#include "event.h"
#include "local_rand.h"
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include <TDatabasePDG.h>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <ranges>
#include <smear.h>
#include <stdexcept>

#include "decaygen.h"

size_t NeutrinoEvent::count_in(int id) const { return in.count(id); }
size_t NeutrinoEvent::count_out(int id) const { return out.count(id); }
size_t NeutrinoEvent::count_post(int id) const { return post.count(id); }
size_t NeutrinoEvent::count_det(int id) const { return in_detector.count(id); }

void NeutrinoEvent::add_in(int id, const ROOT::Math::PxPyPzEVector &particle) {
  in.insert({id, particle});
}
void NeutrinoEvent::add_out(int id, const ROOT::Math::PxPyPzEVector &particle) {
  out.insert({id, particle});
}
void NeutrinoEvent::add_post(int id,
                             const ROOT::Math::PxPyPzEVector &particle) {
  post.insert({id, particle});
  ids_post.insert(id);
}
const std::set<int> &NeutrinoEvent::get_ids_post() const { return ids_post; }
// const std::set<int> &NeutrinoEvent::get_ids_det() const { return
// ids_detector; }

const ROOT::Math::PxPyPzEVector &NeutrinoEvent::get_leading(int id) const {
  const ROOT::Math::PxPyPzEVector *leading = nullptr;
  for (const auto &p : post_range(id)) {
    if (leading == nullptr || p.second.P() > leading->P()) {
      leading = &p.second;
    }
  }
  if (leading == nullptr) {
    throw std::runtime_error("No leading particle found");
  }
  return *leading;
}

const std::pair<ROOT::Math::PxPyPzEVector, ROOT::Math::PxPyPzEVector> &
NeutrinoEvent::get_leading_det(int id) const {

  const std::pair<ROOT::Math::PxPyPzEVector, ROOT::Math::PxPyPzEVector>
      *leading = nullptr;
  for (const auto &p : det_range(id)) {
    if (leading == nullptr || p.second.second.P() > leading->second.P()) {
      leading = &p.second;
    }
  }
  if (leading == nullptr) {
    throw std::runtime_error("No leading particle found");
  }
  return *leading;
}
auto &db = *TDatabasePDG::Instance();
std::string NeutrinoEvent::get_channelname() const {
  std::string channel_name;
  std::stringstream ss{};
  for (const auto &id : ids_post) {
    auto count = count_post(id);
    if (id == -11) {
      continue;
    }
    if (count > 1) {
      ss << count;
    }
    ss << db.GetParticle(id)->GetName() << " ";
  }
  channel_name = ss.str();
  channel_name = channel_name.substr(0, channel_name.size() - 1);
  return channel_name;
}

std::string NeutrinoEvent::get_channelname_no_nucleon() const {
  auto is_this_transparent = is_transparent();
  if (is_this_transparent) {
    return "transparent";
  }
  std::string channel_name_no_nucleon;
  std::stringstream ss{};
  ss << "e+ ";
  for (const auto &id : ids_post) {
    auto count = count_post(id);
    if (id == 2212 || id == 2112 || id == -11) {
      continue;
    }
    if (count > 1) {
      ss << count;
    }
    ss << db.GetParticle(id)->GetName() << " ";
  }
  channel_name_no_nucleon = ss.str();
  channel_name_no_nucleon =
      channel_name_no_nucleon.substr(0, channel_name_no_nucleon.size() - 1);
  return channel_name_no_nucleon;
}

ROOT::Math::PxPyPzEVector construct_4D(const auto &vec, double E) {
  double px = vec.x();
  double py = vec.y();
  double pz = vec.z();
  return ROOT::Math::PxPyPzEVector(px, py, pz, E);
}

void NeutrinoEvent::finalize_and_decay_in_detector() {
  n_michel_electrons = 0;
  for (const auto &[pdg, momentum_raw] : post) {
    auto decayed_particles =
        EvtGenInterface::get_instance().decay({pdg, momentum_raw});
    for (auto &&[pdg, momentum] : decayed_particles) {
      auto stg = GetSmearStrategy(pdg);
      auto smared_momentum = stg ? stg->do_smearing(momentum) : momentum;
      pair_momentum_t momentum_pair{momentum, smared_momentum};
      in_detector.insert({pdg, momentum_pair});
      if (std::abs(pdg) == 11 ||
          pdg == 22) { // e+ , e- , gamma -> shower-like ring
        rings_in_detector.emplace_back(momentum_pair, pdg, true);
      }
      if (std::abs(pdg) == 13 &&
          momentum.P() > 0.118) { // mu+ , mu- -> track-like ring, with Michel
        rings_in_detector.emplace_back(momentum_pair, pdg, false);
        n_michel_electrons++;
      }
      if (pdg == 211 &&
          momentum.P() > 0.156) { // pi+ -> track-like ring, with Michel
        rings_in_detector.emplace_back(momentum_pair, pdg, false);
        n_michel_electrons++;
      }
      if (pdg == -211 &&
          momentum.P() > 0.156) { // pi- -> track-like ring, without Michel
        rings_in_detector.emplace_back(momentum_pair, pdg, false);
      }
      if (std::abs(pdg) == 321) {
        // will be rings from kaon decay.... but ignore for now
        n_michel_electrons++;
      }
    }
  }

  post_process_rings_in_detector();
}
void NeutrinoEvent::post_process_rings_in_detector() {
  // loop over all rings, merge rings that are within 15 degrees
  // and mark the merged rings to be swiped
  std::vector<size_t> disjoint_set(rings_in_detector.size());
  auto find_root = [&](this auto &&self, size_t id) -> size_t {
    if (disjoint_set[id] != id) {
      disjoint_set[id] = self(disjoint_set[id]);
    }
    return disjoint_set[id];
  };
  std::iota(disjoint_set.begin(), disjoint_set.end(), 0);
  for (auto &&[id1, ring1] : rings_in_detector | std::views::enumerate) {
    for (auto &&[id2, ring2] : rings_in_detector | std::views::enumerate |
                                   std::views::drop(id1 + 1)) {
      const auto &&ring1_dir_true = ring1.m_pair.first.Vect().Unit();
      const auto &&ring2_dir_true = ring2.m_pair.first.Vect().Unit();
      double cos_angle = ring1_dir_true.Dot(ring2_dir_true);
      const double cos_threshold = std::cos(15.0 * M_PI / 180.0);
      if (cos_angle > cos_threshold) {
        // merge these two rings
        auto root1 = find_root(id1);
        auto root2 = find_root(id2);
        if (root1 != root2) {
          disjoint_set[root2] = root1;
        }
      }
    }
  }

  for (auto &&[id, ring] : rings_in_detector | std::views::enumerate) {
    auto root = find_root(id);
    auto &root_ring = rings_in_detector[root];
    if (root != id) {
      // merge to root
      root_ring.m_pair =
          root_ring.m_pair + ring.m_pair; // sum the momentum pairs
      root_ring.is_shower |= ring.is_shower;

      ring.to_remove = true;
    }
  }

  for (auto &ring : rings_in_detector) {
    if (ring.m_pair.first.P() < 0.030) {
      ring.to_remove = true;
    }
  }

  // finally remove the marked rings
  std::ranges::remove_if(rings_in_detector,
                         [](const RingInfo &ring) { return ring.to_remove; });
}

bool NeutrinoEvent::is_transparent() const {
  // this checks if the after FSI is "almost" the same as before FSI
  // there should be no particle produced or absorbed in the nucleus
  // the momentum should be almost the same (1e-3 GeV/c)
  if (out.size() != post.size()) {
    return false;
  }
  // each event that seen in post must be in out
  // and the momentum must be almost the same
  return std::ranges::all_of(post, [this](const auto &p) {
    auto range = out_range(p.first);
    return std::ranges::any_of(range, [&p](const auto &q) {
      auto diff_p = p.second - q.second;
      return diff_p.P() < 1e-3;
    });
  });
}

RecResult NeutrinoEvent::Rec_lpi_event(bool is_mu_pi) const {
  RecResult result{};
  switch (rings_in_detector.size()) {
  case 3: {
    size_t best_lepton_candidate = 0;
    double best_pi0_candidate_m{};
    for (size_t i = 0; i < 3; ++i) {
      // if on epi mode, the lepton must be a muon
      // aka. non-shower ring
      if (is_mu_pi && !rings_in_detector[i].is_shower) {
        continue;
      }
      auto sum_p4_rec_pi0 = ROOT::Math::PxPyPzEVector{};
      for (size_t j = 0; j < 3; ++j) {
        if (j != i) {
          sum_p4_rec_pi0 += rings_in_detector[j].m_pair.second;
        }
      }
      if (std::abs(sum_p4_rec_pi0.M() - 0.135) <
          std::abs(best_pi0_candidate_m - 0.135)) {
        best_pi0_candidate_m = sum_p4_rec_pi0.M();
        best_lepton_candidate = i;
      }
    }
    auto &gamma1 = rings_in_detector[(best_lepton_candidate + 1) % 3];
    auto &gamma2 = rings_in_detector[(best_lepton_candidate + 2) % 3];
    result.lepton = rings_in_detector[best_lepton_candidate];
    result.leading_gamma =
        gamma1.m_pair.second.P() > gamma2.m_pair.second.P() ? gamma1 : gamma2;
    result.subleading_gamma =
        gamma1.m_pair.second.P() > gamma2.m_pair.second.P() ? gamma2 : gamma1;
    result.rec_pi0 = gamma1.m_pair + gamma2.m_pair;
  } break;
  case 2: {
    // now we have 2 rings, usually one is e+/mu+ another is gamma
    // in case of epi, it does not matter which is which
    size_t lepton_index = 0;
    if (is_mu_pi && rings_in_detector[0].is_shower) {
      lepton_index = 1;
    }
    result.lepton = rings_in_detector[lepton_index];
    result.rec_pi0 = rings_in_detector[1 - lepton_index].m_pair;
    result.leading_gamma = rings_in_detector[1 - lepton_index];
  } break;

  default:
    throw std::runtime_error("Rec_lpi_event: number of rings is not 2 or 3");
  }

  return result;
}

momentum_pair operator+(const momentum_pair &a, const momentum_pair &b) {
  return {a.first + b.first, a.second + b.second};
}
