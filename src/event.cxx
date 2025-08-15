#include "event.h"
#include "local_rand.h"
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include <TDatabasePDG.h>
#include <generator>
#include <smear.h>
#include <stdexcept>

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

const ROOT::Math::PxPyPzEVector &NeutrinoEvent::get_leading_det(int id) const {
  const ROOT::Math::PxPyPzEVector *leading = nullptr;
  for (const auto &p : det_range(id)) {
    if (leading == nullptr || p.second.P() > leading->P()) {
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

std::vector<std::pair<int, ROOT::Math::PxPyPzEVector>>
gen_decay(const std::pair<int, ROOT::Math::PxPyPzEVector> &to_decay) {
  std::vector<std::pair<int, ROOT::Math::PxPyPzEVector>> decay_products;
  auto &&[pdg, momentum] = to_decay;
  switch (pdg) {
  case 111: {
    // double mass_of_pion = momentum.M();
    double mass_of_pi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    double energy = mass_of_pi0 / 2.;
    auto boost_vec = momentum.BoostToCM();
    // the boost from CMS to LAB frame ( so the - sign )
    auto the_boost = ROOT::Math::Boost(-boost_vec);
    double costheta_in_CMS = get_thread_local_random().Uniform(0, 1);
    double theta = acos(costheta_in_CMS);
    double phi_in_CMS = get_thread_local_random().Uniform(0, 2 * M_PI);
    ROOT::Math::Polar3DVector first_photon_momentum{energy, theta, phi_in_CMS};
    auto another_photon_momentum = -first_photon_momentum;

    auto first_photon = the_boost(construct_4D(first_photon_momentum, energy));
    auto second_photon =
        the_boost(construct_4D(another_photon_momentum, energy));

    decay_products.push_back({22, first_photon});
    decay_products.push_back({22, second_photon});
  } break;

  default:
    decay_products.emplace_back(to_decay);
  }
  return decay_products;
}

void NeutrinoEvent::finalize_and_decay_in_detector() {
  for (auto &&[pdg, momentum_raw] : post) {
    auto decayed_particles = gen_decay({pdg, momentum_raw});
    for (auto &&[pdg, momentum] : decayed_particles) {
      before_smear.insert({pdg, momentum});
      auto stg = GetSmearStrategy(pdg);
      auto smared_momentum = stg ? stg->do_smearing(momentum) : momentum;
      in_detector.insert({pdg, smared_momentum});
    }
  }
}
