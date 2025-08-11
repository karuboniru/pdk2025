#include "event.h"
#include <TDatabasePDG.h>
#include <stdexcept>

size_t NeutrinoEvent::count_in(int id) const { return in.count(id); }
size_t NeutrinoEvent::count_out(int id) const { return out.count(id); }
size_t NeutrinoEvent::count_post(int id) const { return post.count(id); }

void NeutrinoEvent::add_in(int id, const TLorentzVector &particle) {
  in.insert({id, particle});
}
void NeutrinoEvent::add_out(int id, const TLorentzVector &particle) {
  out.insert({id, particle});
}
void NeutrinoEvent::add_post(int id, const TLorentzVector &particle) {
  post.insert({id, particle});
  ids_post.insert(id);
}
const std::set<int> &NeutrinoEvent::get_ids_post() const { return ids_post; }

const TLorentzVector &NeutrinoEvent::get_leading(int id) const {
  const TLorentzVector *leading = nullptr;
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
