#pragma once

#include <Math/LorentzVector.h>
#include <Math/Vector4Dfwd.h>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

using pair_momentum_t =
    std::pair<ROOT::Math::PxPyPzEVector, ROOT::Math::PxPyPzEVector>;
using momentum_pair = pair_momentum_t;

// struct momentum_pair {
//   ROOT::Math::PxPyPzEVector true_momentum;
//   ROOT::Math::PxPyPzEVector momentum;
// };

momentum_pair operator+(const momentum_pair &a, const momentum_pair &b);

struct RingInfo {
  momentum_pair m_pair;
  int from_pdg{};
  bool is_shower{};
};

struct RecResult {
  RingInfo lepton;
  RingInfo leading_gamma;
  std::optional<RingInfo> subleading_gamma;
  std::optional<momentum_pair> rec_pi0;
};

template <typename U> class equal_range_iterable {
public:
  template <typename T>
  equal_range_iterable(T &&m, int id) : range(m.equal_range(id)) {}

  [[nodiscard]] auto &&begin() { return range.first; }

  [[nodiscard]] auto &&end() { return range.second; }

  [[nodiscard]] auto size() const {
    return std::distance(range.first, range.second);
  }

private:
  decltype(std::declval<U>().equal_range(std::declval<int>())) range;
};

template <typename U>
equal_range_iterable(U &&, int) -> equal_range_iterable<U>;



class NeutrinoEvent {
public:
  [[nodiscard]] size_t count_in(int) const;
  [[nodiscard]] size_t count_out(int) const;
  [[nodiscard]] size_t count_post(int) const;
  [[nodiscard]] size_t count_det(int) const;

  void add_in(int, const ROOT::Math::PxPyPzEVector &);
  void add_out(int, const ROOT::Math::PxPyPzEVector &);
  void add_post(int, const ROOT::Math::PxPyPzEVector &);

  [[nodiscard]] auto in_range(int id) const {
    return equal_range_iterable(in, id);
  }
  [[nodiscard]] auto out_range(int id) const {
    return equal_range_iterable(out, id);
  }
  [[nodiscard]] auto post_range(int id) const {
    return equal_range_iterable(post, id);
  }
  [[nodiscard]] auto det_range(int id) const {
    return equal_range_iterable(in_detector, id);
  }

  [[nodiscard]] auto &get_in() const { return in; }
  [[nodiscard]] auto &get_out() const { return out; }
  [[nodiscard]] auto &get_post() const { return post; }
  [[nodiscard]] auto &get_det() const { return in_detector; }

  const std::set<int> &get_ids_post() const;
  // const std::set<int> &get_ids_det() const;

  const ROOT::Math::PxPyPzEVector &get_leading(int id) const;
  const pair_momentum_t &get_leading_det(int id) const;

  [[nodiscard]] std::string get_channelname_no_nucleon() const;
  [[nodiscard]] std::string get_channelname() const;

  void finalize_and_decay_in_detector();

  bool is_transparent() const;

  RecResult Rec_lpi_event(bool is_mu_pi = false) const;

  size_t count_rings_in_detector() const { return rings_in_detector.size(); }
  size_t count_shower_rings_in_detector() const {
    return std::count_if(rings_in_detector.begin(), rings_in_detector.end(),
                         [](const RingInfo &ring) { return ring.is_shower; });
  }
  size_t get_n_michel_electrons() const { return n_michel_electrons; }

private:
  std::unordered_multimap<int, ROOT::Math::PxPyPzEVector> in, out, post;
  std::unordered_multimap<int, pair_momentum_t> in_detector;
  std::set<int> ids_post;
  std::vector<RingInfo> rings_in_detector;
  size_t n_michel_electrons{0};
};
