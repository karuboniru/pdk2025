#pragma once

#include <Math/LorentzVector.h>
#include <Math/Vector4Dfwd.h>
#include <set>
#include <unordered_map>
#include <utility>

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

using pair_momentum_t = std::pair<ROOT::Math::PxPyPzEVector, ROOT::Math::PxPyPzEVector>;

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

private:
  std::unordered_multimap<int, ROOT::Math::PxPyPzEVector> in, out, post;
  std::unordered_multimap<
      int, pair_momentum_t>
      in_detector;
  std::set<int> ids_post;
};
