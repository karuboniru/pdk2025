#pragma once

#include <TLorentzVector.h>
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

class NeutrinoEvent {
public:
  [[nodiscard]] size_t count_in(int) const;
  [[nodiscard]] size_t count_out(int) const;
  [[nodiscard]] size_t count_post(int) const;

  void add_in(int, const TLorentzVector &);
  void add_out(int, const TLorentzVector &);
  void add_post(int, const TLorentzVector &);

  [[nodiscard]] auto in_range(int id) const {
    return equal_range_iterable(in, id);
  }
  [[nodiscard]] auto out_range(int id) const {
    return equal_range_iterable(out, id);
  }
  [[nodiscard]] auto post_range(int id) const {
    return equal_range_iterable(post, id);
  }

  [[nodiscard]] auto &get_in() const { return in; }
  [[nodiscard]] auto &get_out() const { return out; }
  [[nodiscard]] auto &get_post() const { return post; }

  const std::set<int> &get_ids_post() const;

  const TLorentzVector &get_leading(int id) const;

  [[nodiscard]] std::string get_channelname_no_nucleon() const;
  [[nodiscard]] std::string get_channelname() const;

private:
  std::unordered_multimap<int, TLorentzVector> in, out, post;
  std::set<int> ids_post;
};
