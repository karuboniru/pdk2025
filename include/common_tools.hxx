#pragma once
#include <array>
#include <cmath>
#include <span>

template <std::size_t N> double square_sum_array(const std::array<double, N> &arr) {
  double sum = 0.0;
  for (const auto &val : arr) {
    sum += val * val;
  }
  return sum;
}

template <std::size_t N> double abs_sum_array(const std::array<double, N> &arr) {
  double sum = 0.0;
  for (const auto &val : arr) {
    sum += std::abs(val);
  }
  return sum;
}

using para_view_t = std::span<const double>;

constexpr double chi2(double value, double center, double sigma) {
  return ((value - center) * (value - center)) / (sigma * sigma);
}
