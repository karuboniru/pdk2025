#pragma once
#include <array>
#include <cmath>
#include <span>
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h>

template <std::size_t N>
double square_sum_array(const std::array<double, N> &arr) {
  double sum = 0.0;
  for (const auto &val : arr) {
    sum += val * val;
  }
  return sum;
}

template <std::size_t N>
double abs_sum_array(const std::array<double, N> &arr) {
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

size_t guess_nproc_from_env();

double rayleigh_log_likelihood(double value, double sigma);

double rayleigh_log_likelihood_normalized(double value, double sigma);

using momentum_t = ROOT::Math::PxPyPzEVector;
double angle_between(const momentum_t &p4_1, const momentum_t &p4_2);
