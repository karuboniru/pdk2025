#pragma once

#include "local_rand.h"
#include "common_tools.hxx"
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <ranges>

template <class T>
class SingleKFLagMulDoubleMin final : public ROOT::Minuit2::FCNBase, public T {
public:
  static constexpr size_t lagrange_parmaeter_count_per_constrain = 1;
  static constexpr size_t constrain_count = T::constrain_count;
  static constexpr size_t total_lagrange_parameter_count =
      lagrange_parmaeter_count_per_constrain * constrain_count;

  class LagrangeMinimizerInner : public ROOT::Minuit2::FCNBase {
  public:
    template <typename U>
    LagrangeMinimizerInner(const SingleKFLagMulDoubleMin<T> *fcn_, U &&lambda_)
        : fcn(fcn_), lambda(lambda_) {}
    LagrangeMinimizerInner(const SingleKFLagMulDoubleMin<T> *fcn_,
                           const std::vector<double> &lambda_)
        : fcn(fcn_), lambda() {
      if (lambda_.size() != constrain_count) {
        throw std::runtime_error("LagrangeMinimizerInner: invalid lambda size");
      }
      for (size_t i = 0; i < constrain_count; ++i) {
        lambda[i] = lambda_[i];
      }
    }

    [[nodiscard]] double
    operator()(const std::vector<double> &kine_paras) const override {
      if (kine_paras.size() != T::kinematic_parameter_count) {
        throw std::runtime_error(
            "LagrangeMinimizerInner: invalid parameter size");
      }
      auto constrains = fcn->get_constrain(kine_paras);
      double penalty = fcn->get_parameter_penalty(kine_paras);
      for (const auto &[lambda, constrain] :
           std::views::zip(lambda, constrains)) {
        penalty += lambda * constrain + 0.5 * constrain * constrain;
      }
      return penalty;
    }

    [[nodiscard]] double Up() const override { return 1.0; }

  private:
    const SingleKFLagMulDoubleMin<T> *fcn;
    std::array<double, constrain_count> lambda;
  };

  using T::get_constrain;
  using T::get_measured_from_parameters;
  using T::get_parameter_penalty;
  using T::T;

  static auto get_kinematic_parameters(std::span<const double> all_parameters) {
    return all_parameters;
  }

  std::optional<std::vector<double>>
  argmin_x_lagrangian(const std::vector<double> &lambdas) const {
    LagrangeMinimizerInner inner_fcn(this, lambdas);
    ROOT::Minuit2::MnUserParameters upar = T::generate_initial_parameters({});
    ROOT::Minuit2::MnMigrad migrad(inner_fcn, upar);
    ROOT::Minuit2::FunctionMinimum min = migrad();
    if (min.IsValid()) {
      return min.UserParameters().Params();
    }
    return std::nullopt;
  }

  double operator()(const std::vector<double> &params) const override {
    auto argminx = argmin_x_lagrangian(params);
    if (!argminx) {
      return 1e20;
    }
    auto constrains = get_constrain((std::span<const double>(*argminx)));
    return abs_sum_array(constrains);
  }

  double Up() const override { return 1.0; }

  ROOT::Minuit2::MnUserParameters generate_initial_parameters() const {
    ROOT::Minuit2::MnUserParameters upar;
    auto &rand = get_thread_local_random();
    for (size_t i = 0; i < constrain_count; ++i) {
      upar.Add("lambda_" + std::to_string(i), rand.Gaus(0.0, 1.0), 0.1);
    }
    return upar;
  }

  static std::optional<std::tuple<typename T::answer_t, double>>
  do_kinematics_fit(const T::answer_t &input) {
    constexpr double tol = 1e-4;
    SingleKFLagMulDoubleMin<T> fcn(input);
    auto upar = fcn.generate_initial_parameters();
    for (size_t iter{};; iter++) {
      ROOT::Minuit2::MnMigrad migrad(fcn, upar);
      ROOT::Minuit2::FunctionMinimum min = migrad();
      if (min.IsValid()) {
        auto params = min.UserParameters().Params();
        auto argminx = fcn.argmin_x_lagrangian(params);
        if (!argminx) {
          // bad fit, should not happen
          throw std::runtime_error(
              "SingleKFLagMulDoubleMin: invalid argminx after valid min");
        }
        auto kinematics_params = (std::span<const double>(*argminx));
        auto constrains = fcn.get_constrain(kinematics_params);
        // if converged, return result
        if (double penalty = fcn.get_parameter_penalty(kinematics_params);
            std::sqrt(square_sum_array(constrains)) < tol &&
            penalty < T::chi2_cut) {
          return std::make_tuple(
              fcn.get_measured_from_parameters(kinematics_params), penalty);
        }
      } else {
        // bad fit, regenerate initial parameters
        // std::println("Bad fit in iteration {}", iter);
        upar = fcn.generate_initial_parameters();
      }
      if (iter > 16) {
        // std::println("Failed to converge after {} iterations", iter);
        return std::nullopt;
      }
    }
  }
};
