#pragma once

#include "common_tools.hxx"

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <ranges>

template <class T, double hyper_parameter_initial_mu = 1.0,
          double hyper_parameter_initial_lambda = 0.0,
          double hyper_parameter_increase_factor = 1.5>
class SingleKFLagMul final : public ROOT::Minuit2::FCNBase, public T {
  static constexpr size_t lagrange_parmaeter_count_per_constrain = 2;
  static constexpr size_t constrain_count = T::constrain_count;

  struct hyper_parameters_t {
    double lambda = hyper_parameter_initial_lambda;
    double mu = hyper_parameter_initial_mu;
  };

public:
  std::array<hyper_parameters_t, constrain_count> hyper_parameters;

  using T::get_constrain;
  using T::get_measured_from_parameters;
  using T::get_parameter_penalty;
  using T::T;

  double operator()(const std::vector<double> &params) const override {
    if (params.size() != T::kinematic_parameter_count) {
      throw std::runtime_error("SingleKFLagMul: invalid parameter size");
    }
    auto constraints = T::get_constrain(params);
    double penalty = T::get_parameter_penalty(params);
    for (auto &&[lag_paras, constraint] :
         std::views::zip(hyper_parameters, constraints)) {
      penalty += (lag_paras.lambda * constraint) +
                 (0.5 * lag_paras.mu * constraint * constraint);
    }
    return penalty;
  }

  [[nodiscard]] double Up() const override { return 1.0; }

  [[nodiscard]] ROOT::Minuit2::MnUserParameters
  generate_initial_parameters() const {
    ROOT::Minuit2::MnUserParameters upar;
    return T::generate_initial_parameters(upar);
  }

  static std::optional<std::tuple<typename T::answer_t, double>>
  do_kinematics_fit(const T::answer_t &input) {
    constexpr double tol = 1e-5;
    SingleKFLagMul<T> fcn(input);
    auto upar = fcn.generate_initial_parameters();
    unsigned int default_strategy = 1;
    for (size_t iter{};; iter++) {
      ROOT::Minuit2::MnMigrad migrad(
          fcn, upar, ROOT::Minuit2::MnStrategy{default_strategy});
      ROOT::Minuit2::FunctionMinimum min = migrad();
      if (min.IsValid()) {
        auto params = min.UserParameters().Params();
        auto constrains = fcn.get_constrain(params);
        // if converged, return result
        if (double penalty = fcn.get_parameter_penalty(params);
            std::sqrt(square_sum_array(constrains)) < tol &&
            penalty < T::chi2_cut) {
          return std::make_tuple(fcn.get_measured_from_parameters(params),
                                 penalty);
        }
        for (auto &&[hyper_paras, constrain] :
             std::views::zip(fcn.hyper_parameters, constrains)) {
          // update hyper parameters
          hyper_paras.lambda += hyper_paras.mu * constrain;
          hyper_paras.mu *= hyper_parameter_increase_factor;
        }
      } else {
        // bad fit, regenerate initial parameters
        default_strategy = std::min(2U, default_strategy + 1);
        fcn.hyper_parameters = {}; // reset hyper parameters
        upar = fcn.generate_initial_parameters();
      }
      if (iter > 256) {
        return std::nullopt;
      }
    }
  }
};
