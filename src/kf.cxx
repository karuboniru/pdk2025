#include "kf.h"
#include "local_rand.h"
#include "smear.h"
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <cstddef>
#include <optional>
#include <print>
#include <ranges>

double rayleigh_log_likelihood(double value, double sigma) {
  if (value < 0 || sigma <= 0) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::log(value) - ((value * value) / (2 * sigma * sigma));
}

double rayleigh_log_likelihood_normalized(double value, double sigma) {
  return -2. * (rayleigh_log_likelihood(value, sigma) -
                rayleigh_log_likelihood(sigma, sigma));
}

ROOT::Math::PxPyPzEVector
change_vector(const ROOT::Math::PxPyPzEVector &original, double dx, double dy,
              double dz) {
  ROOT::Math::XYZVector vec3d{original.x() + dx, original.y() + dy,
                              original.z() + dz};
  double energy = vec3d.R();
  return ROOT::Math::PxPyPzEVector{vec3d.x(), vec3d.y(), vec3d.z(), energy};
}

constexpr double chi2(double value, double center, double sigma) {
  return ((value - center) * (value - center)) / (sigma * sigma);
}

template <size_t N> double square_sum_array(const std::array<double, N> &arr) {
  double sum = 0.0;
  for (const auto &val : arr) {
    sum += val * val;
  }
  return sum;
}

using para_view_t = std::span<const double>;

class Problem6D {
public:
  static constexpr size_t kinematic_parameter_count = 6;
  static constexpr size_t constrain_count = 3;
  using answer_t = std::array<momentum_t, 2>;
  Problem6D(const answer_t &measured_) : measured(measured_) {
    auto gamma_smear_strategy = GetSmearStrategy(22);
    if (!gamma_smear_strategy) {
      throw std::runtime_error("No smear strategy for gamma (PDG 22)");
    }
    for (size_t i = 0; i < 2; ++i) {
      double e = measured[i].E();
      sigma_angle[i] = gamma_smear_strategy->get_sigma_angle(e) * M_PI / 180.0;
      sigma_momentum[i] = gamma_smear_strategy->get_sigma_energy(e);
    }
  }

  std::array<double, constrain_count>
  get_constrain(const para_view_t &params_kin) const {
    std::array<double, constrain_count> constrains;
    auto fitted = get_measured_from_parameters(params_kin);
    auto pi0_momentum = fitted[0] + fitted[1];
    constrains[0] = pi0_momentum.M() - 0.134977;
    auto fitted_pi0_momentum = (fitted[0].Vect() + fitted[1].Vect()).Unit();
    auto raw_pi0_momentum = (measured[0] + measured[1]).Vect().Unit();
    constrains[1] = fitted_pi0_momentum.x() - raw_pi0_momentum.x();
    constrains[2] = fitted_pi0_momentum.y() - raw_pi0_momentum.y();
    return constrains;
  }

  double get_parameter_penalty(const para_view_t &params_kin) const {
    auto fitted = get_measured_from_parameters(params_kin);
    double penalty = 0.0;
    for (auto &&[orig, fit, s_ang, s_mom] :
         std::views::zip(measured, fitted, sigma_angle, sigma_momentum)) {
      // angle penalty
      ROOT::Math::XYZVector orig_dir = orig.Vect().Unit();
      ROOT::Math::XYZVector fitted_dir = fit.Vect().Unit();
      double angle_diff = std::acos(orig_dir.Dot(fitted_dir));
      penalty += chi2(angle_diff, 0.0, s_ang);

      // momentum penalty
      double orig_mom = orig.P();
      double fitted_mom = fit.P();
      penalty += chi2(fitted_mom / orig_mom, 1., s_mom);
    }
    return penalty;
  }

  auto generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
    auto &local_rand = get_thread_local_random();
    for (size_t i = 0; i < kinematic_parameter_count; ++i) {
      from.Add(std::format("dx{}", i), local_rand.Gaus(0, 0.05), 0.001);
    }
    return from;
  }

  answer_t get_measured_from_parameters(const para_view_t &params_kin) const {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error(
          "Problem6D::get_measured_from_parameters: invalid parameter size");
    }
    answer_t ret;
    for (auto &&[target, source, p3] :
         std::views::zip(ret, measured, params_kin | std::views::chunk(3))) {
      double dx = p3[0];
      double dy = p3[1];
      double dz = p3[2];
      target = change_vector(source, dx, dy, dz);
    }
    return ret;
  }

private:
  answer_t measured;
  std::array<double, 2> sigma_angle{};
  std::array<double, 2> sigma_momentum{};
};

template <class T>
class SingleKFLagMul final : public ROOT::Minuit2::FCNBase, public T {
  static constexpr size_t lagrange_parmaeter_count_per_constrain = 2;
  static constexpr size_t constrain_count = T::constrain_count;
  static constexpr size_t total_lagrange_parameter_count =
      lagrange_parmaeter_count_per_constrain * constrain_count;

public:
  using T::get_constrain;
  using T::get_measured_from_parameters;
  using T::get_parameter_penalty;
  using T::T;

  static auto get_kinematic_parameters(std::span<const double> all_parameters) {
    return all_parameters.subspan(total_lagrange_parameter_count);
  }

  double operator()(const std::vector<double> &params) const override {
    if (params.size() !=
        total_lagrange_parameter_count + T::kinematic_parameter_count) {
      throw std::runtime_error("SingleKFLagMul: invalid parameter size");
    }
    auto span_params =
        get_kinematic_parameters(std::span<const double>(params));
    auto constraints = T::get_constrain(span_params);
    double penalty = T::get_parameter_penalty(span_params);
    for (auto &&[lag_paras, constraint] : std::views::zip(
             std::views::chunk(params, lagrange_parmaeter_count_per_constrain),
             constraints)) {
      double lagrange_multiplier = lag_paras[0];
      double multiplyer = lag_paras[1];
      penalty += lagrange_multiplier * constraint +
                 0.5 * multiplyer * constraint * constraint;
    }
    return penalty;
  }

  double Up() const override { return 1.0; }

  ROOT::Minuit2::MnUserParameters generate_initial_parameters() const {
    ROOT::Minuit2::MnUserParameters upar;
    for (size_t i = 0; i < constrain_count; ++i) {
      upar.Add("lambda_" + std::to_string(i), 0.0);
      upar.Add("mu_" + std::to_string(i), 1.0);
    }
    return T::generate_initial_parameters(upar);
  }

  static std::optional<std::tuple<typename T::answer_t, double>>
  do_kinematics_fit(const T::answer_t &input) {
    constexpr double tol = 1e-5;
    SingleKFLagMul<T> fcn(input);
    auto upar = fcn.generate_initial_parameters();
    for (size_t iter{};; iter++) {
      ROOT::Minuit2::MnMigrad migrad(fcn, upar, ROOT::Minuit2::MnStrategy{2});
      ROOT::Minuit2::FunctionMinimum min = migrad();
      if (min.IsValid()) {
        auto params = min.UserParameters().Params();
        auto kinematics_params =
            get_kinematic_parameters(std::span<const double>(params));
        auto constrains = fcn.get_constrain(kinematics_params);
        // if converged, return result
        if (double penalty = fcn.get_parameter_penalty(kinematics_params);
            std::sqrt(square_sum_array(constrains)) < tol && penalty < 9) {
          return std::make_tuple(
              fcn.get_measured_from_parameters(kinematics_params), penalty);
        }
        for (auto &&[lag_paras, constraint] : std::views::zip(
                 std::views::chunk(params,
                                   lagrange_parmaeter_count_per_constrain),
                 constrains)) {
          auto &lambda = lag_paras[0];
          auto &mu = lag_paras[1];
          lambda += mu * constraint;
          mu *= 1.5;
        }
        for (size_t i = 0; i < total_lagrange_parameter_count; ++i) {
          upar.SetValue(i, params[i]);
        }
      } else {
        // bad fit, regenerate initial parameters
        upar = fcn.generate_initial_parameters();
      }
      if (iter > 128) {
        return std::nullopt;
      }
    }
  }
};

using KFPi0Solver = SingleKFLagMul<Problem6D>;
std::optional<std::tuple<std::array<momentum_t, 2>, double>>
kf_pi0(const std::array<momentum_t, 2> &gammas) {
  return KFPi0Solver::do_kinematics_fit(gammas);
}

class Problem3D {
public:
  static constexpr size_t kinematic_parameter_count = 3;
  static constexpr size_t constrain_count = 1;
  using answer_t = gamma_dof;
  Problem3D(double e1, double e2, double oa) : measured(e1, e2, oa) {}
  Problem3D(const answer_t &measured_) : measured(measured_) {}

  std::array<double, constrain_count>
  get_constrain(const para_view_t &params_kin) const {
    std::array<double, constrain_count> constrains;
    double m2 =
        2 * params_kin[0] * params_kin[1] * (1 - std::cos(params_kin[2]));
    constexpr double pi0_mass = 134.9768 / 1000.; // GeV/c^2
    constexpr double pi0_mass2 = pi0_mass * pi0_mass;
    constrains[0] = m2 - pi0_mass2;
    return constrains;
  }

  double get_parameter_penalty(const para_view_t &params_kin) const {
    constexpr double sigma_e_1 = 0.01734; // sigma of (smeared - true)
    constexpr double sigma_e_2 = 0.04961; // sigma of (smeared - true) / true
    constexpr double sigma_oa = 0.09274;  // sigma of (smeared - true) in rad
    double penalty = 0.0;
    penalty += chi2(params_kin[0], measured.E1, sigma_e_1);
    penalty += chi2(params_kin[1] / measured.E2, 1., sigma_e_2);
    penalty += chi2(params_kin[2], measured.theta_oa, sigma_oa);
    return penalty;
  }

  auto generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
    auto &local_rand = get_thread_local_random();
    from.Add("E1", local_rand.Gaus(measured.E1, 0.1 * measured.E1),
             0.1 * measured.E1);
    from.Add("E2", local_rand.Gaus(measured.E2, 0.1 * measured.E2),
             0.1 * measured.E2);
    from.Add("OA", local_rand.Gaus(measured.theta_oa, 0.1 * measured.theta_oa),
             0.1 * measured.theta_oa);
    return from;
  }

  static answer_t get_measured_from_parameters(const para_view_t &params_kin) {
    return gamma_dof{params_kin[0], params_kin[1], params_kin[2]};
  }

private:
  gamma_dof measured;
};

gamma_dof from_2_gamma(const std::array<momentum_t, 2> &gammas) {
  double e1 = gammas[0].E();
  double e2 = gammas[1].E();
  ROOT::Math::XYZVector v1 = gammas[0].Vect().Unit();
  ROOT::Math::XYZVector v2 = gammas[1].Vect().Unit();
  double cos_oa = v1.Dot(v2);
  cos_oa = std::clamp(cos_oa, -1.0, 1.0);
  double oa = std::acos(cos_oa);
  return gamma_dof{e1, e2, oa};
}

using KFPi0Solver3D = SingleKFLagMul<Problem3D>;
std::optional<std::tuple<gamma_dof, double>>
kf_pi0_3D(const std::array<momentum_t, 2> &gammas) {
  return KFPi0Solver3D::do_kinematics_fit(from_2_gamma(gammas));
}
