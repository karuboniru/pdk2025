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
  return std::log(value / (sigma * sigma)) -
         (value * value) / (2 * sigma * sigma);
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

class SingleKFLagMul final : public ROOT::Minuit2::FCNBase {
public:
  static constexpr size_t lagrange_parmaeter_count_per_constrain = 2;
  static constexpr size_t constrain_count = 3;
  static constexpr size_t lagrange_parmaeter_count =
      lagrange_parmaeter_count_per_constrain * constrain_count;
  static constexpr size_t motion_particle_count = 2;

  std::array<momentum_t, 2>
  get_fitted_gammas(const std::vector<double> &params) const {
    std::array<momentum_t, 2> ret{};
    for (size_t i = 0; i < 2; ++i) {
      double dx = params[lagrange_parmaeter_count + i * 3 + 0];
      double dy = params[lagrange_parmaeter_count + i * 3 + 1];
      double dz = params[lagrange_parmaeter_count + i * 3 + 2];
      ret[i] = change_vector(gammas_[i], dx, dy, dz);
    }
    return ret;
  }

  double get_pi0_mass(const std::vector<double> &params_kin) const {
    auto fitted_gammas = get_fitted_gammas(params_kin);
    auto pi0_momentum = fitted_gammas[0] + fitted_gammas[1];
    return pi0_momentum.M();
  }

  std::array<double, constrain_count>
  get_constrain(const std::vector<double> &params_kin) const {
    // return get_pi0_mass(params_kin) - 0.134977;
    double constrain_0 = get_pi0_mass(params_kin) - 0.134977;
    // return {constrain_0};
    auto fitted_gammas = get_fitted_gammas(params_kin);
    auto fitted_pi0_momentum =
        (fitted_gammas[0] + fitted_gammas[1]).Vect().Unit();
    auto raw_pi0_momentum = (gammas_[0] + gammas_[1]).Vect().Unit();
    double constrain_1 = fitted_pi0_momentum.x() - raw_pi0_momentum.x();
    double constrain_2 = fitted_pi0_momentum.y() - raw_pi0_momentum.y();
    // double constrain_3 = fitted_pi0_momentum.z() - raw_pi0_momentum.z();
    return {constrain_0, constrain_1, constrain_2};
  }

  double get_parameter_penalty(const std::vector<double> &params_kin) const {
    double penalty = 0.0;
    auto fitted_gammas = get_fitted_gammas(params_kin);
    for (size_t i = 0; i < motion_particle_count; ++i) {
      auto &original = gammas_[i];
      auto &fitted = fitted_gammas[i];
      // angle penalty
      ROOT::Math::XYZVector orig_dir = original.Vect().Unit();
      ROOT::Math::XYZVector fitted_dir = fitted.Vect().Unit();
      double angle_diff = std::acos(orig_dir.Dot(fitted_dir));
      // penalty += chi2(angle_diff, 0.0, sigma_angle[i]);
      penalty += rayleigh_log_likelihood_normalized(angle_diff,
      sigma_angle[i]);

      // momentum penalty
      double orig_mom = original.P();
      double fitted_mom = fitted.P();
      penalty += chi2(fitted_mom / orig_mom, 1., sigma_momentum[i]);
    }
    return penalty;
  }

  SingleKFLagMul(const std::array<momentum_t, 2> &gammas) : gammas_(gammas) {
    auto gamma_smear_strategy = GetSmearStrategy(22);
    if (!gamma_smear_strategy) {
      throw std::runtime_error("No smear strategy for gamma (PDG 22)");
    }
    for (auto &&[p4, sigma_ang, sigma_mom] :
         std::views::zip(gammas_, sigma_angle, sigma_momentum)) {
      double e = p4.E();
      sigma_ang = gamma_smear_strategy->get_sigma_angle(e) * M_PI / 180.0;
      sigma_mom = gamma_smear_strategy->get_sigma_energy(e);
    }
  }

  double operator()(const std::vector<double> &params) const override {
    // total parameters:
    // lagrange multipliers: lagrange_parmaeter_count
    // kinematic parameters: 2 * 3 (dx, dy, dz for each gamma)
    if (params.size() !=
        lagrange_parmaeter_count + 2 * 3 /* 2 gammas, each with dx,dy,dz */) {
      throw std::runtime_error("SingleKFLagMul: invalid parameter size");
    }
    auto constraint = get_constrain(params);
    double penalty = get_parameter_penalty(params);
    for (auto &&[p, c] : std::views::zip(
             params | std::views::chunk(lagrange_parmaeter_count_per_constrain),
             constraint)) {
      auto lambda = p[0];
      auto mu = p[1];
      penalty += lambda * c + mu * c * c;
    }
    return penalty;
  }

  double Up() const override { return 1.0; }

private:
  std::array<momentum_t, motion_particle_count> gammas_;
  std::array<double, motion_particle_count> sigma_angle{};
  std::array<double, motion_particle_count> sigma_momentum{};
};

ROOT::Minuit2::MnUserParameters construct_initial_parameters() {
  auto &rand = get_thread_local_random();
  ROOT::Minuit2::MnUserParameters upar;

  for (size_t i = 0; i < SingleKFLagMul::constrain_count; ++i) {
    upar.Add(std::format("lambda{}", i), 0.0);
    upar.Fix(std::format("lambda{}", i));
    upar.Add(std::format("mu{}", i), 0.1);
    upar.Fix(std::format("mu{}", i));
  }
  // gamma 1 parameters
  for (size_t i = 0; i < SingleKFLagMul::motion_particle_count; ++i) {
    upar.Add(std::format("dx{}", i), rand.Gaus(0, 0.1), 0.5);
    upar.Add(std::format("dy{}", i), rand.Gaus(0, 0.1), 0.5);
    upar.Add(std::format("dz{}", i), rand.Gaus(0, 0.1), 0.5);
  }

  return upar;
}

template <size_t N> double norm_array(const std::array<double, N> &arr) {
  double sum = 0.0;
  for (const auto &v : arr) {
    sum += v * v;
  }
  return std::sqrt(sum);
}

std::optional<std::tuple<std::array<momentum_t, 2>, double>>
kf_pi0(const std::array<momentum_t, 2> &gammas) {
  auto &rand = get_thread_local_random();
  // double sigma_angle = 0.01;     // in radian
  // double sigma_magnitude = 0.01; // fractional
  SingleKFLagMul fcn(gammas);
  auto upar = construct_initial_parameters();
  constexpr double tol = 1e-5;
  for (size_t iter{};; iter++) {
    ROOT::Minuit2::MnMigrad migrad(fcn, upar, ROOT::Minuit2::MnStrategy{2});
    ROOT::Minuit2::FunctionMinimum min = migrad();
    if (min.IsValid()) {
      auto params = min.UserParameters().Params();
      auto constraint = fcn.get_constrain(params);
      if (double penalty = fcn.get_parameter_penalty(params);
          norm_array(constraint) < tol && penalty < 9) {
        return std::make_tuple(fcn.get_fitted_gammas(params), penalty);
      }
      for (size_t id{}; id < SingleKFLagMul::constrain_count; id++) {
        auto lambdapos =
            (id * SingleKFLagMul::lagrange_parmaeter_count_per_constrain) + 0;
        auto mupos =
            (id * SingleKFLagMul::lagrange_parmaeter_count_per_constrain) + 1;
        double lambda = params[lambdapos];
        double mu = params[mupos];
        lambda += 2 * mu * constraint[id];
        mu *= 2;
        upar.SetValue(lambdapos, lambda);
        upar.SetValue(mupos, mu);
      }
    } else {
      // re-initialize parameters and try again
      upar = construct_initial_parameters();
    }
    if (iter > 128) {
      // std::println("KF did not converge within {} iterations, died at lambda
      // = "
      //              "{}, mu = {}",
      //              iter, upar.Value("lambda"), upar.Value("mu"));
      return std::nullopt;
    }
  }
}
