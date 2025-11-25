#include "kf.h"
#include "local_rand.h"
#include "smear.h"
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <print>
#include <ranges>

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
  static constexpr size_t lagrange_parmaeter_count = 2;

public:
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

  double get_constrain(const std::vector<double> &params_kin) const {
    return get_pi0_mass(params_kin) - 0.134977;
  }

  double get_parameter_penalty(const std::vector<double> &params_kin) const {
    double penalty = 0.0;
    auto fitted_gammas = get_fitted_gammas(params_kin);
    for (size_t i = 0; i < 2; ++i) {
      auto &original = gammas_[i];
      auto &fitted = fitted_gammas[i];
      // angle penalty
      ROOT::Math::XYZVector orig_dir = original.Vect().Unit();
      ROOT::Math::XYZVector fitted_dir = fitted.Vect().Unit();
      double angle_diff = std::acos(orig_dir.Dot(fitted_dir));
      penalty += chi2(angle_diff, 0.0, sigma_angle[i]);
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
    double lambda = params[0];
    double mu = params[1];
    double constraint = get_constrain(params);
    double penalty = get_parameter_penalty(params);
    double fcn_value =
        penalty + lambda * constraint + mu * constraint * constraint;
    return fcn_value;
  }

  double Up() const override { return 1.0; }

private:
  std::array<momentum_t, 2> gammas_;
  std::array<double, 2> sigma_angle{};
  std::array<double, 2> sigma_momentum{};
};

ROOT::Minuit2::MnUserParameters construct_initial_parameters() {
  auto &rand = get_thread_local_random();
  ROOT::Minuit2::MnUserParameters upar;
  // lagrange multipliers
  upar.Add("lambda", 0.0); // not being fit but iterated
  upar.Fix("lambda");
  upar.Add("mu", 1); // not being fit but iterated
  upar.Fix("mu");
  // gamma 1 parameters
  upar.Add("dx1", rand.Gaus(0,0.1), 0.5);
  upar.Add("dy1", rand.Gaus(0,0.1), 0.5);
  upar.Add("dz1", rand.Gaus(0,0.1), 0.5);
  // gamma 2 parameters
  upar.Add("dx2", rand.Gaus(0,0.1), 0.5);
  upar.Add("dy2", rand.Gaus(0,0.1), 0.5);
  upar.Add("dz2", rand.Gaus(0,0.1), 0.5);
  return upar;
}

std::optional<std::array<momentum_t, 2>>
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
      double constraint = fcn.get_constrain(params);
      if (double penalty = fcn.get_parameter_penalty(params);
          std::abs(constraint) < tol && penalty < 9) {
        return fcn.get_fitted_gammas(params);
      }
      // update lagrange multipliers
      double lambda = params[0];
      double mu = params[1];
      lambda += 2 * mu * constraint;
      mu *= 2;
      upar.SetValue("lambda", lambda);
      upar.SetValue("mu", mu);
    } else {
      // re-initialize parameters and try again
      upar = construct_initial_parameters();
    }
    if (iter > 128) {
      std::println("KF did not converge within {} iterations, died at lambda = "
                   "{}, mu = {}",
                   iter, upar.Value("lambda"), upar.Value("mu"));
      return std::nullopt;
    }
  }
}
