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
change_vector(const ROOT::Math::PxPyPzEVector &original, double theta,
              double phi, double mag) {
  auto unit_vect_of_target = original.Vect().Unit();
  double angle_of_orig_to_z = -std::acos(unit_vect_of_target.z());
  // the rotation that rotate z axis to the target
  ROOT::Math::XYZVector rotate_vec{unit_vect_of_target.y(),
                                   -unit_vect_of_target.x(), 0};
  const ROOT::Math::AxisAngle the_rec_rot(rotate_vec, angle_of_orig_to_z);

  // consider the vector before rotate is along z-axis, and an another vector
  // next to it given by (theta,phi), the the rotation that moves z to the
  // vector before rotate will also move (theta,phi) to the desired direction
  const ROOT::Math::XYZVector to_rot{
      ROOT::Math::Polar3D<double>{original.P(), theta, phi}};
  auto new_space_vec = the_rec_rot(to_rot);

  // scale to the desired magnitude, and construct the final 4D vector
  // assuming gamma (massless) so E = |p|
  auto result = ROOT::Math::PxPyPzEVector(
      new_space_vec.x() * mag, new_space_vec.y() * mag, new_space_vec.z() * mag,
      original.P() * mag);
  return result;
}

constexpr double chi2(double value, double center, double sigma) {
  return ((value - center) * (value - center)) / (sigma * sigma);
}

class SingleKFLagMul final : public ROOT::Minuit2::FCNBase {
  static constexpr size_t lagrange_parmaeter_count = 2;

public:
  std::array<momentum_t, 2>
  get_fitted_gammas(const std::vector<double> &params) const {
    double theta1 = params[0 + lagrange_parmaeter_count];
    double phi1 = params[1 + lagrange_parmaeter_count];
    double mag1 = params[2 + lagrange_parmaeter_count];
    double theta2 = params[3 + lagrange_parmaeter_count];
    double phi2 = params[4 + lagrange_parmaeter_count];
    double mag2 = params[5 + lagrange_parmaeter_count];
    std::array<momentum_t, 2> ret{};
    ret[0] = change_vector(gammas_[0], theta1, phi1, mag1);
    ret[1] = change_vector(gammas_[1], theta2, phi2, mag2);
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
    double theta1 = params_kin[0 + lagrange_parmaeter_count];
    double mag1 = params_kin[2 + lagrange_parmaeter_count];
    double theta2 = params_kin[3 + lagrange_parmaeter_count];
    double mag2 = params_kin[5 + lagrange_parmaeter_count];
    penalty +=
        chi2(theta1, 0.0, sigma_angle[0]) + chi2(mag1, 1.0, sigma_momentum[0]);
    penalty +=
        chi2(theta2, 0.0, sigma_angle[1]) + chi2(mag2, 1.0, sigma_momentum[1]);
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
  upar.Add("theta1", rand.Gaus(0, M_PI / 120), 0.1);
  upar.Add("phi1", rand.Uniform(0, 2 * M_PI), 0.1);
  upar.Add("mag1", 1.0, 0.01);
  // gamma 2 parameters
  upar.Add("theta2", rand.Gaus(0, M_PI / 120), 0.1);
  upar.Add("phi2", rand.Uniform(0, 2 * M_PI), 0.1);
  upar.Add("mag2", 1.0, 0.01);
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
          std::abs(constraint) < tol ) {
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
