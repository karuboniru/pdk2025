#include "kf.h"
#include "local_rand.h"
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

// vary lambda, for each lambda do kinematic fit
// until find the best lambda that minimize the total cost
class SingleKFLagMul final : public ROOT::Minuit2::FCNBase {
public:
  // real fit to kinematic with given fixed lambda
  class SingleKFLagMulKinematic final : public ROOT::Minuit2::FCNBase {
  public:
    SingleKFLagMulKinematic(const SingleKFLagMul *parent, double lambda)
        : parent_(parent), lambda_(lambda) {}
    std::array<momentum_t, 2>
    get_fitted_gammas(const std::vector<double> &params) const {
      double theta1 = params[0];
      double phi1 = params[1];
      double mag1 = params[2];
      double theta2 = params[3];
      double phi2 = params[4];
      double mag2 = params[5];
      std::array<momentum_t, 2> ret{};
      ret[0] = change_vector(parent_->gammas_[0], theta1, phi1, mag1);
      ret[1] = change_vector(parent_->gammas_[1], theta2, phi2, mag2);
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
      double theta1 = params_kin[0];
      double mag1 = params_kin[2];
      double theta2 = params_kin[3];
      double mag2 = params_kin[5];
      double sigma_ang_ = parent_->sigma_ang_;
      double sigma_mag_ = parent_->sigma_mag_;
      // phi don't have penalty
      penalty += (theta1 * theta1) / (sigma_ang_ * sigma_ang_);
      penalty += (theta2 * theta2) / (sigma_ang_ * sigma_ang_);
      penalty += ((mag1 - 1.0) * (mag1 - 1.0)) / (sigma_mag_ * sigma_mag_);
      penalty += ((mag2 - 1.0) * (mag2 - 1.0)) / (sigma_mag_ * sigma_mag_);
      return penalty;
    }

    double operator()(const std::vector<double> &params) const override {
      double constrain = get_constrain(params);
      double penalty = get_parameter_penalty(params);

      return (lambda_ * constrain) + penalty;
    }

    double Up() const override { return 1.0; }

  private:
    const SingleKFLagMul *parent_;
    double lambda_;
  };

  SingleKFLagMul(const std::array<momentum_t, 2> &gammas, double sigma_ang,
                 double sigma_mag)
      : gammas_(gammas), sigma_ang_(sigma_ang), sigma_mag_(sigma_mag) {}

  auto do_minimization(double lambda) const {
    SingleKFLagMulKinematic kinematic_fcn(this, lambda);
    auto &rand = get_thread_local_random();
    for (size_t attempt = 0;; ++attempt) {
      ROOT::Minuit2::MnUserParameters upar_kin;
      upar_kin.Add("theta1", sigma_ang_ / 2 * rand.Gaus(1., 0.2), sigma_ang_, 0,
                   M_PI);
      upar_kin.Add("phi1", 0.0 + rand.Gaus(0, M_PI), 0.3, -M_PI, M_PI);
      upar_kin.Add("mag1", 1.0 * rand.Gaus(1., 0.1), sigma_mag_ / 2, 0.5, 2.0);
      upar_kin.Add("theta2", sigma_ang_ / 2 * rand.Gaus(1., 0.2), sigma_ang_, 0,
                   M_PI);
      upar_kin.Add("phi2", 0.0 + rand.Gaus(0, M_PI), 0.3, -M_PI, M_PI);
      upar_kin.Add("mag2", 1.0 * rand.Gaus(1., 0.1), sigma_mag_ / 2, 0.5, 2.0);
      ROOT::Minuit2::MnMigrad migrad_kin(kinematic_fcn, upar_kin);
      ROOT::Minuit2::FunctionMinimum min_kin = migrad_kin(1000000);
      if (min_kin.IsValid() || attempt == 16) {
        if (attempt == 16 && !min_kin.IsValid()) {
          std::println("kin: KF kinematic fit failed to converge @ lambda = {}",
                       lambda);
          if (min_kin.HesseFailed()) {
            std::println("kin: Hesse failed!");
          }
          if (min_kin.HasReachedCallLimit()) {
            std::println("kin: Reached call limit!");
          }
          if (min_kin.IsAboveMaxEdm()) {
            std::println("kin: Above max edm!");
          }
        }
        return min_kin;
      }
    }
  }

  double operator()(const std::vector<double> &params) const override {
    double lambda = params[0];
    auto min_kin = do_minimization(lambda);
    return min_kin.Fval();
  }

  double Up() const override { return 1.0; }

private:
  std::array<momentum_t, 2> gammas_;
  double sigma_ang_;
  double sigma_mag_;
};

std::optional<std::array<momentum_t, 2>>
kf_pi0(const std::array<momentum_t, 2> &gammas) {
  auto &rand = get_thread_local_random();
  double sigma_angle = 1;    // in radian
  double sigma_magnitude = 1; // fractional
  SingleKFLagMul fcn(gammas, sigma_angle, sigma_magnitude);
  for (size_t attempt{};; attempt++) {
    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("lambda", rand.Gaus(0, 1000), 10, -1e5, 1e5);
    ROOT::Minuit2::MnSimplex migrad(fcn, upar);
    ROOT::Minuit2::FunctionMinimum min = migrad(4000, 1.);
    if (!min.IsValid()) {
      if (attempt > 16) {
        if (min.HesseFailed()) {
          std::println("lambda: Hesse failed!");
        }
        if (min.HasReachedCallLimit()) {
          std::println("lambda: Reached call limit!");
        }
        if (min.IsAboveMaxEdm()) {
          std::println("lambda: Above max edm!");
        }
        return std::nullopt;
      }
      continue;
    }
    const auto &result_params = min.UserParameters().Params();
    double lambda = result_params[0];

    SingleKFLagMul::SingleKFLagMulKinematic kinematic_fcn(&fcn, lambda);
    auto min_kin = fcn.do_minimization(lambda);
    if (!min_kin.IsValid()) {
      if (attempt > 16) {
        std::println("Final kinematic fit failed, should not happen");
        return std::nullopt;
      }
      continue;
    }
    if (auto diff =
            kinematic_fcn.get_constrain(min_kin.UserParameters().Params());
        std::abs(diff) > 0.0001) {
      if (attempt > 16) {
        std::println(
            "Final kinematic fit constraint not satisfied, bad fit, diff = {}",
            diff);
        return std::nullopt;
      }
      continue;
    }
    const auto &kin_params = min_kin.UserParameters().Params();
    auto fitted_gammas = kinematic_fcn.get_fitted_gammas(kin_params);
    return fitted_gammas;
  }
}
