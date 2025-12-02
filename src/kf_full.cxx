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
#include <ranges>
#include <utility>

#include "alm.hxx"

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

class Problem6D {
public:
  static constexpr size_t kinematic_parameter_count = 6;
  static constexpr size_t constrain_count = 3;
  static constexpr double chi2_cut = 16;
  using answer_t = std::array<momentum_t, 2>;
  Problem6D(answer_t measured_) : measured(std::move(measured_)) {
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
    std::array<double, constrain_count> constrains{};
    auto fitted = get_measured_from_parameters(params_kin);
    auto pi0_momentum = fitted[0] + fitted[1];
    constrains[0] = pi0_momentum.M() - 0.134977;
    // reconstructed pi0 follows preserved momentum
    auto orig_pi0_momentum = measured[0] + measured[1];
    auto reconstructed_unit = pi0_momentum.Vect().Unit();
    auto original_unit = orig_pi0_momentum.Vect().Unit();
    constrains[1] =
        reconstructed_unit.X() - original_unit.X(); // px conservation
    constrains[2] =
        reconstructed_unit.Y() - original_unit.Y(); // py conservation

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
      // penalty += rayleigh_log_likelihood_normalized(angle_diff, s_ang);

      // momentum penalty
      double orig_mom = orig.P();
      double fitted_mom = fit.P();
      penalty += chi2(fitted_mom / orig_mom, 1., s_mom);
    }
    return penalty;
  }

  static auto
  generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) {
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

using KFPi0Solver = SingleKFLagMul<Problem6D, 0.5, 0.0, 2.0>;
std::optional<std::tuple<std::array<momentum_t, 2>, double>>
kf_pi0(const std::array<momentum_t, 2> &gammas) {
  return KFPi0Solver::do_kinematics_fit(gammas);
}
