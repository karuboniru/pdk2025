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

class Problem6D {
public:
  static constexpr size_t kinematic_parameter_count = 6;
  static constexpr size_t constrain_count = 1;
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
    // auto orig_pi0_momentum = measured[0] + measured[1];
    // auto reconstructed_unit = pi0_momentum.Vect().Unit();
    // auto original_unit = orig_pi0_momentum.Vect().Unit();
    // constrains[1] =
    //     reconstructed_unit.X() - original_unit.X(); // px conservation
    // constrains[2] =
    //     reconstructed_unit.Y() - original_unit.Y(); // py conservation

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
      // penalty += chi2(angle_diff, s_ang / 1.5, s_ang);
      // penalty += rayleigh_log_likelihood_normalized(angle_diff, s_ang);
      if (angle_diff > s_ang / 1.5)
        penalty += rayleigh_log_likelihood_normalized(angle_diff, s_ang / 1.5);
      else
        penalty += rayleigh_log_likelihood_normalized(s_ang / 1.5, s_ang / 1.5);
      // momentum penalty
      double orig_mom = orig.P();
      double fitted_mom = fit.P();
      penalty += chi2(fitted_mom / orig_mom, 1., s_mom);
    }
    return penalty;
  }

  auto generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
    // auto &local_rand = get_thread_local_random();
    auto gamma_smear_strategy = GetSmearStrategy(22);
    for (auto &&[id, particle_measured] : measured | std::views::enumerate) {
      auto smeared = gamma_smear_strategy->do_smearing(particle_measured);
      auto x = smeared.Vect().X();
      auto dx = x - particle_measured.Vect().X();
      auto y = smeared.Vect().Y();
      auto dy = y - particle_measured.Vect().Y();
      auto z = smeared.Vect().Z();
      auto dz = z - particle_measured.Vect().Z();
      from.Add(std::format("x{}", id), x,
               std::abs(dx) > 1e-6 ? std::abs(dx) : 0.01);
      from.Add(std::format("y{}", id), y,
               std::abs(dy) > 1e-6 ? std::abs(dy) : 0.01);
      from.Add(std::format("z{}", id), z,
               std::abs(dz) > 1e-6 ? std::abs(dz) : 0.01);
    }
    return from;
  }

  static answer_t get_measured_from_parameters(const para_view_t &params_kin) {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error(
          "Problem6D::get_measured_from_parameters: invalid parameter size");
    }
    answer_t ret;
    for (auto &&[target, p3] :
         std::views::zip(ret, params_kin | std::views::chunk(3))) {
      double x = p3[0];
      double y = p3[1];
      double z = p3[2];
      ROOT::Math::XYZVector vec3d{x, y, z};
      double energy = vec3d.R();
      target =
          ROOT::Math::PxPyPzEVector{vec3d.x(), vec3d.y(), vec3d.z(), energy};
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

using KFPi0Solver = SingleKFLagMul<Problem6D, 0.5, 0.0, 2.0>;
std::optional<std::tuple<std::array<momentum_t, 3>, double>>
kf_pi0(const std::array<momentum_t, 3> &lepton_gammas) {
  auto gammas = std::array<momentum_t, 2>{lepton_gammas[1], lepton_gammas[2]};
  auto result = KFPi0Solver::do_kinematics_fit(gammas);
  if (!result) {
    return std::nullopt;
  }
  auto [fitted_gammas, chi2] = *result;
  return std::make_tuple(std::array<momentum_t, 3>{lepton_gammas[0],
                                                   fitted_gammas[0],
                                                   fitted_gammas[1]},
                         chi2);
}
