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

class Problem9D {
public:
  static constexpr size_t particle_count = 3;         // lepton, gamma1, gamma2
  static constexpr size_t parameter_per_particle = 3; // x,y,z shifts
  static constexpr size_t kinematic_parameter_count =
      particle_count * parameter_per_particle;
  static constexpr size_t constrain_count = 2;
  static constexpr double chi2_cut = +std::numeric_limits<double>::infinity();
  using answer_t = std::array<momentum_t, 3>; // lepton, gamma1, gamma2
  Problem9D(answer_t measured_) : measured(std::move(measured_)) {
    auto gamma_smear_strategy = GetSmearStrategy(22);
    auto positron_smear_strategy = GetSmearStrategy(-11);
    if (!gamma_smear_strategy) {
      throw std::runtime_error("No smear strategy for gamma (PDG 22)");
    }
    if (!positron_smear_strategy) {
      throw std::runtime_error("No smear strategy for positron (PDG -11)");
    }
    // positron
    sigma_angle[0] = positron_smear_strategy->get_sigma_angle(measured[0].E()) *
                     M_PI / 180.0;
    sigma_momentum[0] =
        positron_smear_strategy->get_sigma_energy(measured[0].E());
    // gammas
    for (size_t i = 1; i < 3; ++i) {
      double e = measured[i].E();
      sigma_angle[i] = gamma_smear_strategy->get_sigma_angle(e) * M_PI / 180.0;
      sigma_momentum[i] = gamma_smear_strategy->get_sigma_energy(e);
    }
  }

  static std::array<double, constrain_count>
  get_constrain(const para_view_t &params_kin) {
    std::array<double, constrain_count> constrains{};
    auto fitted = get_measured_from_parameters(params_kin);
    auto pi0_momentum = fitted[1] + fitted[2];
    // pi0 on -> mass constraint
    const double pi0_mass = 0.134977; // GeV/c^2
    constrains[0] = pi0_momentum.M() - pi0_mass;
    // reconstructed pi0 follows preserved momentum
    auto epi_system = fitted[0] + pi0_momentum;
    // constrains[1] = epi_system.Px() - 0.0;
    // constrains[2] = epi_system.Py() - 0.0;
    // constrains[3] = epi_system.Pz() - 0.0;
    // epi on -> proton on shell
    double m_proton = 0.938272;
    constrains[1] = epi_system.M() - m_proton;
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
    for (auto &&[id, particle_measured] : measured | std::views::enumerate) {
      from.Add(std::format("x{}", id), local_rand.Gaus(0.0, 0.01), 0.01);
      from.Add(std::format("y{}", id), local_rand.Gaus(0.0, 0.01), 0.01);
      from.Add(std::format("z{}", id), local_rand.Gaus(0.0, 0.01), 0.01);
    }
    return from;
  }

  static answer_t get_measured_from_parameters(const para_view_t &params_kin) {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error("Invalid parameter size");
    }
    answer_t ret{};
    for (auto &&[p3, mass, target] :
         std::views::zip(params_kin | std::views::chunk(parameter_per_particle),
                         masses, ret)) {
      auto E = std::sqrt(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2] +
                         mass * mass);
      target = momentum_t(p3[0], p3[1], p3[2], E);
    }

    return ret;
  }

private:
  answer_t measured;
  std::array<double, particle_count> sigma_angle{};
  std::array<double, particle_count> sigma_momentum{};
  constexpr static std::array<double, particle_count> masses{0.000511, 0.0,
                                                             0.0};
};

using KFPi0Solver = SingleKFLagMul<Problem9D, 0.5, 0.0, 2.0>;
std::optional<std::tuple<std::array<momentum_t, 3>, double>>
kf_pi0_full(const std::array<momentum_t, 3> &system) {
  return KFPi0Solver::do_kinematics_fit(system);
}
