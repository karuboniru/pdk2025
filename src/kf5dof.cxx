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

#include "alm.hxx"
#include "common_tools.hxx"
#include "doublemin.hxx"

class Problem5D {
public:
  static constexpr size_t kinematic_parameter_count = 5;
  static constexpr size_t constrain_count = 1;
  static constexpr double chi2_cut = 16;
  using answer_t = proton_dof;
  Problem5D(double e1, double e2, double oa) : measured(e1, e2, oa) {}
  Problem5D(const answer_t &measured_) : measured(measured_) {}

  static std::array<double, constrain_count>
  get_constrain(const para_view_t &params_kin) {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error(
          "Problem5D::get_constrain: invalid parameter size");
    }
    std::array<double, constrain_count> constrains{};
    double m2 =
        2 * params_kin[0] * params_kin[1] * (1 - std::cos(params_kin[2]));
    constexpr double pi0_mass = 134.9768 / 1000.; // GeV/c^2
    constexpr double pi0_mass2 = pi0_mass * pi0_mass;
    constrains[0] = m2 - pi0_mass2;
    return constrains;
  }

  static constexpr double sigma_e_1 = 0.01734; // sigma of (smeared - true)
  static constexpr double sigma_e_2 =
      0.07625; // sigma of (smeared - true) / true
  static constexpr double sigma_e_lepton = 0.0167;
  static constexpr double sigma_oa = 0.1066; // sigma of (smeared - true) in rad
  static constexpr double sigma_pi0_lepton =
      0.045; // sigma of (smeared - true) in rad

  [[nodiscard]] double
  get_parameter_penalty(const para_view_t &params_kin) const {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error(
          "Problem5D::get_parameter_penalty: invalid parameter size");
    }

    double penalty = 0.0;
    penalty += chi2(params_kin[0], measured.E1, sigma_e_1);
    penalty += chi2(params_kin[1], measured.E2, sigma_e_2 * measured.E2);
    penalty += chi2(params_kin[2], measured.Elepton, sigma_e_lepton);
    penalty += chi2(params_kin[3], measured.theta_2gamma, sigma_oa);
    penalty += chi2(params_kin[4], measured.theta_lepton_pi0, sigma_pi0_lepton);
    return penalty;
  }

  [[nodiscard]] auto
  generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
    auto &rand = get_thread_local_random();
    from.Add("E1", measured.E1 * (1 + rand.Gaus(0, sigma_e_1)), 0.001);
    from.Add("E2", measured.E2 * (1 + rand.Gaus(0, sigma_e_2)), 0.001);
    from.Add("Elepton", measured.Elepton, 0.001);
    from.Add("OA", measured.theta_2gamma + rand.Gaus(0, sigma_oa), 0.001);
    from.Add("Theta_lepton_pi0", measured.theta_lepton_pi0, 0.001);
    return from;
  }

  static answer_t get_measured_from_parameters(const para_view_t &params_kin) {
    return proton_dof{params_kin[0], params_kin[1], params_kin[2],
                      params_kin[3], params_kin[4]};
  }

private:
  proton_dof measured;
};

proton_dof from_2_gamma(const std::array<momentum_t, 3> &system) {
  proton_dof ret{};
  ret.Elepton = system[0].E();
  ret.E1 = system[1].E();
  ret.E2 = system[2].E();
  ret.theta_2gamma = angle_between(system[1], system[2]);
  auto pi0 = system[1] + system[2];
  ret.theta_lepton_pi0 = angle_between(system[0], pi0);
  return ret;
}

using KFPi0Solver3DDM = SingleKFLagMulDoubleMin<Problem5D>;
std::optional<std::tuple<proton_dof, double>>
kf_pi0_3D_DM(const std::array<momentum_t, 3> &system) {}
