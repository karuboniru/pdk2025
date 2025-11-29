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
#include "doublemin.hxx"

class Problem3D {
public:
  static constexpr size_t kinematic_parameter_count = 3;
  static constexpr size_t constrain_count = 1;
  static constexpr double chi2_cut = 16;
  using answer_t = gamma_dof;
  Problem3D(double e1, double e2, double oa) : measured(e1, e2, oa) {}
  Problem3D(const answer_t &measured_) : measured(measured_) {}

  static std::array<double, constrain_count>
  get_constrain(const para_view_t &params_kin) {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error(
          "Problem3D::get_constrain: invalid parameter size");
    }
    std::array<double, constrain_count> constrains{};
    double m2 =
        2 * params_kin[0] * params_kin[1] * (1 - std::cos(params_kin[2]));
    constexpr double pi0_mass = 134.9768 / 1000.; // GeV/c^2
    constexpr double pi0_mass2 = pi0_mass * pi0_mass;
    constrains[0] = m2 - pi0_mass2;
    return constrains;
  }

  [[nodiscard]] double
  get_parameter_penalty(const para_view_t &params_kin) const {
    if (params_kin.size() != kinematic_parameter_count) {
      throw std::runtime_error(
          "Problem3D::get_parameter_penalty: invalid parameter size");
    }
    constexpr double sigma_e_1 = 0.01734; // sigma of (smeared - true)
    constexpr double sigma_e_2 = 0.07625; // sigma of (smeared - true) / true
    constexpr double sigma_oa = 0.09274;  // sigma of (smeared - true) in rad
    double penalty = 0.0;
    penalty += chi2(params_kin[0], measured.E1, sigma_e_1);
    penalty += chi2(params_kin[1] / measured.E2, 1., sigma_e_2);
    // penalty += chi2(params_kin[2], measured.theta_oa, sigma_oa);
    return penalty;
  }

  [[nodiscard]] auto
  generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
    auto &rand = get_thread_local_random();
    from.Add("E1", rand.Gaus(measured.E1, 0.01), 0.01734);
    from.Add("E2", rand.Gaus(measured.E2, 0.04961 * measured.E2),
             0.04961 * measured.E2);
    // from.Add("E2", measured.E2);
    from.Add("OA", measured.theta_oa);
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
kf_pi0_3D_1(const std::array<momentum_t, 2> &gammas) {
  return std::nullopt;
  return KFPi0Solver3D::do_kinematics_fit(from_2_gamma(gammas));
}

using KFPi0Solver3DDM = SingleKFLagMulDoubleMin<Problem3D>;
std::optional<std::tuple<gamma_dof, double>>
kf_pi0_3D(const std::array<momentum_t, 2> &gammas) {
  // ScopedTimer timer("KF Pi0 3D DM");
  return std::nullopt;
  return KFPi0Solver3DDM::do_kinematics_fit(from_2_gamma(gammas));
}
