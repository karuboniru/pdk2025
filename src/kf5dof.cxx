#include "kf.h"
#include "local_rand.h"
#include "smear.h"
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <cstddef>
#include <cstdlib>
#include <optional>

#include "alm.hxx"
#include "common_tools.hxx"
#include "doublemin.hxx"

class Problem5D {
public:
  static constexpr size_t kinematic_parameter_count = 5;
  static constexpr size_t constrain_count = 2;
  static constexpr double chi2_cut = +std::numeric_limits<double>::infinity();
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
    double m_proton = 0.9382720813; // GeV/c^2
    double m_pi0 = 0.134977;        // GeV/c^2

    auto dof = get_measured_from_parameters(params_kin);
    constrains[0] = dof.calc_mass_pi0() - m_pi0;
    constrains[1] = dof.calc_mass_proton() - m_proton;

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

  // [[nodiscard]] auto
  // generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
  //   auto &rand = get_thread_local_random();
  //   from.Add("E1", measured.E1 * (1 + rand.Gaus(0, sigma_e_1)), 0.001);
  //   from.Add("E2", measured.E2 * (1 + rand.Gaus(0, sigma_e_2)), 0.001);
  //   from.Add("Elepton", measured.Elepton, 0.001);
  //   from.Add("OA", measured.theta_2gamma + rand.Gaus(0, sigma_oa), 0.001);
  //   from.Add("Theta_lepton_pi0", measured.theta_lepton_pi0 , 0.001);
  //   return from;
  // }
  [[nodiscard]] auto
  generate_initial_parameters(ROOT::Minuit2::MnUserParameters from) const {
    from.Add("E1", measured.E1);
    from.Add("E2", measured.E2);
    from.Add("Elepton", measured.Elepton, 0.001);
    from.Add("OA", measured.theta_2gamma);
    from.Add("Theta_lepton_pi0", measured.theta_lepton_pi0);
    return from;
  }

  static answer_t get_measured_from_parameters(const para_view_t &params_kin) {
    return proton_dof{params_kin[0], params_kin[1], params_kin[2],
                      params_kin[3], params_kin[4]};
  }

private:
  proton_dof measured;
};

proton_dof from_system(const std::array<momentum_t, 3> &system) {
  proton_dof ret{};
  ret.Elepton = system[0].E();
  ret.E1 = system[1].E();
  ret.E2 = system[2].E();
  ret.theta_2gamma = angle_between(system[1], system[2]);
  auto pi0 = system[1] + system[2];
  ret.theta_lepton_pi0 = angle_between(system[0], pi0);
  return ret;
}

#include <print>
using KFPi0Solver5D = SingleKFLagMul<Problem5D>;
std::optional<std::tuple<proton_dof, double>>
kf_pi0_5D_ALM(const std::array<momentum_t, 3> &system) {
  auto smear_stra = GetSmearStrategy(-11);
  auto smeared_lepton = smear_stra->do_smearing(system[0]);
  auto smeared_lepton_E = smeared_lepton.E();

  auto smear_stra_gamma = GetSmearStrategy(22);
  auto smeared_gamma_1 = smear_stra_gamma->do_smearing(system[1]);
  auto smeared_gamma_1_E = smeared_gamma_1.E();
  // auto new_lepton_p3 = system[0].Vect().Unit() * smeared_lepton_p;
  // auto new_lepton = ROOT::Math::PxPyPzMVector(
  //     new_lepton_p3.X(), new_lepton_p3.Y(), new_lepton_p3.Z(),
  //     system[0].M());
  // auto before_fit = system;
  // before_fit[0] = new_lepton;
  auto dof = from_system(system);
  // dof.E1 = smeared_gamma_1_E;
  dof.Elepton = smeared_lepton_E;
  auto res =  KFPi0Solver5D::do_kinematics_fit(dof);
  return res;
}
