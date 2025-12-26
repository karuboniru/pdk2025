#pragma once
#include <Math/LorentzVector.h>
#include <Math/Vector4Dfwd.h>
#include <TObject.h>
#include <algorithm>
#include <optional>

using momentum_t = ROOT::Math::PxPyPzEVector;

struct EventRec {
public:
  momentum_t lepton;
  momentum_t gamma1;
  momentum_t gamma2;
  bool has_gamma2;
  bool is_valid;

  EventRec();
  EventRec(momentum_t lepton, momentum_t gamma1,
           const std::optional<momentum_t> &gamma2);
  EventRec(const EventRec &);
  EventRec &operator=(const EventRec &);
  EventRec(EventRec &&) noexcept;
  EventRec &operator=(EventRec &&) noexcept;
  virtual ~EventRec();

  ClassDef(EventRec, 2)
};

struct gamma_dof {
  double E1, E2;
  double theta_oa;
};


struct proton_dof {
  double E1, E2, Elepton;
  double theta_2gamma, theta_lepton_pi0;

  [[nodiscard]] double calc_mass_pi0() const {
    double m_pi0_sq = 2 * E1 * E2 * (1 - cos(theta_2gamma));
    return m_pi0_sq > 0 ? sqrt(m_pi0_sq) : 0.0;
  }

  // 2. Pi0 动量计算 (逻辑正确，但建议直接用能量算)
  // P_pi0^2 = E_pi0^2 - m_pi0^2
  // E_pi0 = E1 + E2
  [[nodiscard]] double calc_momentum_pi0() const {
    double m_pi0 = calc_mass_pi0();
    double E_pi0 = E1 + E2;
    double p_pi0_sq = E_pi0 * E_pi0 - m_pi0 * m_pi0;
    return p_pi0_sq > 0 ? sqrt(p_pi0_sq) : 0.0;
  }

  // 3. 质子动量计算 (修正了符号错误)
  // P_proton^2 = P_lep^2 + P_pi^2 + 2*P_lep*P_pi*cos(theta)
  [[nodiscard]] double calc_momentum_proton() const {
    double p_pi0 = calc_momentum_pi0();
    // 假设 Elepton ≈ Plepton (超相对论近似)
    double p_lepton = Elepton;

    // 修正：向量加法应使用 + 2 * ...
    double p_proton_squared = p_lepton * p_lepton + p_pi0 * p_pi0 +
                              2 * p_lepton * p_pi0 * cos(theta_lepton_pi0);

    return sqrt(p_proton_squared);
  }

  // 4. 质子质量计算 (修正了能量计算和符号错误)
  // m_p^2 = E_p^2 - P_p^2
  [[nodiscard]] double calc_mass_proton() const {
    // 修正：质子能量 = Lepton能量 + Pi0总能量 (不是静止质量!)
    double E_pi0 = E1 + E2;
    double E_proton = Elepton + E_pi0;

    // 获取修正后的动量平方
    double p_pi0 = calc_momentum_pi0();
    double p_lepton = Elepton;
    double p_proton_squared = p_lepton * p_lepton + p_pi0 * p_pi0 +
                              2 * p_lepton * p_pi0 * cos(theta_lepton_pi0);

    double m_proton_squared = E_proton * E_proton - p_proton_squared;
    return m_proton_squared > 0 ? sqrt(m_proton_squared) : 0.0;
  }
};
