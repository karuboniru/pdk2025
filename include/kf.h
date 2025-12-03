#pragma once

#include "data.h"
#include "event.h"

std::optional<std::tuple<std::array<momentum_t, 2>, double>>
kf_pi0(const std::array<momentum_t, 2> &gammas);
std::optional<std::tuple<std::array<momentum_t, 3>, double>>
kf_pi0(const std::array<momentum_t, 3> &lepton_gammas);

std::optional<std::tuple<gamma_dof, double>>
kf_pi0_3D_DM(const std::array<momentum_t, 2> &gammas);

std::optional<std::tuple<gamma_dof, double>>
kf_pi0_3D_ALM(const std::array<momentum_t, 2> &gammas);

std::optional<std::tuple<std::array<momentum_t, 3>, double>>
kf_pi0_full(const std::array<momentum_t, 3> &system);
