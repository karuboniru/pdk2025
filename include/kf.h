#pragma once

#include "event.h"
#include "data.h"

std::optional<std::tuple<std::array<momentum_t, 2>, double>>
kf_pi0(const std::array<momentum_t, 2> &gammas);


std::optional<std::tuple<gamma_dof, double>>
kf_pi0_3D(const std::array<momentum_t, 2> &gammas);

std::optional<std::tuple<gamma_dof, double>>
kf_pi0_3D_1(const std::array<momentum_t, 2> &gammas);
