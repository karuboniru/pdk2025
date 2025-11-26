#pragma once

#include "event.h"


std::optional<std::tuple<std::array<momentum_t, 2>, double>>
kf_pi0(const std::array<momentum_t, 2> &gammas);
