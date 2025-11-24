#pragma once

#include "event.h"


std::optional<std::array<momentum_t, 2>>
kf_pi0(const std::array<momentum_t, 2> &gammas);
