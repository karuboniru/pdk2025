#include "data.h"
EventRec::EventRec() : lepton{}, gamma1{}, gamma2{}, has_gamma2(false),
                   is_valid(false) {}
EventRec::EventRec(const EventRec &) = default;
EventRec &EventRec::operator=(const EventRec &) = default;
EventRec::EventRec(EventRec &&) noexcept = default;
EventRec &EventRec::operator=(EventRec &&) noexcept = default;
EventRec::~EventRec() = default;

EventRec::EventRec(momentum_t lepton, momentum_t gamma1,
                   const std::optional<momentum_t> &gamma2)
    : lepton(std::move(lepton)), gamma1(std::move(gamma1)),
      gamma2(gamma2.value_or(momentum_t{})), has_gamma2(gamma2.has_value()),
      is_valid(true) {}
