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
