#pragma once

#include <Math/Vector4Dfwd.h>

class TRandom3;
TRandom3 &get_thread_local_random();

// defines the interface for all smearing strategies
class ISmearStrategy {
public:
  ISmearStrategy() = default;
  ISmearStrategy(const ISmearStrategy &) = default;
  ISmearStrategy(ISmearStrategy &&) = default;
  ISmearStrategy &operator=(const ISmearStrategy &) = default;
  ISmearStrategy &operator=(ISmearStrategy &&) = default;

  virtual ~ISmearStrategy() = default;

  [[nodiscard]] virtual ROOT::Math::PxPyPzEVector
  do_smearing(ROOT::Math::PxPyPzEVector vec) const = 0;
};

ISmearStrategy *GetSmearStrategy(int pdg_particle);

void initializeGaussianSmearStrategy();
