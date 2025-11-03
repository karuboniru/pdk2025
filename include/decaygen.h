#pragma once

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include <memory>
#include <utility>
#include <vector>

using particle = std::pair<int, ROOT::Math::PxPyPzEVector>;

class EvtGenInterface {
public:
  static EvtGenInterface &get_instance();

  std::vector<particle> decay(
      const particle &to_decay);

private:
  EvtGenInterface();
  ~EvtGenInterface()= default;
  EvtGenInterface(const EvtGenInterface&) = delete;
  EvtGenInterface& operator=(const EvtGenInterface&) = delete;
  EvtGenInterface(EvtGenInterface&&) = delete;
  EvtGenInterface& operator=(EvtGenInterface&&) = delete;

  std::unique_ptr<EvtRandomEngine> random_engine;
  EvtGen evtgen;
};
