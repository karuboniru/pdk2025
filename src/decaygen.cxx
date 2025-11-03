#include "decaygen.h"
#include "local_rand.h"

#include <mutex>

class MyRandomEngine final : public EvtRandomEngine {
public:
  double random() override {
    return get_thread_local_random().Uniform(0.0, 1.0);
  }
};

EvtGenInterface &EvtGenInterface::get_instance() {
  static EvtGenInterface instance;
  return instance;
}

EvtGenInterface::EvtGenInterface()
    : random_engine(std::make_unique<MyRandomEngine>()),
      evtgen(DATA_PATH "/decay/decay.dec",
             "/cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc15-opt/share/"
             "EvtGen/evt.pdl",
             random_engine.get()) {}

std::vector<particle> EvtGenInterface::decay(const particle &to_decay) {
  std::vector<particle> decay_products;
  auto &&[pdg, momentum] = to_decay;
  if (auto parentId = EvtPDL::evtIdFromStdHep(pdg); parentId.getId() != -1) {
    EvtVector4R p4_lab(momentum.e(), momentum.x(), momentum.y(), momentum.z());
    EvtParticle *parent = EvtParticleFactory::particleFactory(parentId, p4_lab);
    {
      static std::mutex mtx;
      std::lock_guard<std::mutex> lock(mtx);
      evtgen.generateDecay(parent);
    }
    auto push_daug = [&](this auto &&self, EvtParticle *p) {
      if (!p)
        return;
      if (p->getNDaug() == 0) {
        EvtVector4R p4 = p->getP4Lab(); // 4-momentum in lab frame
        decay_products.emplace_back(
            EvtPDL::getStdHep(p->getId()),
            ROOT::Math::PxPyPzEVector{p4.get(1), p4.get(2), p4.get(3),
                                      p4.get(0)});
      } else {
        for (int i = 0; i < p->getNDaug(); ++i)
          self(p->getDaug(i));
      }
    };
    push_daug(parent);
    parent->deleteTree();
  }
  if (decay_products.empty()) {
    decay_products.emplace_back(to_decay);
  }
  return decay_products;
}
