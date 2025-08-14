#include "smear.h"
#include "local_rand.h"
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Rtypes.h>
#include <TRandom.h>
#include <TVector3.h>
#include <cmath>
#include <map>
#include <memory>

static std::map<int, std::unique_ptr<ISmearStrategy>> smear_strategies;

// smear: the input PxPyPzEVector to be smeared
// angle: the angle in degrees to smear the vector by
ROOT::Math::PxPyPzEVector smear_angle(const ROOT::Math::PxPyPzEVector &original,
                                      double angle) {
  double angle_rad = angle * M_PI / 180.0;
  auto phi = get_thread_local_random().Uniform(0, 2 * M_PI);

  // ROOT::Math::AxisAngle rotate_back
  auto unit_vect_of_target = original.Vect().Unit();

  // rotate_vec is normalized dot product
  // of space_of_raw and z axis, with angle between them
  // which does the job of rotating the z axis to space_of_raw
  // thus adding the extra smearing around the axis
  // to the smearing around space_of_raw
  double angle_of_orig_to_z = -std::acos(unit_vect_of_target.z());
  ROOT::Math::XYZVector rotate_vec{unit_vect_of_target.y(),
                                   -unit_vect_of_target.x(), 0};
  const ROOT::Math::AxisAngle the_rec_rot(rotate_vec, angle_of_orig_to_z);
  // const ROOT::Math::Polar3D<double> to_rot{original.P(), angle_rad, phi};
  const ROOT::Math::XYZVector to_rot{
      ROOT::Math::Polar3D<double>{original.P(), angle_rad, phi}};
  auto new_space_vec = the_rec_rot(to_rot);
  auto result = ROOT::Math::PxPyPzEVector(new_space_vec.x(), new_space_vec.y(),
                                          new_space_vec.z(), original.E());

  return result;
}

// scale momentum but left mass & direction unchanged
ROOT::Math::PxPyPzEVector
smear_momentum(const ROOT::Math::PxPyPzEVector &original, double scale) {
  auto p3_original = original.Vect();
  auto p3_new = p3_original * scale;
  return ROOT::Math::PxPyPzEVector{ROOT::Math::PxPyPzMVector(
      p3_new.x(), p3_new.y(), p3_new.z(), original.M())};
}

class SmearGaussianDirection : public ISmearStrategy {
public:
  SmearGaussianDirection(double sigma) : m_sigma(sigma) {}

  SmearGaussianDirection(const SmearGaussianDirection &) = default;
  SmearGaussianDirection(SmearGaussianDirection &&) = default;
  SmearGaussianDirection &operator=(const SmearGaussianDirection &) = default;
  SmearGaussianDirection &operator=(SmearGaussianDirection &&) = default;

  virtual ~SmearGaussianDirection() = default;

  ROOT::Math::PxPyPzEVector do_smearing(ROOT::Math::PxPyPzEVector vec) const {
    double angle = get_thread_local_random().Gaus(0, m_sigma);
    return smear_angle(vec, angle);
  }

private:
  double m_sigma;
};

class SmearRayleighDirection : public ISmearStrategy {
public:
  SmearRayleighDirection(double sigma) : m_sigma(sigma) {}

  SmearRayleighDirection(const SmearRayleighDirection &) = default;
  SmearRayleighDirection(SmearRayleighDirection &&) = default;
  SmearRayleighDirection &operator=(const SmearRayleighDirection &) = default;
  SmearRayleighDirection &operator=(SmearRayleighDirection &&) = default;

  virtual ~SmearRayleighDirection() = default;

  ROOT::Math::PxPyPzEVector do_smearing(ROOT::Math::PxPyPzEVector vec) const {
    double cdf = get_thread_local_random().Uniform(0, 1);
    double angle = m_sigma * std::sqrt(-2.0 * std::log(cdf));
    return smear_angle(vec, angle);
  }

private:
  double m_sigma;
};

template <class AngularSmear> class SmearMomentumAndDir : public AngularSmear {
public:
  template <typename... Args>
  SmearMomentumAndDir(double momentum_frac, Args &&...args)
      : AngularSmear(std::forward<Args>(args)...),
        m_momentum_frac(momentum_frac) {}
  SmearMomentumAndDir(const SmearMomentumAndDir &) = default;
  SmearMomentumAndDir(SmearMomentumAndDir &&) = default;
  SmearMomentumAndDir &operator=(const SmearMomentumAndDir &) = default;
  SmearMomentumAndDir &operator=(SmearMomentumAndDir &&) = default;
  virtual ~SmearMomentumAndDir() = default;

  ROOT::Math::PxPyPzEVector do_smearing(ROOT::Math::PxPyPzEVector vec) const {
    auto smeared_vec = AngularSmear::do_smearing(vec);
    double momentum_scaling =
        get_thread_local_random().Gaus(1, m_momentum_frac);
    return smear_momentum(smeared_vec, momentum_scaling);
  }

private:
  double m_momentum_frac;
};

ISmearStrategy *GetSmearStrategy(int pdg_particle) {
  auto it = smear_strategies.find(pdg_particle);
  if (it != smear_strategies.end()) {
    return it->second.get();
  }
  return nullptr;
}

void initializeGaussianSmearStrategy() {
  smear_strategies[11] =
      std::make_unique<SmearMomentumAndDir<SmearRayleighDirection>>(0.06, 2.91);
  smear_strategies[-11] =
      std::make_unique<SmearMomentumAndDir<SmearRayleighDirection>>(0.06, 2.91);
  smear_strategies[13] =
      std::make_unique<SmearMomentumAndDir<SmearRayleighDirection>>(0.02, 1.75);
  smear_strategies[22] =
      std::make_unique<SmearMomentumAndDir<SmearRayleighDirection>>(0.06, 1.75);
}
