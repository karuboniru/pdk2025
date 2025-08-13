#include "smear.h"
#include "local_rand.h"
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
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

class SmearGaussianDirection final : public ISmearStrategy {
public:
  SmearGaussianDirection(double sigma) : m_sigma(sigma) {}

  SmearGaussianDirection(const SmearGaussianDirection &) = default;
  SmearGaussianDirection(SmearGaussianDirection &&) = default;
  SmearGaussianDirection &operator=(const SmearGaussianDirection &) = default;
  SmearGaussianDirection &operator=(SmearGaussianDirection &&) = default;

  virtual ~SmearGaussianDirection() final = default;

  ROOT::Math::PxPyPzEVector
  do_smearing(ROOT::Math::PxPyPzEVector vec) const final {
    double angle = get_thread_local_random().Gaus(0, m_sigma);
    return smear_angle(vec, angle);
  }

private:
  double m_sigma;
};

ISmearStrategy *GetSmearStrategy(int pdg_particle) {
  auto it = smear_strategies.find(pdg_particle);
  if (it != smear_strategies.end()) {
    return it->second.get();
  }
  return nullptr;
}

void initializeGaussianSmearStrategy() {
  smear_strategies[11] = std::make_unique<SmearGaussianDirection>(2.91);
  smear_strategies[13] = std::make_unique<SmearGaussianDirection>(1.75);
  smear_strategies[22] = std::make_unique<SmearGaussianDirection>(1.75);
}
