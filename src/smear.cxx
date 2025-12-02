#include "smear.h"
#include "local_rand.h"
#include <Math/AxisAngle.h>
#include <Math/Polar3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Rtypes.h>
#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TSpline.h>
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <generator>
#include <map>
#include <memory>

std::generator<std::pair<double, double>>
parse_response_file(std::string filepath) {
  std::ifstream infile(filepath);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file: " + filepath);
  }
  double energy, value;
  for (std::string line; std::getline(infile, line);) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    if (!(iss >> energy >> value)) {
      throw std::runtime_error("Error parsing line: " + line);
    }
    co_yield std::make_pair(energy, value);
  }
  infile.close();
}

class smear_angle_generation {
public:
  static smear_angle_generation &
  get_instance(const std::string &filepath = DATA_PATH "/11/angular_raw") {
    static smear_angle_generation instance(filepath);
    return instance;
  }

  double get_smeared_angle() const {
    double cdf = get_thread_local_random().Uniform(0, 1);
    return InverseCDFSpline.Eval(cdf);
  }

private:
  smear_angle_generation(std::string filepath) {
    constexpr size_t num_points = 500;
    std::vector<double> x_vals, y_vals;
    for (const auto &[x, y] : parse_response_file(filepath)) {
      x_vals.push_back(x);
      y_vals.push_back(y);
    }
    // create a spline from the data
    auto min = x_vals.front();
    auto max = x_vals.back();
    TSpline3 spline("response_spline", x_vals.data(), y_vals.data(),
                    x_vals.size());
    TF1 angle_smear_func(
        "angle_smear_func",
        [&spline](double *x, double *p) { return spline.Eval(x[0]); }, min, max,
        0);
    std::array<double, num_points> x_range{}, cdf_values{};
    double step = (max - min) / (num_points - 1);
    for (size_t i = 0; i < num_points; ++i) {
      x_range[i] = x_vals.front() + (i * step);
      cdf_values[i] = angle_smear_func.Integral(min, x_range[i]);
    }
    // normalize CDF values
    double total_integral = cdf_values[num_points - 1];
    for (auto &val : cdf_values) {
      val /= total_integral;
    }
    InverseCDFSpline = TSpline3("inverse_cdf_spline", cdf_values.data(),
                                x_range.data(), num_points);
  }
  TSpline3 InverseCDFSpline;
};

TSpline3 build_spline_from_file(std::string filepath) {
  std::vector<double> x_vals, y_vals;
  for (const auto &[x, y] : parse_response_file(filepath)) {
    x_vals.push_back(x / 1000.); // convert MeV to GeV
    y_vals.push_back(y);
  }
  return TSpline3("response_spline", x_vals.data(), y_vals.data(),
                  x_vals.size());
}

static std::map<int, std::unique_ptr<ISmearStrategy>> smear_strategies;

// smear: the input PxPyPzEVector to be smeared
// angle: the angle in degrees to smear the vector by
ROOT::Math::PxPyPzEVector smear_angle(const ROOT::Math::PxPyPzEVector &original,
                                      double angle) {
  if (angle == 0.0) {
    return original;
  }
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

  [[nodiscard]] ROOT::Math::PxPyPzEVector
  do_smearing(ROOT::Math::PxPyPzEVector vec) const override {
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
    // double cdf = get_thread_local_random().Uniform(0, 1);
    // double angle = m_sigma * std::sqrt(-2.0 * std::log(cdf));
    auto angle = smear_angle_generation::get_instance().get_smeared_angle() /
                 2.91 * m_sigma;
    // cap angle at 25
    angle = std::min(angle, 25.0);
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
    momentum_scaling = std::max(0.0, momentum_scaling);
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

class SplineBasedSmear final : public ISmearStrategy {
public:
  SplineBasedSmear(const TSpline3 &angle_spline,
                   const TSpline3 &momentum_spline, double scale_ang = 1.,
                   double scale_mom = 1.)
      : m_angle_spline(angle_spline), m_momentum_spline(momentum_spline),
        m_scale_angle_smear(scale_ang), m_scale_momentum_smear(scale_mom) {}
  SplineBasedSmear(const SplineBasedSmear &) = default;
  SplineBasedSmear(SplineBasedSmear &&) = default;
  SplineBasedSmear &operator=(const SplineBasedSmear &) = default;
  SplineBasedSmear &operator=(SplineBasedSmear &&) = default;
  virtual ~SplineBasedSmear() = default;

  [[nodiscard]] ROOT::Math::PxPyPzEVector
  do_smearing(ROOT::Math::PxPyPzEVector vec) const override {
    double energy = vec.E();
    double angle_sigma = m_angle_spline.Eval(energy) * m_scale_angle_smear;
    double momentum_frac =
        m_momentum_spline.Eval(energy) / 100. * m_scale_momentum_smear;

    auto smeared_vec = SmearRayleighDirection{angle_sigma}.do_smearing(vec);
    double momentum_scaling = get_thread_local_random().Gaus(1, momentum_frac);

    momentum_scaling = std::max(0.0, momentum_scaling);
    return smear_momentum(smeared_vec, momentum_scaling);
  }

  [[nodiscard]] double get_sigma_energy(double energy) const override {
    double momentum_frac =
        m_momentum_spline.Eval(energy) / 100. * m_scale_momentum_smear;
    return momentum_frac;
  }

  [[nodiscard]] double get_sigma_angle(double energy) const override {
    double angle_sigma = m_angle_spline.Eval(energy) * m_scale_angle_smear;
    return angle_sigma;
  }

private:
  TSpline3 m_angle_spline;
  TSpline3 m_momentum_spline;
  double m_scale_angle_smear;
  double m_scale_momentum_smear;
};

void initializeGaussianSmearStrategy() {
  auto ang_file = DATA_PATH "/11/angular";
  auto mom_file = DATA_PATH "/11/momentum";
  auto ang_spline = build_spline_from_file(ang_file);
  auto mom_spline = build_spline_from_file(mom_file);
  double scale_ang = 1.8;
  double scale_mom = 0.8;
  smear_strategies[11] = std::make_unique<SplineBasedSmear>(
      ang_spline, mom_spline, scale_ang, scale_mom);
  smear_strategies[-11] = std::make_unique<SplineBasedSmear>(
      ang_spline, mom_spline, scale_ang, scale_mom);
  smear_strategies[13] = std::make_unique<SplineBasedSmear>(
      ang_spline, mom_spline, scale_ang, scale_mom);
  smear_strategies[-13] = std::make_unique<SplineBasedSmear>(
      ang_spline, mom_spline, scale_ang, scale_mom);
  // smear_strategies[22] =
  //     std::make_unique<SplineBasedSmear>(ang_spline, mom_spline, 0.78, 1.25);
  smear_strategies[22] =
      std::make_unique<SplineBasedSmear>(ang_spline, mom_spline, 0.9, 1.3);
}
