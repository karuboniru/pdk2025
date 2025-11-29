#include "Math/Vector3Dfwd.h"
#include "Math/Vector4Dfwd.h"
#include "TROOT.h"
#include <Math/DisplacementVector3D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TDatabasePDG.h>
#include <array>
#include <cmath>
#include <random>

double get_momentum_in_rest_frame(double in_mass, double out1_mass,
                                  double out2_mass) {
  using std::pow;
  using std::sqrt;
  const auto var = pow(in_mass, 4) + pow(out1_mass, 4) + pow(out2_mass, 4) -
                   2 * pow(in_mass, 2) * pow(out1_mass, 2) -
                   2 * pow(in_mass, 2) * pow(out2_mass, 2) -
                   2 * pow(out1_mass, 2) * pow(out2_mass, 2);
  if (var > 0)
    return sqrt(var) / (2 * in_mass);
  return -1;
}

int main(int argc, char **agrv) {
  constexpr double mass_proton = 0.9382720813; // GeV/c^2
  constexpr double mass_pion0 = 0.1349768;     // GeV
  const auto pdg_lepton = argc > 1 ? std::atoi(agrv[1]) : -11;
  const auto mass_lepton =
      TDatabasePDG::Instance()->GetParticle(pdg_lepton)->Mass();
  const std::string lepton_name =
      TDatabasePDG::Instance()->GetParticle(pdg_lepton)->GetName();
  const size_t event_count = argc > 2 ? std::atoi(agrv[2]) : 1000000;
  const int cluster_size = event_count / 96;
  const std::string output_filename =
      argc > 3
          ? std::string(agrv[3])
          : std::format("proton_decay_{}pi.{}.root", lepton_name, event_count);

  ROOT::RDF::RSnapshotOptions options;
  options.fAutoFlush = cluster_size;

  ROOT::EnableImplicitMT();

  const auto momentum_electron =
      get_momentum_in_rest_frame(mass_proton, mass_pion0, mass_lepton);

  using p3d = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>;
  using p4d = ROOT::Math::PxPyPzEVector;

  ROOT::RDataFrame{event_count}
      .Define("rand_dir",
              []() {
                thread_local std::mt19937 gen(std::random_device{}());
                auto phi = std::uniform_real_distribution<>(0, 2 * M_PI)(gen);
                auto costheta = std::uniform_real_distribution<>(-1, 1)(gen);
                return ROOT::Math::DisplacementVector3D<
                    ROOT::Math::Cartesian3D<double>>(
                    ROOT::Math::Polar3DVector{1, std::acos(costheta), phi});
              })
      .Define("p4_pion",
              [&](const p3d &p) {
                auto pion_momentum = p * momentum_electron;
                return p4d{ROOT::Math::PxPyPzMVector{
                    pion_momentum.x(), pion_momentum.y(), pion_momentum.z(),
                    mass_pion0}};
              },
              {"rand_dir"})
      .Define("p4_electron",
              [&](const p3d &p) {
                auto electron_momentum = p * (-momentum_electron);
                return p4d{ROOT::Math::PxPyPzMVector{
                    electron_momentum.x(), electron_momentum.y(),
                    electron_momentum.z(), mass_lepton}};
              },
              {"rand_dir"})
      .Define("nparticles", []() { return 5; })
      .Define("P",
              [&](p4d &p_pion, p4d &p_e) {
                p4d p_proton{0, 0, 0, mass_proton};
                std::array<double, 20> ret;
                auto start_index = 0;
                auto fill = [&](const p4d &p) {
                  ret[start_index] = p.E();
                  ret[start_index + 1] = p.Px();
                  ret[start_index + 2] = p.Py();
                  ret[start_index + 3] = p.Pz();
                  start_index += 4;
                };
                fill(p_proton);
                fill(p_pion);
                fill(p_e);
                fill(p_pion);
                fill(p_e);
                return ret;
              },
              {"p4_pion", "p4_electron"})
      .Define(
          "pdg",
          [=]() { return std::array<int, 5>{2212, 111, pdg_lepton, 111, -11}; })
      .Define("status", []() { return std::array<int, 5>{0, 1, 1, 2, 2}; })
      .Snapshot(
          "outtree",
          output_filename,
          {"nparticles", "P", "pdg", "status"}, options);
}
