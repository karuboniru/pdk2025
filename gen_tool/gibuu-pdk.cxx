#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TVector3.h>
#include <array>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <ranges>
#include <sstream>
#include <string>

double get_momentum_in_rest_frame(double in_mass, double out1_mass,
                                  double out2_mass) {
  const auto var = pow(in_mass, 4) + pow(out1_mass, 4) + pow(out2_mass, 4) -
                   2 * pow(in_mass, 2) * pow(out1_mass, 2) -
                   2 * pow(in_mass, 2) * pow(out2_mass, 2) -
                   2 * pow(out1_mass, 2) * pow(out2_mass, 2);
  if (var > 0)
    return sqrt(var) / (2 * in_mass);
  else
    return -1;
}

int main(int argc, char const *argv[]) {
  using json = nlohmann::json;
  json j;
  {
    if (argc != 2) {
      std::cout << "Usage: " << argv[0] << " <json file>" << std::endl;
      return 1;
    }
    std::ifstream i(argv[1]);
    i >> j;
  }
  std::ifstream gibuu_in(j["in"].get<std::string>());
  std::ofstream gibuu_out(j["out"].get<std::string>(), std::ios::trunc);
  std::ofstream out_detail(j["out_detail"].get<std::string>(), std::ios::trunc);
  std::string_view spliter{"|"};
  auto particles = j["particles"];
  auto mass1 = TDatabasePDG::Instance()
                   ->GetParticle(particles[0]["pdgid"].get<int>())
                   ->Mass();
  auto mass2 = TDatabasePDG::Instance()
                   ->GetParticle(particles[1]["pdgid"].get<int>())
                   ->Mass();
  size_t id{1};
  for (std::string in; std::getline(gibuu_in, in);) {
    if (in[0] == 'E')
      continue;
    std::vector<std::string> splited{};
    boost::split(splited, in, boost::is_any_of(spliter));
    {
      auto &&prop = splited[1];
      auto pos = prop.find_first_of(':');
      int id, charge;
      std::stringstream ss{prop.substr(pos + 1)};
      ss >> id >> charge;
      // std::cout << id << " " << charge << std::endl;
      if (charge == 0)
        continue;
    }
    auto &&p4 = splited[2];
    auto &&x3 = splited[3];
    std::stringstream ss_p4{p4};
    TLorentzVector p4_v{};
    ss_p4 >> p4_v[3] >> p4_v[0] >> p4_v[1] >> p4_v[2];
    auto mass = p4_v.M();
    if (mass < mass1 + mass2) {
      std::cout << mass << " < " << mass1 << " + " << mass2 << std::endl;
      continue;
    }
    auto boostv = p4_v.BoostVector();
    auto momentum = get_momentum_in_rest_frame(mass, mass1, mass2);
    if (momentum < 0) {
      std::cerr << "Impossible kin found" << std::endl;
      continue;
    }
    else {
      std::cout << "mass    : " << mass << std::endl;
      std::cout << "Momentum: " << momentum << std::endl;
    }
    TVector3 p{};
    gRandom->Sphere(p[0], p[1], p[2], momentum);
    TLorentzVector p4_v_out1{p, sqrt(p.Mag2() + mass1 * mass1)};
    TLorentzVector p4_v_out2{-p, sqrt(p.Mag2() + mass2 * mass2)};
    p4_v_out1.Boost(boostv);
    p4_v_out2.Boost(boostv);
    out_detail << 2212 << " " << p4_v[0] << " " << p4_v[1] << " " << p4_v[2]
               << " " << p4_v[3] << " " //
               << particles[0]["pdgid"].get<int>() << " " << p4_v_out1[0] << " "
               << p4_v_out1[1] << " " << p4_v_out1[2] << " " << p4_v_out1[3]
               << " " //
               << particles[1]["pdgid"].get<int>() << " " << p4_v_out2[0] << " "
               << p4_v_out2[1] << " " << p4_v_out2[2] << " " << p4_v_out2[3]
               << std::endl;
    if (particles[0]["output"].get<bool>())
      gibuu_out << particles[0]["gibuuid"].get<int>() << " "
                << TDatabasePDG::Instance()
                           ->GetParticle(particles[0]["pdgid"].get<int>())
                           ->Charge() /
                       3
                << " " << mass1 << " " << x3 << " " << p4_v_out1.X() << " "
                << p4_v_out1.Y() << " " << p4_v_out1.Z() << " " << id
                << std::endl;
    if (particles[1]["output"].get<bool>())
      gibuu_out << particles[1]["gibuuid"].get<int>() << " "
                << TDatabasePDG::Instance()
                           ->GetParticle(particles[1]["pdgid"].get<int>())
                           ->Charge() /
                       3
                << " " << mass2 << " " << x3 << " " << p4_v_out2.X() << " "
                << p4_v_out2.Y() << " " << p4_v_out2.Z() << " " << id
                << std::endl;
    id++;
  }
  return 0;
}
