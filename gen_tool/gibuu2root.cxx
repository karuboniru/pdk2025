// gibuu2root: Convert GiBUU FinalEvents.dat output to a ROOT TTree in StdHep
// format.
//
// GiBUU writes one text line per final-state particle with the columns:
//   Run  Event  ID  Charge  perweight  pos(1-3)  mom(0-3)  history
//   production_ID  enu
//
// Lines whose perweight column equals 0 describe the struck nucleon (target
// particle) rather than an outgoing particle; those are stored separately as
// the initial-state nucleon.
//
// The output ROOT file contains a single TTree named "out_tree" with branches:
//   StdHepN        -- number of particles in the event (outgoing + 2 initial)
//   StdHepPdg[N]   -- PDG codes; index 0 = incoming neutrino, 1 = target
//                     nucleon, 2..N-1 = outgoing particles
//   StdHepStatus[N]-- 0 for initial-state particles, 1 for final-state
//   StdHepP4[N][4] -- 4-momenta (E, px, py, pz) in GeV
//   StdHepX4[N][4] -- 4-position (t, x, y, z) in fm
//   weight         -- GiBUU perturbative weight (perweight) for the event
//   channel        -- GiBUU production_ID (interaction channel code)
//
// Usage:
//   gibuu2root <input_file> <output_file> [neutrino_pdg_id]
//
// neutrino_pdg_id defaults to 14 (nu_mu). Use 12 for nu_e, -14 for anti-nu_mu,
// -12 for anti-nu_e.
//
// The GiBUU-to-PDG ID mapping is read from DATA_PATH/gibuudata.dat (the same
// mapping file used by pert_to_root and the rest of the analysis pipeline).

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <array>
#include <boost/functional/hash.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

class gibuu_to_pdg {
private:
  std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>
      datamap;
  gibuu_to_pdg(std::string mapfile) {
    std::ifstream ifs(mapfile);
    for (std::string line; std::getline(ifs, line);) {
      std::istringstream iss(line);
      int pdg, gibuu, charge;
      iss >> pdg >> gibuu >> charge;
      datamap[{gibuu, charge}] = pdg;
      if (charge != 0)
        datamap[{gibuu, -charge}] = -pdg;
      datamap[{-gibuu, -charge}] = -pdg;
    }
  }

public:
  static gibuu_to_pdg &get_instance(const std::string &mapfile = "") {
    static gibuu_to_pdg instance{mapfile};
    return instance;
  }
  int get_pdg(int gibuu, int charge) { return datamap[{gibuu, charge}]; }
};

auto get_pdg(int gibuu, int charge) {
  auto res = gibuu_to_pdg::get_instance().get_pdg(gibuu, charge);
  if (res == 0) {
    std::cout << "Error: no PDG code for gibuu=" << gibuu
              << " charge=" << charge << std::endl;
  }
  return res;
}

class event {
public:
  int Run;
  int Event;
  std::vector<int> ID;
  std::vector<int> Charge;
  double perweight;
  std::vector<std::array<double, 3>> position;
  std::vector<std::array<double, 4>> momentum;
  std::array<double, 4> momentum_in{};
  std::array<double, 4> position_in{};
  int ID_in{};
  int Charge_in{};
  std::vector<int> history;
  int production_ID;
  double enu;

  bool initialized{false};
  int particle_count{0};

  std::optional<event> add_line(std::string line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
    if (tokens.size() != 15) {
      std::cout << "Error: line has " << tokens.size() << " tokens, expected 15"
                << std::endl;
      std::cout << line << std::endl;
      return std::nullopt;
    }
    if (initialized) {
      if (Run != std::stoi(tokens[0]) || Event != std::stoi(tokens[1])) {
        event next;
        next.add_line(line);
        return next;
      }
    } else {
      Run = std::stoi(tokens[0]);
      Event = std::stoi(tokens[1]);
      production_ID = std::stoi(tokens[13]);
      perweight = std::stod(tokens[4]);
      enu = std::stod(tokens[14]);
      initialized = true;
    }
    if (std::stod(tokens[4]) == 0) {
      momentum_in[0] = std::stod(tokens[8]);
      momentum_in[1] = std::stod(tokens[9]);
      momentum_in[2] = std::stod(tokens[10]);
      momentum_in[3] = std::stod(tokens[11]);
      position_in[0] = std::stod(tokens[5]);
      position_in[1] = std::stod(tokens[6]);
      position_in[2] = std::stod(tokens[7]);

      ID_in = std::stoi(tokens[2]);
      Charge_in = std::stoi(tokens[3]);
      return std::nullopt;
    }
    ID.emplace_back(std::stoi(tokens[2]));
    Charge.emplace_back(std::stoi(tokens[3]));
    position.emplace_back(std::array<double, 3>{
        std::stod(tokens[5]), std::stod(tokens[6]), std::stod(tokens[7])});
    momentum.emplace_back(
        std::array<double, 4>{std::stod(tokens[8]), std::stod(tokens[9]),
                              std::stod(tokens[10]), std::stod(tokens[11])});
    history.emplace_back(std::stoi(tokens[12]));

    particle_count++;
    return std::nullopt;
  }
};

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0]
              << "<input_file> <output_file> [neutrino id (default 14)]"
              << std::endl;
    return 1;
  }
  // std::string map_file = argv[1];
  std::string map_file = BASEPATH "/data/gibuudata.dat";
  gibuu_to_pdg::get_instance(map_file);
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  int Nuid = argc > 4 ? std::stoi(argv[3]) : 14;
  std::ifstream input(input_file);
  auto out_tree = std::make_unique<TTree>("out_tree", "out_tree");
  constexpr int max_n = 256;
  int StdHepN;
  int StdHepPdg[max_n];
  int StdHepStatus[max_n];
  double StdHepP4[max_n][4];
  double StdHepX4[max_n][4];
  double weight{};
  int channel{};
  // set address for tree
  out_tree->Branch("StdHepN", &StdHepN, "StdHepN/I");
  out_tree->Branch("StdHepPdg", StdHepPdg, "StdHepPdg[StdHepN]/I");
  out_tree->Branch("StdHepStatus", StdHepStatus, "StdHepStatus[StdHepN]/I");
  out_tree->Branch("StdHepP4", StdHepP4, "StdHepP4[StdHepN][4]/D");
  out_tree->Branch("StdHepX4", StdHepX4, "StdHepX4[StdHepN][4]/D");
  out_tree->Branch("weight", &weight, "weight/D");
  out_tree->Branch("channel", &channel, "channel/I");
  // fill tree

  auto fill = [&](event &e) {
    channel = e.production_ID;
    weight = e.perweight;
    StdHepN = e.particle_count + 2;
    // set incoming nu
    StdHepPdg[0] = Nuid;
    StdHepStatus[0] = 0;
    StdHepP4[0][3] = e.enu;
    StdHepP4[0][0] = 0;
    StdHepP4[0][1] = 0;
    StdHepP4[0][2] = e.enu;
    StdHepX4[0][0] = 0;
    StdHepX4[0][1] = 0;
    StdHepX4[0][2] = 0;
    StdHepX4[0][3] = 0;
    // set target nucl
    StdHepPdg[1] = get_pdg(e.ID_in, e.Charge_in);
    StdHepStatus[1] = 0;
    StdHepP4[1][3] = e.momentum_in[0];
    StdHepP4[1][0] = e.momentum_in[1];
    StdHepP4[1][1] = e.momentum_in[2];
    StdHepP4[1][2] = e.momentum_in[3];
    StdHepX4[1][0] = 0;
    StdHepX4[1][1] = e.position_in[0];
    StdHepX4[1][2] = e.position_in[1];
    StdHepX4[1][3] = e.position_in[2];

    // set outgoing particles
    for (int i = 0; i < e.particle_count; i++) {
      StdHepPdg[i + 2] = get_pdg(e.ID[i], e.Charge[i]);
      StdHepStatus[i + 2] = 1;
      StdHepP4[i + 2][3] = e.momentum[i][0];
      StdHepP4[i + 2][0] = e.momentum[i][1];
      StdHepP4[i + 2][1] = e.momentum[i][2];
      StdHepP4[i + 2][2] = e.momentum[i][3];
      StdHepX4[i + 2][0] = 0;
      StdHepX4[i + 2][1] = e.position[i][0];
      StdHepX4[i + 2][2] = e.position[i][1];
      StdHepX4[i + 2][3] = e.position[i][2];
    }
    out_tree->Fill();
  };

  // iterate over the input file line by line
  // per each line we have:
  //  1:Run  2:Event  3:ID 4:Charge
  //  5:perweight   6:position(1)   7:position(2)
  //  8:position(3)   9:momentum(0)  10:momentum(1)
  //  11:momentum(2)  12:momentum(3)
  // 13:history  14:production_ID  15:enu
  for (auto [line, e] = std::tuple<std::string, event>{};
       std::getline(input, line);) {
    if (line[0] == '#')
      continue;
    auto next = e.add_line(line);
    if (next.has_value()) {
      fill(e);
      e = next.value();
    }
  }
  // write tree to file
  TFile f(output_file.c_str(), "RECREATE");
  out_tree->Write();
  f.Close();
  return 0;
}
