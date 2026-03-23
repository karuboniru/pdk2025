#include <TBranch.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <array>
#include <boost/functional/hash.hpp>
#include <cassert>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

class gibuu_to_pdg {
private:
  std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>
      map;
  gibuu_to_pdg() {
    std::ifstream ifs(DATA_PATH "/gibuudata.dat");
    for (std::string line; std::getline(ifs, line);) {
      std::istringstream iss(line);
      int pdg, gibuu, charge;
      iss >> pdg >> gibuu >> charge;
      map[{gibuu, charge}] = pdg;
      if (charge != 0)
        map[{gibuu, -charge}] = -pdg;
    }
  }

public:
  static gibuu_to_pdg &get_instance() {
    static gibuu_to_pdg instance;
    return instance;
  }
  int get_pdg(int gibuu, int charge) { return map[{gibuu, charge}]; }
};

auto get_pdg(int gibuu, int charge) {
  return gibuu_to_pdg::get_instance().get_pdg(gibuu, charge);
}

int main(const int argc, const char **argv) {
  // gInterpreter->GenerateDictionary(
  //     "std::vector<std::pair<int,TLorentzVector> >", "vector");
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
              << " <PertParticles_Init_mom.dat> <PertParticles_Final_mom.dat> "
                 "<output_file>"
              << std::endl;
    return 1;
  }
  const std::filesystem::path input_file_init(argv[1]);
  const std::filesystem::path input_file_final(argv[2]);
  const std::filesystem::path output_file(argv[3]);
  if (!std::filesystem::exists(input_file_init)) {
    std::cerr << "File " << input_file_init << " does not exist" << std::endl;
    return 1;
  }
  if (!std::filesystem::exists(input_file_final)) {
    std::cerr << "File " << input_file_final << " does not exist" << std::endl;
    return 1;
  }
  if (!output_file.parent_path().empty())
    std::filesystem::create_directory(output_file.parent_path());
  auto tree = std::make_unique<TTree>("outtree", "outtree");
  int nparticles{};
  double P[100][4];
  int pdg[100];
  int status[100];
  tree->Branch("nparticles", &nparticles, "nparticles/I");
  tree->Branch("P", P, "P[nparticles][4]/D");
  tree->Branch("pdg", pdg, "pdg[nparticles]/I");
  tree->Branch("status", status, "status[nparticles]/I");

  // int pdgids{};
  TLorentzVector *p4_init{};
  // TLorentzVector *final_partcles;

  std::deque<std::pair<int, TLorentzVector>> particles_init;
  std::deque<std::vector<int>> particle_pdgids;
  std::deque<std::vector<TLorentzVector>> particle_p4;

  std::cout << "Reading " << input_file_init << " and " << input_file_final
            << std::endl;
  std::ifstream input_file_init_stream(input_file_init);
  std::ifstream input_file_final_stream(input_file_final);
  // for (std::string line; std::getline(input_file_init_stream, line);) {

  // }
  auto finish_event = [&]() {
    auto lepton_to_final_state = [&] {
      auto abs_pdg = abs(pdg[nparticles]);
      if (status[nparticles] == 2 && abs_pdg > 10 && abs_pdg < 17) {
        status[nparticles + 1] = 1;
        memcpy(P[nparticles + 1], P[nparticles], sizeof(double) * 4);
        pdg[nparticles + 1] = pdg[nparticles];
        nparticles++;
      }
    };
    std::string line;
    if (!std::getline(input_file_init_stream, line)) {
      return false;
    }
    std::istringstream iss{line};
    status[nparticles] = 0;
    iss >> pdg[nparticles];
    iss >> P[nparticles][1] >> P[nparticles][2] >> P[nparticles][3] >>
        P[nparticles][0];
    nparticles++;
    status[nparticles] = 2;
    iss >> pdg[nparticles];
    iss >> P[nparticles][1] >> P[nparticles][2] >> P[nparticles][3] >>
        P[nparticles][0];
    lepton_to_final_state();
    nparticles++;
    status[nparticles] = 2;
    iss >> pdg[nparticles];
    iss >> P[nparticles][1] >> P[nparticles][2] >> P[nparticles][3] >>
        P[nparticles][0];
    lepton_to_final_state();
    nparticles++;
    tree->Fill();
    nparticles = 0;
    return true;
  };
  int currentid{1};
  int final_event_id{24000};
  for (std::string line; std::getline(input_file_final_stream, line);) {
    if (line.find("#") != std::string::npos) {
      continue;
    }
    std::istringstream iss(line);
    int number, eventid, gibuuid, charge;
    double p0, p1, p2, p3;
    iss >> number >> eventid >> gibuuid >> charge >> p0 >> p1 >> p2 >> p3;
    if (eventid != currentid) {
      // finish_event();
      for (; currentid != eventid; currentid++) {
        finish_event();
      }
    }
    P[nparticles][0] = p0;
    P[nparticles][1] = p1;
    P[nparticles][2] = p2;
    P[nparticles][3] = p3;
    pdg[nparticles] = get_pdg(gibuuid, charge);
    status[nparticles] = 1;
    nparticles++;
  }
  while (finish_event())
    ;

  tree->SaveAs(output_file.c_str());
}
