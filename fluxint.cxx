#include <TF1.h>
#include <TSpline.h>
#include <fstream>
#include <generator>
#include <print>
#include <sstream>
#include <string>

std::generator<std::pair<double, double>>
parse_flux_file(std::string filepath) {
  std::ifstream infile(filepath);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open file: " + filepath);
  }
  for (std::string line; std::getline(infile, line);) {
    double energy{}, value{};
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    if (!(iss >> energy >> value)) {
      throw std::runtime_error("Error parsing line: " + line);
    }
    co_yield std::make_pair(energy, value);
  }
}

int main(int argc, char **argv) {
  auto flux_file = argv[1];
  std::vector<double> energies, values;
  for (const auto &[energy, value] : parse_flux_file(flux_file)) {
    energies.push_back(energy);
    values.push_back(value);
  }
  auto emin = energies.front();
  auto emax = 20.;
  TSpline3 flux_spline("flux_spline", energies.data(), values.data(),
                       energies.size());
  TF1 flux_func(
      "flux_func",
      [&flux_spline, emin, emax](double *x, double *) {
        double energy = x[0];
        if (energy < emin || energy > emax) {
          return 0.0;
        }
        return flux_spline.Eval(energy);
      },
      emin, emax, 0);
  auto fluxint = flux_func.Integral(emin, emax);
  std::println("Flux integral from {} to {} is {}", emin, emax, fluxint);
}
