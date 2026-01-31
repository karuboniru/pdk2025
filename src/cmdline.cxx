#include "cmdline.h"

#include <boost/program_options.hpp>
#include <iostream>

bool nofsi = false;
bool is_mupi = false;
bool do_n_tagging = false;
double n_tagging_eff = 0.0;
double unit = 1.0;
bool external_capture_count = false;

configuration parse_command_line(int argc, char **argv) {
  namespace po = boost::program_options;
  configuration config;

  po::options_description desc("Allowed options");
  po::positional_options_description pos_desc;
  // clang-format off
  desc.add_options()
    ("help,h", "Display this help message")
    ("input,i", po::value<std::vector<std::string>>(&config.input_files)->required(), "Input ROOT files")
    ("output,o", po::value<std::string>(&config.output_file)->required(), "Output ROOT file")
    ("genie-mode,g", po::bool_switch(&config.genie_mode)->default_value(false), "Enable GENIE mode")
    ("no-fsi,n", po::bool_switch(&nofsi)->default_value(false), "Disable FSI simulation")
    ("mu-pi-mode,m", po::bool_switch(&is_mupi)->default_value(false), "Enable mu-pi mode")
    ("n-tagging,t", po::bool_switch(&do_n_tagging)->default_value(false), "Enable neutron tagging")
    ("n-tagging-eff,e", po::value<double>(&n_tagging_eff)->default_value(0.25), "Neutron tagging efficiency (0.0 - 1.0)")
    ("unit,u", po::value<double>(&unit)->default_value(1.0), "Unit conversion factor")
    ("external-capture-count,c", po::bool_switch(&external_capture_count)->default_value(false), "Use external neutron capture count")
    ;
  // clang-format on
  pos_desc.add("input", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .positional(pos_desc)
                  .run(),
              vm);
    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "Error: " << e.what() << "\n";
    std::cerr << desc << "\n";
    exit(EXIT_FAILURE);
  }

  if (vm.contains("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }

  return config;
}
