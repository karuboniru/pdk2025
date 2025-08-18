#include "cmdline.h"

#include <boost/program_options.hpp>
#include <iostream>

configuration parse_command_line(int argc, char **argv) {
  namespace po = boost::program_options;
  configuration config;

  po::options_description desc("Allowed options");
  po::positional_options_description pos_desc;
  // clang-format off
  desc.add_options()
    ("help,h", "Display this help message")
    ("input,i", po::value<std::vector<std::string>>(&config.input_files)->required(), "Input ROOT files")
    ("output,o", po::value<std::string>(&config.output_file)->required(), "Output ROOT file");
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
