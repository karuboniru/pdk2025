#pragma once
#include <string>
#include <vector>

struct configuration {
  std::vector<std::string> input_files;
  std::string output_file;
  bool genie_mode = false;
};

extern bool nofsi;
extern bool is_mupi;
extern bool do_n_tagging;
extern double n_tagging_eff;

configuration parse_command_line(int argc, char **argv);
