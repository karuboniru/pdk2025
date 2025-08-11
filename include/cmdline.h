#pragma once
#include <string>
#include <vector>

struct configuration {
  std::vector<std::string> input_files;
  std::string output_file;
};

configuration parse_command_line(int argc, char **argv);
