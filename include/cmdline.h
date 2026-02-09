#pragma once
#include <string>
#include <vector>

struct configuration {
  std::vector<std::string> input_files;
  std::vector<std::string> input_corr;
  std::string output_file;
  bool genie_mode = false;
};

extern bool nofsi;
extern bool is_mupi;
extern bool do_n_tagging;
extern double n_tagging_eff;
extern double unit;
extern bool external_capture_count;

configuration parse_command_line(int argc, char **argv);

#include <format>

// Declaration to be placed in header file:
template <> struct std::formatter<configuration> {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

  auto format(const configuration &config, std::format_context &ctx) const {
    return std::format_to(ctx.out(),
                          "configuration {{\n"
                          "  input_files: [{}],\n"
                          "  input_corr: [{}],\n"
                          "  output_file: \"{}\",\n"
                          "  genie_mode: {}\n"
                          "}}",
                          fmt_vector(config.input_files),
                          fmt_vector(config.input_corr), config.output_file,
                          config.genie_mode);
  }

private:
  static std::string fmt_vector(const std::vector<std::string> &vec) {
    if (vec.empty())
      return "";
    std::string result;
    for (size_t i = 0; i < vec.size(); ++i) {
      if (i > 0)
        result += ", ";
      result += std::format("\"{}\"", vec[i]);
    }
    return result;
  }
};
