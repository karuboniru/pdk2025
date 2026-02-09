#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RVec.hxx>
#include <TDatabasePDG.h>
#include <TMemFile.h>
#include <TROOT.h>
#include <array>
#include <format>
#include <print>
#include <ranges>
#include <regex>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPie.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "EvtTracker2event.h"
#include "ROOT/RResultPtr.hxx"
#include "cmdline.h"
#include "commondefine.h"
#include "data.h"
#include "decaygen.h"
#include "event.h"
#include "kf.h"
#include "smear.h"

ROOT::RDF::RResultPtr<TH1D> make_plot(auto df, ROOT::RDF::TH1DModel model,
                                      const std::string &varname,
                                      const std::string &prefix = "") {
  model.fName = prefix + varname;
  model.fTitle = prefix + varname;
  return df.Histo1D(model, varname);
}

double angle_between(const momentum_t &p4_1, const momentum_t &p4_2) {
  auto dot_product = p4_1.Vect().Dot(p4_2.Vect());
  auto magnitude_product = p4_1.P() * p4_2.P();
  auto cos_angle = dot_product / magnitude_product;
  // Clamp the value to the valid range for acos
  cos_angle = std::clamp(cos_angle, -1.0, 1.0);
  return std::acos(cos_angle);
}

template <typename CharT>
struct std::formatter<ROOT::Math::PxPyPzEVector, CharT> {
  // Reuse the formatter for a double
  std::formatter<double, CharT> double_fmt;

  // Parse format specification and forward it to the double formatter
  constexpr auto parse(std::basic_format_parse_context<CharT> &ctx) {
    return double_fmt.parse(ctx);
  }

  template <typename FormatContext>
  auto format(const ROOT::Math::PxPyPzEVector &p, FormatContext &ctx) const {
    auto out = ctx.out();

    // Write '('
    *out++ = CharT('(');

    // Format each component with the stored double formatter
    out = double_fmt.format(p.x(), ctx);
    *out++ = CharT(',');
    out = double_fmt.format(p.y(), ctx);
    *out++ = CharT(',');
    out = double_fmt.format(p.z(), ctx);
    *out++ = CharT(',');
    out = double_fmt.format(p.t(), ctx);

    // Write ')'
    *out++ = CharT(')');

    return out;
  }
};

int main(int argc, char **argv) {
  TH1::AddDirectory(false);
  auto [input_files, input_corr, output_path, _] = parse_command_line(argc, argv);
  if (!output_path.ends_with(".root")) {
    output_path += ".root";
  }
  std::string basename = output_path;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms;
  ROOT::RDataFrame{"outtree", input_files}
      .Filter([](bool x) { return x; }, {"nrings_cut"}, "2 or 3 rings topology")
      .Filter([](size_t x) { return x == 3; }, {"nrings"}, "3 rings topology")
      .Filter([](bool x) { return x; }, {"shower_ring_cut"},
              "0 (epi) or 1(mupi) shower-like ring")
      .Filter([](bool x) { return x; }, {"nmichel_electrons_cut"},
              "0 (epi) or 1(mupi) michel electron")
      .Filter([](const EventRec &kf) { return kf.has_gamma2 && kf.is_valid; },
              {"kf"}, "low chi2")
      .Range(0, 20)
      .Foreach(
          [](const EventRec &truth, const EventRec &smeared, const EventRec &kf,
             double kf_chi2, double truth_chi2) {
            std::println(
                "Truth:    gamma1 {:.4f}, gamma2 {:.4f}, lepton {:.4f}",
                truth.gamma1, truth.gamma2, truth.lepton);
            std::println(
                "Smeared:  gamma1 {:.4f}, gamma2 {:.4f}, lepton {:.4f}",
                smeared.gamma1, smeared.gamma2, smeared.lepton);
            std::println(
                "KF:       gamma1 {:.4f}, gamma2 {:.4f}, lepton {:.4f}",
                kf.gamma1, kf.gamma2, kf.lepton);
            std::println("kf_chi2: {:.4f}, truth_chi2 {:.4f}\n", kf_chi2, truth_chi2);
          },
          {"truth", "smeared", "kf", "kf_chi2", "truth_chi2"});

  return 0;
}
