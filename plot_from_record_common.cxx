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

int main(int argc, char **argv) {
  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  auto [input_files, output_path, _] = parse_command_line(argc, argv);
  if (!output_path.ends_with(".root")) {
    output_path += ".root";
  }
  std::string basename = output_path;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms{};
  auto df_all = FilterTrackedRDF{ROOT::RDataFrame{"outtree", input_files}}
                    .SetWeightColumnName("weight");
  histograms.emplace_back(
      make_plot(df_all, {"nrings", "nrings", 20, -0.5, 19.5}, "nrings"));
  ROOT::RDF::TH1DModel inv_mass_model{"inv_mass_epip_system", "inv mass", 400,
                                      0.0, 1.2};
  ROOT::RDF::TH1DModel momentum_model{"inv_mass_epip_system", "inv mass", 400,
                                      0.0, 0.8};
  ROOT::RDF::TH1DModel angle_model{"inv_mass_epip_system", "inv mass", 400, 0.0,
                                   M_PI};
  ROOT::RDF::TH1DModel angle_diff_model{"inv_mass_epip_system", "inv mass", 400,
                                        0.0, 0.5};
  ROOT::RDF::TH1DModel angle_diff_model_1{"inv_mass_epip_system", "inv mass",
                                          400, -0.5, 0.5};
  ROOT::RDF::TH1DModel momentum_ratio_model{"momentum_ratio_epip_system",
                                            "momentum ratio", 400, 0.0, 2.0};
  std::vector<std::string> to_snapshot{
      "raw_proton_momentum", "raw_mass_proton", "raw_pi0_before_fsi",
      "channel_name",
      // cut and topology info
      "nrings", "nshower_rings", "nmichel_electrons", "nrings_cut",
      "shower_ring_cut", "nmichel_electrons_cut",
      // dummy weight
      "weight", "kf_chi2", "truth_chi2"};
  auto topo_cut_pass =
      df_all
          .FilterTracked([](bool x) { return x; }, {"nrings_cut"},
                         "2 or 3 rings topology")
          .FilterTracked([](size_t x) { return x == 3; }, {"nrings"},
                         "3 rings topology")
          .FilterTracked([](bool x) { return x; }, {"shower_ring_cut"},
                         "0 (epi) or 1(mupi) shower-like ring")
          .FilterTracked([](bool x) { return x; }, {"nmichel_electrons_cut"},
                         "0 (epi) or 1(mupi) michel electron")
          .Redefine("kf",
                    [](const EventRec &smeared, const EventRec &kf) {
                      if (kf.is_valid) {
                        return kf;
                      }
                      return smeared;
                    },
                    {"smeared", "kf"});

  auto confs = std::to_array({"truth", "smeared", "kf"});
  auto systems = std::to_array({"pi0_system", "epi_system", "primary_lepton",
                                "sublead_photon", "lead_photon"});
  ROOT::RDF::RNode df_with_define = topo_cut_pass;
  for (auto &&conf : confs) {
    df_with_define =
        df_with_define
            .Define(std::format("{}_pi0_system", conf),
                    [](const EventRec &rec) { return rec.gamma1 + rec.gamma2; },
                    {conf})
            .Define(std::format("{}_primary_lepton", conf),
                    [](const EventRec &rec) { return rec.lepton; }, {conf})
            .Define(std::format("{}_epi_system", conf),
                    [](const EventRec &rec) {
                      return rec.lepton + rec.gamma1 + rec.gamma2;
                    },
                    {conf})
            .Define(std::format("{}_lead_photon", conf),
                    [](const EventRec &rec) { return rec.gamma1; }, {conf})
            .Define(std::format("{}_sublead_photon", conf),
                    [](const EventRec &rec) { return rec.gamma2; }, {conf});

    for (auto &&[id, system1] : systems | std::views::enumerate) {
      for (auto &&system2 : systems | std::views::drop(id + 1)) {
        // angle between different particles in the same reconstruction
        auto system1_full_name = std::format("{}_{}", conf, system1);
        auto system2_full_name = std::format("{}_{}", conf, system2);
        auto angle_var_name =
            std::format("angle_{}_{}_{}", conf, system1, system2);
        df_with_define =
            df_with_define.Define(angle_var_name, angle_between,
                                  {system1_full_name, system2_full_name});
        histograms.emplace_back(
            make_plot(df_with_define, angle_model, angle_var_name));
        to_snapshot.emplace_back(angle_var_name);
      }
    }
  }

  for (auto &&[named_define, system_name] :
       std::views::cartesian_product(confs, systems)) {
    auto system_full_name = std::format("{}_{}", named_define, system_name);
    auto mass_var_name = system_full_name + "_m";
    auto momentum_var_name = system_full_name + "_p";
    auto energy_var_name = system_full_name + "_E";
    df_with_define =
        df_with_define
            .Define(mass_var_name, [](const momentum_t &p4) { return p4.M(); },
                    {system_full_name})
            .Define(momentum_var_name,
                    [](const momentum_t &p4) { return p4.P(); },
                    {system_full_name})
            .Define(energy_var_name,
                    [](const momentum_t &p4) { return p4.E(); },
                    {system_full_name});
    to_snapshot.emplace_back(energy_var_name);
    to_snapshot.emplace_back(mass_var_name);
    to_snapshot.emplace_back(momentum_var_name);
    histograms.emplace_back(
        make_plot(df_with_define, inv_mass_model, mass_var_name));
    histograms.emplace_back(
        make_plot(df_with_define, momentum_model, momentum_var_name));
    histograms.emplace_back(
        make_plot(df_with_define, momentum_model, energy_var_name));
  }

  for (auto &&[id1, conf1] : confs | std::views::enumerate) {
    for (auto &&conf2 : confs | std::views::drop(id1 + 1)) {
      // angle and momentum ratio between same particles in different
      // reconstructions
      for (auto &&system_name : systems) {
        auto system1_full_name = std::format("{}_{}", conf1, system_name);
        auto system2_full_name = std::format("{}_{}", conf2, system_name);
        auto angle_var_name =
            std::format("angle_{}_{}_{}", conf1, conf2, system_name);
        auto momentum_ratio_var_name =
            std::format("momentum_ratio_{}_{}_{}", conf1, conf2, system_name);
        auto energy_ratio_var_name =
            std::format("energy_ratio_{}_{}_{}", conf1, conf2, system_name);
        df_with_define =
            df_with_define
                .Define(angle_var_name, angle_between,
                        {system1_full_name, system2_full_name})
                .Define(momentum_ratio_var_name,
                        [](const momentum_t &p4_1, const momentum_t &p4_2) {
                          return p4_2.P() / p4_1.P();
                        },
                        {system1_full_name, system2_full_name})
                .Define(energy_ratio_var_name,
                        [](const momentum_t &p4_1, const momentum_t &p4_2) {
                          return p4_2.E() / p4_1.E();
                        },
                        {system1_full_name, system2_full_name});
        histograms.emplace_back(
            make_plot(df_with_define, angle_diff_model, angle_var_name));
        histograms.emplace_back(make_plot(df_with_define, momentum_ratio_model,
                                          momentum_ratio_var_name));
        histograms.emplace_back(make_plot(df_with_define, momentum_ratio_model,
                                          energy_ratio_var_name));
        to_snapshot.emplace_back(energy_ratio_var_name);
        to_snapshot.emplace_back(momentum_ratio_var_name);
        to_snapshot.emplace_back(angle_var_name);
      }
    }
  }

  try {
    df_with_define = df_with_define.Define(
        "alm5d_protonp",
        [](const proton_dof &p) { return p.calc_momentum_proton(); },
        {

            "kf_5d"});
    histograms.emplace_back(
        make_plot(df_with_define, momentum_model, "alm5d_protonp", ""));
    to_snapshot.emplace_back("alm5d_protonp");

    df_with_define = df_with_define.Define(
        "diff_angle",
        [](double angle1, double angle2) { return angle1 - angle2; },
        {"angle_smeared_pi0_system_primary_lepton",
         "angle_truth_pi0_system_primary_lepton"});
    histograms.emplace_back(
        make_plot(df_with_define, angle_diff_model_1, "diff_angle", ""));
    to_snapshot.emplace_back("diff_angle");
  } catch (const std::exception &e) {
    std::println("Exception caught: {}", e.what());
  }

  // df_with_define =
  //     df_with_define
  //         .Define("dme1", [](const gamma_dof &kf_3d) { return kf_3d.E1; },
  //                 {"kf_3d"})
  //         // .Define("dme2", "kf_3d.E2")
  //         // .Define("alme1", "kf_3d_1.E1")
  //         // .Define("alme2", "kf_3d_1.E2");
  //         .Define("dme2", [](const gamma_dof &kf_3d) { return kf_3d.E2; },
  //                 {"kf_3d"})
  //         .Define("alme1", [](const gamma_dof &kf_3d_1) { return kf_3d_1.E1;
  //         },
  //                 {"kf_3d_1"})
  //         .Define("alme2", [](const gamma_dof &kf_3d_1) { return kf_3d_1.E2;
  //         },
  //                 {"kf_3d_1"});

  // histograms.emplace_back(make_plot(
  //     df_with_define
  //         .Define("e_pi0_3dof",
  //                 [](const double e1, const double e2) { return e1 + e2; },
  //                 {"E1", "E2"})
  //         .Filter([](double e_pi0_3dof) { return e_pi0_3dof > 0.0; },
  //                 {"e_pi0_3dof"}, "e_pi0_3dof positive"),
  //     momentum_model, "e_pi0_3dof", ""));

  // histograms.emplace_back(
  //     make_plot(df_with_define, angle_model, "theta_oa", ""));
  // histograms.emplace_back(make_plot(
  //     df_with_define
  //         .Filter([](double e1) { return e1 > 0.0; }, {"E1"}, "E1 positive")
  //         .Define("p_pi0_3dof",
  //                 [](const double e1, const double e2, const double oa) {
  //                   return std::sqrt(e1 * e1 + e2 * e2 +
  //                                    2 * e1 * e2 * std::cos(oa));
  //                 },
  //                 {"E1", "E2", "theta_oa"}),
  //     momentum_model, "p_pi0_3dof", ""));

  df_with_define.Snapshot("outtree", (basename + "_with_vars.root"),
                          to_snapshot);

  TFile fout{output_path.c_str(), "RECREATE"};
  for (auto &hist : histograms) {
    hist->SetDirectory(&fout);
    hist->Write();
    hist->SetDirectory(nullptr);
  }
  fout.Close();
  std::println("Wrote output to {}", output_path);
  return 0;
}
