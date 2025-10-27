#include "commondefine.h"
#include "event.h"

auto iter_pair(auto &iter_range) {
  auto en = iter_range | std::views::enumerate;
  return std::views::cartesian_product(en, en) |
         std::views::filter([](const auto &pair) {
           auto [first, second] = pair;
           return std::get<0>(first) < std::get<0>(second);
         }) |
         std::views::transform([](const auto &pair) {
           auto [first, second] = pair;
           return std::make_pair(std::get<1>(first), std::get<1>(second));
         });
}

std::tuple<ROOT::RDF::RNode, std::vector<std::string>, std::vector<std::string>,
           std::vector<std::string>>
DefineForEPi(ROOT::RDF::RNode all_with_vars_in) {
  std::tuple<ROOT::RDF::RNode, std::vector<std::string>,
             std::vector<std::string>, std::vector<std::string>>
      ret{all_with_vars_in, {}, {}, {}};
  auto &[all_with_vars, to_snapshot, mass_list, p_list] = ret;

  constexpr double to_deg = 180. / M_PI;

  auto name_p4_pairs =
      std::to_array({"electron", "lead_photon", "sublead_photon", "pi0_system",
                     "epi_system"});

  for (const auto &name : name_p4_pairs) {
    all_with_vars = all_with_vars
                        .Define(std::format("true_{}_p", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.first.P();
                                },
                                {name})
                        .Define(std::format("true_{}_m", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.first.M();
                                },
                                {name})
                        .Define(std::format("true_{}_theta", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.first.Theta() * to_deg;
                                },
                                {name})
                        .Define(std::format("true_{}_phi", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.first.Phi() * to_deg;
                                },
                                {name})
                        .Define(std::format("smared_{}_p", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.second.P();
                                },
                                {name})
                        .Define(std::format("smared_{}_m", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.second.M();
                                },
                                {name})
                        .Define(std::format("smared_{}_theta", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.second.Theta() * to_deg;
                                },
                                {name})
                        .Define(std::format("smared_{}_phi", name),
                                [](const pair_momentum_t &p4_pair) {
                                  return p4_pair.second.Phi() * to_deg;
                                },
                                {name});
    for (const auto &suffix : std::to_array({"p", "m", "theta", "phi"})) {
      to_snapshot.push_back(std::format("true_{}_{}", name, suffix));
      to_snapshot.push_back(std::format("smared_{}_{}", name, suffix));
    }
    mass_list.push_back(std::format("true_{}_m", name));
    mass_list.push_back(std::format("smared_{}_m", name));
    p_list.push_back(std::format("true_{}_p", name));
    p_list.push_back(std::format("smared_{}_p", name));
  }

  for (auto &&[var1, var2] : iter_pair(name_p4_pairs)) {
    all_with_vars =
        all_with_vars
            .Define(std::format("true_{}_{}_angle", var1, var2),
                    [](const pair_momentum_t &p4_pair1,
                       const pair_momentum_t &p4_pair2) {
                      return std::acos(p4_pair1.first.Vect().Unit().Dot(
                                 p4_pair2.first.Vect().Unit())) *
                             to_deg;
                    },
                    {var1, var2})
            .Define(std::format("smared_{}_{}_angle", var1, var2),
                    [](const pair_momentum_t &p4_pair1,
                       const pair_momentum_t &p4_pair2) {
                      return std::acos(p4_pair1.second.Vect().Unit().Dot(
                                 p4_pair2.second.Vect().Unit())) *
                             to_deg;
                    },
                    {var1, var2});
    to_snapshot.push_back(std::format("true_{}_{}_angle", var1, var2));
    to_snapshot.push_back(std::format("smared_{}_{}_angle", var1, var2));
  }

  return ret;
}

std::array<ROOT::RDF::RNode, 3> FilterSignalKinematics(ROOT::RDF::RNode df) {
  auto signal =
      df.Filter(
            [](double rec_m_pi0) {
              return rec_m_pi0 > 0.085 && rec_m_pi0 < 0.185;
            },
            {"smared_pi0_system_m"}, "reconstructed pi0 mass cut (86-185 MeV)")
          .Filter(
              [](double rec_m_p) { return rec_m_p > 0.8 && rec_m_p < 1.05; },
              {"smared_epi_system_m"},
              "reconstructed proton mass cut (800-1050 MeV)")
          .Filter([](double rec_p_p) { return rec_p_p < 0.25; },
                  {"smared_epi_system_p"},
                  "reconstructed proton momentum cut (< 250 MeV)");
  auto signal_upper = signal.Filter(
      [](double rec_p_p) { return rec_p_p >= 0.1; }, {"smared_epi_system_p"},
      "reconstructed proton momentum upper region (>= 100 MeV)");
  auto signal_lower = signal.Filter(
      [](double rec_p_p) { return rec_p_p < 0.1; }, {"smared_epi_system_p"},
      "reconstructed proton momentum lower region (< 100 MeV)");
  return {signal, signal_upper, signal_lower};
}

std::tuple<FilterTrackedRDF, std::vector<std::string>, std::vector<std::string>,
           std::vector<std::string>>
DefineForEPi(FilterTrackedRDF all_with_vars) {
  auto res_raw = DefineForEPi(static_cast<ROOT::RDF::RNode>(all_with_vars));
  auto &[all_with_vars_node, to_snapshot, mass_list, p_list] = res_raw;
  FilterTrackedRDF all_with_vars_ft{all_with_vars};
  all_with_vars_ft = all_with_vars_node;
  return {all_with_vars_ft, to_snapshot, mass_list, p_list};
}

std::array<FilterTrackedRDF, 3> FilterSignalKinematics(FilterTrackedRDF df) {
  auto signal =
      df.FilterTracked(
            [](double rec_m_pi0) {
              return rec_m_pi0 > 0.085 && rec_m_pi0 < 0.185;
            },
            {"smared_pi0_system_m"}, "reconstructed pi0 mass cut (85-185 MeV)")
          .FilterTracked(
              [](double rec_m_p) { return rec_m_p > 0.8 && rec_m_p < 1.05; },
              {"smared_epi_system_m"},
              "reconstructed proton mass cut (800-1050 MeV)")
          .FilterTracked([](double rec_p_p) { return rec_p_p < 0.25; },
                         {"smared_epi_system_p"},
                         "reconstructed proton momentum cut (< 250 MeV)");
  auto signal_upper = signal.FilterTracked(
      [](double rec_p_p) { return rec_p_p >= 0.1; }, {"smared_epi_system_p"},
      "reconstructed proton momentum upper region (>= 100 MeV)");
  auto signal_lower = signal.FilterTracked(
      [](double rec_p_p) { return rec_p_p < 0.1; }, {"smared_epi_system_p"},
      "reconstructed proton momentum lower region (< 100 MeV)");
  return {signal, signal_upper, signal_lower};
}

#include <print>
#include <ranges>
void FilterTrackedRDF::Report() {
  for (auto &&[weights, name, count] : std::views::zip(
           m_tracked_weight_sum | std::views::adjacent<2>, m_tracked_filters, m_tracked_count)) {
    auto [before, after] = weights;
    double efficiency = after.GetValue() / before.GetValue();
    std::println(
        "Filter {0}\t: efficiency = {1:.2f}% ({2:.2e} -> {3:.2e})\ttotal "
        "efficiency: {4:.4e} ({5:.2e} -> {3:.2e}), count: {6}",
        name, efficiency * 100., before.GetValue(), after.GetValue(),
        (after.GetValue() / m_initial_weight_sum.GetValue()),
        m_initial_weight_sum.GetValue(), count.GetValue());
  }
}
