#pragma once

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <string>
#include <vector>

class FilterTrackedRDF : public ROOT::RDF::RNode {
private:
  std::vector<std::string> m_tracked_filters;
  std::vector<ROOT::RDF::RResultPtr<double>> m_tracked_weight_sum;
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> m_tracked_count;
  ROOT::RDF::RResultPtr<double> m_initial_weight_sum;
  std::string m_weight_column_name{"weight"};

public:
  using ROOT::RDF::RNode::RNode;
  FilterTrackedRDF() = delete;
  FilterTrackedRDF(const FilterTrackedRDF &) = default;
  FilterTrackedRDF(FilterTrackedRDF &&) = default;
  FilterTrackedRDF &operator=(const FilterTrackedRDF &) = default;
  FilterTrackedRDF &operator=(FilterTrackedRDF &&) = default;
  FilterTrackedRDF &operator=(const ROOT::RDF::RNode &from) {
    ROOT::RDF::RNode::operator=(from);
    return *this;
  }
  FilterTrackedRDF &operator=(ROOT::RDF::RNode &&from) {
    ROOT::RDF::RNode::operator=(std::move(from));
    return *this;
  }

  ~FilterTrackedRDF() = default;

  template <typename T>
  FilterTrackedRDF(T from) : ROOT::RDF::RNode(std::forward<T>(from)) {}

  FilterTrackedRDF SetWeightColumnName(const std::string &name) {
    auto new_node = *this;
    new_node.m_weight_column_name = name;
    new_node.m_initial_weight_sum = this->Sum(name);
    if (new_node.m_tracked_filters.size() != 0) {
      throw std::runtime_error(
          "SetWeightColumnName must be called before any FilterTracked");
    }
    m_tracked_weight_sum.push_back(new_node.m_initial_weight_sum);
    return new_node;
  }

  template <typename F>
  FilterTrackedRDF FilterTracked(F &&func, const std::vector<std::string> &cols,
                                 const std::string &name) {
    FilterTrackedRDF new_node = *this;
    new_node =
        ROOT::RDF::RNode{this->Filter(std::forward<F>(func), cols, name)};
    new_node.m_tracked_filters.push_back(name);
    new_node.m_tracked_weight_sum.push_back(
        new_node.Sum(new_node.m_weight_column_name));
    new_node.m_tracked_count.push_back(new_node.Count());
    return new_node;
  }

  template <typename F>
  FilterTrackedRDF Define(const std::string &name, F &&func,
                          const std::vector<std::string> &cols) {
    FilterTrackedRDF new_node = *this;
    new_node = ROOT::RDF::RNode{static_cast<ROOT::RDF::RNode *>(this)->Define(
        name, std::forward<F>(func), cols)};
    return new_node;
  }

  void Report();
};

std::tuple<ROOT::RDF::RNode, std::vector<std::string>, std::vector<std::string>,
           std::vector<std::string>>
DefineForEPi(ROOT::RDF::RNode all_with_vars);

std::tuple<FilterTrackedRDF, std::vector<std::string>, std::vector<std::string>,
           std::vector<std::string>>
DefineForEPi(FilterTrackedRDF all_with_vars);

std::array<ROOT::RDF::RNode, 3> FilterSignalKinematics(ROOT::RDF::RNode df);
std::array<FilterTrackedRDF, 3> FilterSignalKinematics(FilterTrackedRDF df);
