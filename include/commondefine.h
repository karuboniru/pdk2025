#pragma once

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDataFrame.hxx>

std::tuple<ROOT::RDF::RNode, std::vector<std::string>, std::vector<std::string>,
           std::vector<std::string>>
DefineForEPi(ROOT::RDF::RNode all_with_vars);

std::array<ROOT::RDF::RNode,3> FilterSignalKinematics(ROOT::RDF::RNode df);
