#pragma once

#include "ROOT/RDF/InterfaceUtils.hxx"
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

ROOT::RDF::RNode TrackerPrepare(ROOT::RDF::RNode df);
ROOT::RDF::RNode TrackerPrepareNeutrino(ROOT::RDF::RNode df);
ROOT::RDF::RNode TrackerPrepareGENIE(ROOT::RDF::RNode df);
