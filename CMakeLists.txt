cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Find Boost components
find_package(Boost REQUIRED COMPONENTS program_options)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_LINK_DEPENDS_NO_SHARED ON)

find_package(ROOT REQUIRED)
find_package(EvtGen REQUIRED)
find_package(nlohmann_json REQUIRED)
include_directories(
  ${EVTGEN_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
)

add_compile_definitions(DATA_PATH="${CMAKE_SOURCE_DIR}/data")
add_compile_definitions(EVTGEN_DATA_DIR="${EVTGEN_DATA_DIR}")

option(EVTGEN_THREAD_MODE "Enable EvtGen multi-thread mode" ON)
if(EVTGEN_THREAD_MODE)
  add_compile_definitions(EVTGEN_THREAD_MODE)
endif()

file(
  GLOB
  SOURCES_EVENT
  CONFIGURE_DEPENDS
  src/*.cxx
)
add_library(event SHARED
  ${SOURCES_EVENT}
)

target_link_libraries(event
  PUBLIC
  ROOT::Physics
  ROOT::EG
  ROOT::Minuit2
  boost_program_options
  EvtGen::EvtGen
  ROOT::ROOTDataFrame
)

root_generate_dictionary(
  G__event
  include/data.h
  MODULE event
  LINKDEF include/dataLinkDef.h
)

add_executable(main main.cxx)

# Link libraries
target_link_libraries(main
  ROOT::Physics ROOT::ROOTDataFrame
  event
)

target_link_options(main
  PRIVATE
  -Wl,--as-needed
)

add_executable(neutrinoevent neutrinoevent.cxx)

# Link libraries
target_link_libraries(neutrinoevent
  ROOT::Physics ROOT::ROOTDataFrame
  event
)

target_link_options(neutrinoevent
  PRIVATE
  -Wl,--as-needed
)


add_executable(fluxint fluxint.cxx)

# Link libraries
target_link_libraries(fluxint
  ROOT::Physics
  event
)

target_link_options(fluxint
  PRIVATE
  -Wl,--as-needed
)


add_executable(neutrinoeventGenie neutrinoeventGenie.cxx)

# Link libraries
target_link_libraries(neutrinoeventGenie
  ROOT::Physics ROOT::ROOTDataFrame
  event
)

target_link_options(neutrinoeventGenie
  PRIVATE
  -Wl,--as-needed
)


add_executable(genie_rate genie_rate.cxx)

# Link libraries
target_link_libraries(genie_rate
  ROOT::Physics ROOT::ROOTDataFrame
  boost_program_options
)

target_link_options(genie_rate
  PRIVATE
  -Wl,--as-needed
)

add_executable(freeproton freeproton.cxx)

# Link libraries
target_link_libraries(freeproton
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
)

target_link_options(freeproton
  PRIVATE
  -Wl,--as-needed
)

add_executable(record_pdk record_pdk.cxx)

# Link libraries
target_link_libraries(record_pdk
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)


add_executable(plot_from_record_common plot_from_record_common.cxx)

# Link libraries
target_link_libraries(plot_from_record_common
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)


add_executable(dump dump.cxx)

# Link libraries
target_link_libraries(dump
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)

add_executable(construct_th3d construct_th3d.cxx)

target_link_libraries(construct_th3d
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)


add_executable(slice_event slice_event.cxx)

target_link_libraries(slice_event
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)

add_executable(capture_check capture_check.cxx)

target_link_libraries(capture_check
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)

add_executable(eff_check eff_check.cxx)

target_link_libraries(eff_check
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)

add_executable(si_check si_check.cxx)

target_link_libraries(si_check
  ROOT::Physics
  ROOT::EG
  ROOT::ROOTDataFrame
  event
)

# gen_tool: GiBUU-based nucleon decay event generation pipeline
add_executable(gibuu-pdk gen_tool/gibuu-pdk.cxx)
target_link_libraries(gibuu-pdk
  ROOT::Physics
  ROOT::EG
  Boost::headers
  nlohmann_json::nlohmann_json
)

add_executable(pert_to_root gen_tool/pert_to_root.cxx)
target_link_libraries(pert_to_root
  ROOT::Core
  ROOT::Physics
  ROOT::EG
  ROOT::Tree
  Boost::headers
  nlohmann_json::nlohmann_json
)

# gibuu2root: convert GiBUU FinalEvents.dat to a StdHep-layout ROOT TTree.
# Used in the atmospheric neutrino generation pipeline (see example_atmo/).
add_executable(gibuu2root gen_tool/gibuu2root.cxx)
target_compile_definitions(gibuu2root PRIVATE BASEPATH="${CMAKE_SOURCE_DIR}")
target_link_libraries(gibuu2root
  ROOT::Core
  ROOT::Physics
  ROOT::Tree
  Boost::headers
)
