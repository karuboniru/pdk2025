cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Find Boost components
find_package(Boost REQUIRED COMPONENTS program_options)

# Enable LINK_WHAT_YOU_USE
set(CMAKE_LINK_WHAT_YOU_USE TRUE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED)
find_package(EvtGen REQUIRED)
include_directories(
  ${EVTGEN_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
)

add_compile_definitions(DATA_PATH="${CMAKE_SOURCE_DIR}/data")

file(
  GLOB
  SOURCES_EVENT
  CONFIGURE_DEPENDS
  src/*.cxx
)
add_library(event OBJECT
  ${SOURCES_EVENT}
)

target_link_libraries(event
  PUBLIC
  ROOT::Physics
  ROOT::EG
  boost_program_options
  EvtGen::EvtGen
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
