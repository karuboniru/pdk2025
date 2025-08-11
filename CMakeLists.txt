cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Find Boost components
find_package(Boost REQUIRED COMPONENTS program_options)

# Enable LINK_WHAT_YOU_USE
set(CMAKE_LINK_WHAT_YOU_USE TRUE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED)

file(
  GLOB
  SOURCES_EVENT
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
)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
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
