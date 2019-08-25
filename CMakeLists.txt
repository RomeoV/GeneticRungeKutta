cmake_minimum_required(VERSION 3.10)
project(GeneticRungeKutta)

# set(CMAKE_CXX_FLAGS="-ggdb")
include(ExternalProject)
externalproject_add(fmt-lib
  PREFIX external
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
)
add_subdirectory(build/external/src/fmt)
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
find_package(TBB REQUIRED)

add_library(GeneticScheme
  src/genetic.cpp
  src/vector_operations.cpp
)
target_link_libraries(GeneticScheme
  TBB::tbb
)
target_include_directories(GeneticScheme
  PUBLIC
  src
)
set_target_properties(GeneticScheme
  PROPERTIES
  CXX_STANDARD 17
)

add_executable(run
  src/main.cpp
)
target_link_libraries(run
  GeneticScheme
)
set_target_properties(run
  PROPERTIES
  CXX_STANDARD 17
)

enable_testing()
add_subdirectory(test)
