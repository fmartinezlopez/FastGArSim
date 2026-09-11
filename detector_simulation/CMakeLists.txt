#----------------------------------------------------------------------------
# GArSimulation -- Geant4 simulation of the gaseous argon detector
#
# Built as part of the top-level FastGArSim build, or on its own:
#
#   mkdir build && cd build && cmake .. && make -j4
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16...3.21)
project(GArSimulation CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include(FastGArSimCommon)

fastgarsim_standard_settings()

#----------------------------------------------------------------------------
# Find Geant4, activating all available UI and Vis drivers by default. Set
# WITH_GEANT4_UIVIS to OFF to build a batch-mode-only executable.
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Geom and Gdml are needed by the GeoVis tool; the interpreter autoloads
# them, a compiled tool has to link them.
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Geom Gdml)
include(${ROOT_USE_FILE})
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Shared data types and their dictionaries
fastgarsim_use_common()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${FASTGARSIM_COMMON_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Sources and headers
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# The executable. The simulation data types come from the dictionary library
# built in common/, so the classes are defined in exactly one place.
add_executable(GArSimulation GArSimulation.cc ${sources} ${headers})
target_link_libraries(GArSimulation
                      ${Geant4_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      FastGArSimCommon)
set_target_properties(GArSimulation PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Copy the run macros next to the executable, which expects to find them in the
# working directory
fastgarsim_copy_files(DESTINATION macros
                      FILES macros/init.mac
                            macros/init_base.mac
                            macros/init_gun.mac
                            macros/run.mac
                            macros/nu.mac
                            macros/gun.mac
                            macros/vis.mac)

#----------------------------------------------------------------------------
# Copy the ROOT utility macros to the build directory.
fastgarsim_copy_files(DESTINATION .
                      FILES utils/GeoVis.C)

#----------------------------------------------------------------------------
# One command line tool per macro in utils/, named after the macro. Adding a
# macro there is all it takes to get an executable for it.
fastgarsim_add_macro_apps(DIRECTORY utils
                          LINK FastGArSimCommon)

#----------------------------------------------------------------------------
fastgarsim_register_component(NAME simulation
                              BIN_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                              MACRO_DIRS ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS GArSimulation DESTINATION bin)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/macros/ DESTINATION share/GArSimulation/macros)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/utils/  DESTINATION share/GArSimulation/utils)

#----------------------------------------------------------------------------
# Standalone build: write the environment files here
if(NOT FASTGARSIM_SUPERBUILD)
  fastgarsim_write_env()
  fastgarsim_print_summary()
endif()
