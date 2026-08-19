#----------------------------------------------------------------------------
# GArReconstruction -- modular reconstruction framework
#
# Built as part of the top-level FastGArSim build, or on its own:
#
#   mkdir build && cd build && cmake .. && make -j4
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16...3.21)
project(GArReconstruction CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include(FastGArSimCommon)

fastgarsim_standard_settings()

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore Hist Graf Gpad)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Shared data types and their dictionaries: the simulation types this reads and
# the digitization types it writes both come from common/
fastgarsim_use_common()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${FASTGARSIM_COMMON_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Sources and headers; new modules are picked up automatically
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Dictionary for the reconstruction output types (RecoTrack, RecoCluster, ...),
# which are specific to this package
fastgarsim_add_dictionary(RecoDataDict
                          HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/RecoDataTypes.hh
                          LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/RecoDataTypesLinkDef.h)

#----------------------------------------------------------------------------
# Main executable
add_executable(GArReconstruction GArReconstruction.cc ${sources} ${headers})
target_link_libraries(GArReconstruction
                      ${ROOT_LIBRARIES}
                      RecoDataDictLib
                      FastGArSimCommon)
set_target_properties(GArReconstruction PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Copy the configuration macros next to the executable
fastgarsim_copy_files(DESTINATION macros
                      FILES macros/example.mac
                            macros/ecal_digi.mac
                            macros/tpc_reco.mac)

#----------------------------------------------------------------------------
# Copy the ROOT utility macros to the build directory
fastgarsim_copy_files(DESTINATION .
                      FILES utils/RecoExample.C
                            utils/ECalDigiAnalysis.C
                            utils/TPCRecoAnalysis.C)

#----------------------------------------------------------------------------
# One command line tool per macro in utils/, named after the macro
fastgarsim_add_macro_apps(DIRECTORY utils
                          LINK FastGArSimCommon RecoDataDictLib)

#----------------------------------------------------------------------------
fastgarsim_register_component(NAME reconstruction
                              BIN_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              LIB_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                              LOAD_LIBS RecoDataDict
                              MACRO_DIRS ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS GArReconstruction DESTINATION bin)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/macros/ DESTINATION share/GArReconstruction/macros)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/utils/  DESTINATION share/GArReconstruction/utils)

#----------------------------------------------------------------------------
# Standalone build: write the environment files here
if(NOT FASTGARSIM_SUPERBUILD)
  fastgarsim_write_env()
  fastgarsim_print_summary()
endif()
