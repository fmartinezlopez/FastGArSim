#----------------------------------------------------------------------------
# FastGArSim -- top-level build
#
# Builds the whole chain in one go:
#
#   mkdir build && cd build
#   cmake .. && make -j4
#
# Individual sub-packages can be switched off, either one at a time
#
#   cmake -DBUILD_SIMULATION=OFF ..
#
# or by listing the ones wanted (this overrides the BUILD_* options)
#
#   cmake -DCOMPONENTS="reco,analysis" ..
#
# Each sub-package also still configures on its own, from its own directory.
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16...3.21)
project(FastGArSim VERSION 1.0.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(FastGArSimCommon)

# Tells the sub-packages they are part of this build, and not to set up the
# environment files themselves
set(FASTGARSIM_SUPERBUILD ON)

fastgarsim_standard_settings()

#----------------------------------------------------------------------------
# What to build
option(BUILD_SIMULATION     "Build the Geant4 detector simulation" ON)
option(BUILD_RECONSTRUCTION "Build the reconstruction framework"   ON)
option(BUILD_ANALYSIS       "Build the analysis framework"         ON)
option(WITH_GEANT4_UIVIS    "Build the simulation with Geant4 UI and Vis drivers" ON)

set(COMPONENTS "" CACHE STRING
    "Sub-packages to build, e.g. \"sim,analysis\"; empty means all of them")

if(COMPONENTS)
  string(REPLACE "," ";" _requested "${COMPONENTS}")
  string(REPLACE " " ";" _requested "${_requested}")

  set(BUILD_SIMULATION     OFF)
  set(BUILD_RECONSTRUCTION OFF)
  set(BUILD_ANALYSIS       OFF)

  foreach(_component ${_requested})
    string(TOLOWER "${_component}" _component)
    if(_component MATCHES "^(sim|simulation|detector_simulation)$")
      set(BUILD_SIMULATION ON)
    elseif(_component MATCHES "^(reco|reconstruction)$")
      set(BUILD_RECONSTRUCTION ON)
    elseif(_component MATCHES "^(ana|analysis)$")
      set(BUILD_ANALYSIS ON)
    elseif(_component STREQUAL "all")
      set(BUILD_SIMULATION ON)
      set(BUILD_RECONSTRUCTION ON)
      set(BUILD_ANALYSIS ON)
    else()
      message(FATAL_ERROR
              "Unknown component \"${_component}\". "
              "Valid names: sim, reco, analysis, all.")
    endif()
  endforeach()
endif()

if(NOT BUILD_SIMULATION AND NOT BUILD_RECONSTRUCTION AND NOT BUILD_ANALYSIS)
  message(FATAL_ERROR "Nothing selected to build.")
endif()

#----------------------------------------------------------------------------
# Shared data types and their dictionaries, needed by the simulation and the
# reconstruction and useful to the analysis
add_subdirectory(common)

if(BUILD_SIMULATION)
  add_subdirectory(detector_simulation)
endif()

if(BUILD_RECONSTRUCTION)
  add_subdirectory(reconstruction)
endif()

if(BUILD_ANALYSIS)
  add_subdirectory(analysis)
endif()

#----------------------------------------------------------------------------
# setup.sh and rootlogon.C covering everything that was built
fastgarsim_write_env()
fastgarsim_print_summary()
