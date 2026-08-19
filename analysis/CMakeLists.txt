#----------------------------------------------------------------------------
# GArAnalysis -- compiled analysis framework
#
# Built as part of the top-level FastGArSim build, or on its own:
#
#   mkdir build && cd build && cmake .. && make -j4
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16...3.21)
project(GArAnalysis CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include(FastGArSimCommon)

fastgarsim_standard_settings()

# EG provides TDatabasePDG, which the macros use to name particles. The
# interpreter autoloads it; a compiled tool has to link it.
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist Graf Gpad MathCore Physics EG)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Shared data types. The analysis reads flat ntuples and so does not need the
# dictionaries itself, but macros that open a simulation or reconstruction file
# directly do, and they get them from the generated rootlogon.C.
fastgarsim_use_common()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Sources and headers
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Shared analysis library. Macros load this at run time and derive from
# ana::AnalysisBase; no ROOT dictionary is needed because none of these types
# are written to a TTree.
add_library(GArAnalysis SHARED ${sources} ${headers})
target_link_libraries(GArAnalysis PUBLIC ${ROOT_LIBRARIES})
target_include_directories(GArAnalysis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(GArAnalysis PROPERTIES
                      OUTPUT_NAME GArAnalysis
                      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Copy the analysis macros to the build directory, so that ROOT can be run
# straight from there
file(GLOB _macros RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/macros/*.C)
fastgarsim_copy_files(DESTINATION macros FILES ${_macros})

#----------------------------------------------------------------------------
# One command line tool per analysis macro in macros/, named after the macro.
# A new analysis dropped in there gets an executable with no changes here.
fastgarsim_add_macro_apps(DIRECTORY macros
                          LINK GArAnalysis FastGArSimCommon)

#----------------------------------------------------------------------------
fastgarsim_register_component(NAME analysis
                              BIN_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              LIB_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                              LOAD_LIBS GArAnalysis
                              MACRO_DIRS ${CMAKE_CURRENT_BINARY_DIR}/macros)

install(TARGETS GArAnalysis DESTINATION lib)
install(FILES ${headers} DESTINATION include)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/macros/ DESTINATION share/GArAnalysis/macros)

#----------------------------------------------------------------------------
# Standalone build: write the environment files here
if(NOT FASTGARSIM_SUPERBUILD)
  fastgarsim_write_env()
  fastgarsim_print_summary()
endif()
