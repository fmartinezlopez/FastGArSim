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
# Shared data types and their dictionaries. The analysis reads the simulation
# and reconstruction objects as they were written, so it needs them itself,
# along with the schema record that says what a file holds.
fastgarsim_use_common()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${FASTGARSIM_COMMON_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Sources and headers
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Shared analysis library. The analysis macros are compiled against this at
# run time and derive from ana::AnalysisBase; no ROOT dictionary of its own is
# needed, because none of the types defined here are written to a TTree -- the
# ones that are read from one come from common/, with their dictionaries.
#
# The target is named GArAnalysisLib so that the executable below can be the
# one called GArAnalysis; the library file is libGArAnalysis either way.
add_library(GArAnalysisLib SHARED ${sources} ${headers})
target_link_libraries(GArAnalysisLib PUBLIC ${ROOT_LIBRARIES} FastGArSimCommon)
target_include_directories(GArAnalysisLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(GArAnalysisLib PROPERTIES
                      OUTPUT_NAME GArAnalysis
                      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Copy the analysis macros and their job macros to the build directory. Unlike
# the other sub-packages this is not only a convenience: the analysis macro is
# compiled when the job runs, so it has to be somewhere GArAnalysis can find
# it -- see fastgarsim::FindMacro.
file(GLOB _macros RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     ${CMAKE_CURRENT_SOURCE_DIR}/macros/*.C
     ${CMAKE_CURRENT_SOURCE_DIR}/macros/*.mac)
fastgarsim_copy_files(DESTINATION macros FILES ${_macros})

#----------------------------------------------------------------------------
# The one analysis executable. It takes the analysis macro to run as an
# argument and compiles it with ACLiC at start-up, so adding an analysis is a
# matter of writing a macro -- nothing here changes, and nothing is rebuilt.
#
# It needs to be able to tell ACLiC where the headers are without the
# environment having been set up first, so the directories of this build are
# compiled into it; an installed or relocated copy is found relative to the
# executable instead.
add_executable(GArAnalysis GArAnalysis.cc)
target_link_libraries(GArAnalysis PRIVATE ${ROOT_LIBRARIES} GArAnalysisLib FastGArSimCommon)
target_compile_definitions(GArAnalysis PRIVATE
    FASTGARSIM_ANALYSIS_INCLUDE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/include"
    FASTGARSIM_COMMON_INCLUDE_DIR="${FASTGARSIM_COMMON_INCLUDE_DIR}"
    FASTGARSIM_ANALYSIS_LIB_DIR="$<TARGET_FILE_DIR:GArAnalysisLib>"
    FASTGARSIM_COMMON_LIB_DIR="$<TARGET_FILE_DIR:GArCommon>")
set_target_properties(GArAnalysis PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
fastgarsim_register_component(NAME analysis
                              BIN_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              LIB_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                              LOAD_LIBS GArAnalysis
                              MACRO_DIRS ${CMAKE_CURRENT_BINARY_DIR}/macros)

install(TARGETS GArAnalysisLib DESTINATION lib)
install(TARGETS GArAnalysis DESTINATION bin)
install(FILES ${headers} DESTINATION include)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/macros/ DESTINATION share/GArAnalysis/macros)

#----------------------------------------------------------------------------
# Standalone build: write the environment files here
if(NOT FASTGARSIM_SUPERBUILD)
  fastgarsim_write_env()
  fastgarsim_print_summary()
endif()
