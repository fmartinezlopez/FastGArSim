#----------------------------------------------------------------------------
# Shared data types
#
# The classes written to and read from the ROOT files are defined here, so that
# the simulation, the reconstruction and the analysis all agree on them. Their
# dictionaries are generated once, in this directory: generating them in more
# than one sub-package would give ROOT two dictionaries for the same classes in
# a combined build.
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16...3.21)
project(FastGArSimCommon CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
include(FastGArSimCommon)

fastgarsim_standard_settings()

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Simulation output types (root::Event, root::Particle, root::ECalHit, ...)
fastgarsim_add_dictionary(SimDataDict
                          HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/SimDataTypes.hh
                          LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/SimDataTypesLinkDef.h)

#----------------------------------------------------------------------------
# Digitization output types (digi::TileDigiHit, digi::StripDigiHit, ...)
fastgarsim_add_dictionary(DigiDataDict
                          HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/DigiDataTypes.hh
                          LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/DigiDataTypesLinkDef.h)

#----------------------------------------------------------------------------
# Anything linking this gets the headers with it
add_library(FastGArSimCommon INTERFACE)
target_include_directories(FastGArSimCommon INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(FastGArSimCommon INTERFACE SimDataDictLib DigiDataDictLib)

fastgarsim_register_component(NAME common
                              LIB_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                              LOAD_LIBS SimDataDict DigiDataDict)

#----------------------------------------------------------------------------
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/SimDataTypes.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/DigiDataTypes.hh
        DESTINATION include)
