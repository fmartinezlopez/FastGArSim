#----------------------------------------------------------------------------
# Shared data types and file-level utilities
#
# The classes written to and read from the ROOT files are defined here, so that
# the simulation, the reconstruction and the analysis all agree on them. Their
# dictionaries are generated once, in this directory: generating them in more
# than one sub-package would give ROOT two dictionaries for the same classes in
# a combined build.
#
# This directory also holds the code that works on a file without knowing what
# is in it -- the schema record and the generic flattener -- and the command
# line tools built on them.
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
# Utilities that work on a file without knowing its contents. No dictionary of
# their own: nothing here is written to a TTree as an object, and the schema
# trees are plain std::string branches on purpose.
add_library(GArCommon SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/src/ProductSchema.cc
            ${CMAKE_CURRENT_SOURCE_DIR}/src/TreeFlattener.cc)
target_link_libraries(GArCommon PUBLIC ${ROOT_LIBRARIES})
target_include_directories(GArCommon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(GArCommon PROPERTIES
                      OUTPUT_NAME GArCommon
                      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Anything linking this gets the headers and the dictionaries with it
add_library(FastGArSimCommon INTERFACE)
target_include_directories(FastGArSimCommon INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(FastGArSimCommon INTERFACE SimDataDictLib DigiDataDictLib GArCommon)

#----------------------------------------------------------------------------
# Copy the utility macros next to the library, so ROOT can run them from the
# build directory, and build one command line tool per macro
fastgarsim_copy_files(DESTINATION .
                      FILES utils/MakeNtuple.C
                            utils/DumpSchema.C)

fastgarsim_add_macro_apps(DIRECTORY utils
                          LINK GArCommon SimDataDictLib DigiDataDictLib)

#----------------------------------------------------------------------------
fastgarsim_register_component(NAME common
                              BIN_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              LIB_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                              INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                              LOAD_LIBS SimDataDict DigiDataDict GArCommon
                              MACRO_DIRS ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
install(TARGETS GArCommon DESTINATION lib)
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/SimDataTypes.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/DigiDataTypes.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/ProductSchema.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/TreeFlattener.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/AppOptions.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/MacroPath.hh
        DESTINATION include)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/utils/ DESTINATION share/FastGArSim/utils)
