#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(GArSimulation)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT package with required components
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Generate ROOT dictionary for ROOTDataTypes (without Geant4 dependencies)
ROOT_GENERATE_DICTIONARY(ROOTDataDict
                         ${PROJECT_SOURCE_DIR}/include/ROOTDataTypes.hh
                         LINKDEF ${PROJECT_SOURCE_DIR}/include/ROOTDataTypesLinkDef.h
                         OPTIONS -inlineInputHeader)  # This embeds header content

#----------------------------------------------------------------------------
# Create a shared library with the dictionary
add_library(ROOTDataDictLib SHARED ROOTDataDict.cxx)
target_link_libraries(ROOTDataDictLib ${ROOT_LIBRARIES})
set_target_properties(ROOTDataDictLib PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
                      SUFFIX ".so"
                      OUTPUT_NAME "ROOTDataDict")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
add_executable(GArSimulation GArSimulation.cc ${sources} ${headers} ROOTDataDict.cxx)
target_link_libraries(GArSimulation ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build GArSimulation. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
set(SCRIPTS
  macros/init.mac
  macros/init_base.mac
  macros/init_gun.mac
  macros/run.mac
  macros/nu.mac
  macros/gun.mac
  macros/vis.mac
  )

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Copy also all macros to the build directory
set(MACROS
  utils/EventToNtupleConverter.C
  utils/GeoVis.C
  )

foreach(_macro ${MACROS})
  get_filename_component(_name ${_macro} NAME)
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_macro}
    ${PROJECT_BINARY_DIR}/${_name}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS GArSimulation DESTINATION bin)

#----------------------------------------------------------------------------
# Install dictionary files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/ROOTDataDict_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/libROOTDataDict.rootmap
  DESTINATION bin)

#----------------------------------------------------------------------------
# Install the shared library
install(TARGETS ROOTDataDictLib DESTINATION lib)
