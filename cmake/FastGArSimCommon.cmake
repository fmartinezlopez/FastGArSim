#----------------------------------------------------------------------------
# FastGArSimCommon.cmake
#
# Settings and helpers shared by all FastGArSim sub-packages, so that each of
# them builds the same way whether it is configured on its own or as part of
# the top-level build.
#
# Include it from a sub-package with
#
#   list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
#   include(FastGArSimCommon)
#----------------------------------------------------------------------------

include_guard(GLOBAL)

set(FASTGARSIM_CMAKE_DIR "${CMAKE_CURRENT_LIST_DIR}" CACHE INTERNAL "")
get_filename_component(FASTGARSIM_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(FASTGARSIM_SOURCE_ROOT "${FASTGARSIM_SOURCE_ROOT}" CACHE INTERNAL "")
set(FASTGARSIM_COMMON_DIR "${FASTGARSIM_SOURCE_ROOT}/common" CACHE INTERNAL "")
set(FASTGARSIM_COMMON_INCLUDE_DIR "${FASTGARSIM_COMMON_DIR}/include" CACHE INTERNAL "")

#----------------------------------------------------------------------------
# fastgarsim_standard_settings()
#
# Compiler and build-type defaults. A macro rather than a function, so that the
# variables land in the caller's scope.
#----------------------------------------------------------------------------
macro(fastgarsim_standard_settings)
  # Default build type; without one CMake builds unoptimised code silently
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Build type: Debug, Release, RelWithDebInfo or MinSizeRel" FORCE)
  endif()

  if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
  endif()

  # The UPS/CVMFS ROOT builds hand out Clang-specific flags that GCC rejects.
  # AppleClang accepts them, so they are only stripped for GCC.
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    string(REPLACE "-Qunused-arguments" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REPLACE "-Qunused-arguments" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  endif()

  # The executables link the dictionary libraries, so they need to find them
  # both in the build tree and after `make install`. The relative-lookup token
  # is the one real difference between the two platforms; the library suffix is
  # not, because CMake and ROOT both derive it from
  # CMAKE_SHARED_LIBRARY_SUFFIX (.dylib on macOS, .so on Linux) and so agree on
  # what the .rootmap should name.
  set(CMAKE_SKIP_BUILD_RPATH OFF)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
  if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
  else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
  endif()
endmacro()

#----------------------------------------------------------------------------
# fastgarsim_use_common()
#
# Make the shared data types (and their dictionaries) available. In the
# top-level build common/ has already been added; standalone it is pulled in
# under this package's build directory.
#----------------------------------------------------------------------------
macro(fastgarsim_use_common)
  if(NOT TARGET SimDataDictLib)
    add_subdirectory(${FASTGARSIM_COMMON_DIR} ${CMAKE_CURRENT_BINARY_DIR}/common)
  endif()
endmacro()

#----------------------------------------------------------------------------
# fastgarsim_add_dictionary(<name>
#                           HEADERS <header> [<header> ...]
#                           LINKDEF <linkdef>)
#
# Generate a ROOT dictionary and wrap it in a shared library named
# lib<name>.<platform suffix>, matching the name ROOT writes into the .rootmap
# so that autoloading works. Creates the target <name>Lib.
#----------------------------------------------------------------------------
function(fastgarsim_add_dictionary name)
  cmake_parse_arguments(ARG "" "LINKDEF" "HEADERS" ${ARGN})

  if(NOT ARG_HEADERS OR NOT ARG_LINKDEF)
    message(FATAL_ERROR "fastgarsim_add_dictionary(${name}): HEADERS and LINKDEF are required")
  endif()

  # -inlineInputHeader embeds the header contents, so the dictionary stays
  # usable when the headers are not where they were at build time
  ROOT_GENERATE_DICTIONARY(${name}
                           ${ARG_HEADERS}
                           LINKDEF ${ARG_LINKDEF}
                           OPTIONS -inlineInputHeader)

  add_library(${name}Lib SHARED ${CMAKE_CURRENT_BINARY_DIR}/${name}.cxx)
  target_link_libraries(${name}Lib PUBLIC ${ROOT_LIBRARIES})
  set_target_properties(${name}Lib PROPERTIES
                        OUTPUT_NAME ${name}
                        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # The .pcm and .rootmap have to sit next to the library for ROOT to pick
  # them up, which they do because both land in CMAKE_CURRENT_BINARY_DIR.
  # rootcling names the pcm lib<name>_rdict.pcm, matching the library.
  install(TARGETS ${name}Lib DESTINATION lib)
  install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/lib${name}_rdict.pcm
          ${CMAKE_CURRENT_BINARY_DIR}/lib${name}.rootmap
          DESTINATION lib)
endfunction()

#----------------------------------------------------------------------------
# fastgarsim_copy_files(DESTINATION <dir> FILES <file> [<file> ...])
#
# Copy files into the build tree, so that the executables and ROOT can be run
# straight from there. Paths are relative to the current source directory;
# DESTINATION is relative to the current binary directory.
#----------------------------------------------------------------------------
function(fastgarsim_copy_files)
  cmake_parse_arguments(ARG "" "DESTINATION" "FILES" ${ARGN})

  foreach(_file ${ARG_FILES})
    get_filename_component(_name ${_file} NAME)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_file}
                   ${CMAKE_CURRENT_BINARY_DIR}/${ARG_DESTINATION}/${_name}
                   COPYONLY)
  endforeach()
endfunction()

#----------------------------------------------------------------------------
# _fastgarsim_macro_description(<macro file> <output variable>)
#
# Pull a one-or-two line description out of the macro's header comment, so that
# the generated tool has something useful to print for --help. Prefers the text
# under a "Description:" heading and otherwise takes the first meaningful
# comment line.
#----------------------------------------------------------------------------
function(_fastgarsim_macro_description path out)
  file(STRINGS ${path} _lines LIMIT_COUNT 40)

  get_filename_component(_name ${path} NAME)
  set(_description "")
  set(_in_description FALSE)

  foreach(_line ${_lines})
    # Strip comment markers and surrounding whitespace
    string(REGEX REPLACE "^[ \t]*(/\\*+|\\*+/|\\*|//+)" "" _text "${_line}")
    string(REGEX REPLACE "\\*+/[ \t]*$" "" _text "${_text}")
    string(STRIP "${_text}" _text)

    # Semicolons and quotes would break the generated string literal
    string(REPLACE ";" "," _text "${_text}")
    string(REPLACE "\"" "'" _text "${_text}")

    if(_in_description)
      if(_text STREQUAL "")
        break()
      endif()
      if(_description STREQUAL "")
        set(_description "${_text}")
      else()
        set(_description "${_description}\n${_text}")
      endif()
    elseif(_text MATCHES "^[Dd]escription:")
      set(_in_description TRUE)
      string(REGEX REPLACE "^[Dd]escription:[ \t]*" "" _text "${_text}")
      if(NOT _text STREQUAL "")
        set(_description "${_text}")
      endif()
    elseif(_description STREQUAL "" AND _text MATCHES "[A-Za-z].*[-|] .*[A-Za-z]")
      # A "Name.C - what it does" style first line
      set(_description "${_text}")
    endif()
  endforeach()

  if(_description STREQUAL "")
    set(_description "Command line front end for ${_name}.")
  endif()

  set(${out} "${_description}" PARENT_SCOPE)
endfunction()

#----------------------------------------------------------------------------
# fastgarsim_add_macro_apps(DIRECTORY <dir>
#                           [EXCLUDE <name> ...]
#                           [LINK <target> ...]
#                           [INCLUDE_DIRS <dir> ...])
#
# Build one command line tool per ROOT macro found in <dir>, named after the
# macro with the .C dropped. Dropping a new macro into the directory is all it
# takes to get an executable for it; nothing here needs editing.
#
# The generated main() includes the macro and forwards argv to its entry
# function, which must share the macro's file name -- the convention ROOT
# already requires for `.x Macro.C`. Macros that cannot be compiled, or that
# only make sense interactively, are named in EXCLUDE.
#----------------------------------------------------------------------------
function(fastgarsim_add_macro_apps)
  cmake_parse_arguments(ARG "" "DIRECTORY" "EXCLUDE;LINK;INCLUDE_DIRS" ${ARGN})

  if(NOT ARG_DIRECTORY)
    message(FATAL_ERROR "fastgarsim_add_macro_apps: DIRECTORY is required")
  endif()

  set(_dir ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_DIRECTORY})
  file(GLOB _macros CONFIGURE_DEPENDS ${_dir}/*.C)
  list(SORT _macros)

  set(_built "")
  foreach(_macro ${_macros})
    get_filename_component(FASTGARSIM_APP_NAME ${_macro} NAME_WE)

    if(FASTGARSIM_APP_NAME IN_LIST ARG_EXCLUDE)
      continue()
    endif()

    if(TARGET ${FASTGARSIM_APP_NAME})
      message(WARNING "fastgarsim_add_macro_apps: a target named "
                      "${FASTGARSIM_APP_NAME} already exists, skipping ${_macro}")
      continue()
    endif()

    _fastgarsim_macro_description(${_macro} FASTGARSIM_APP_DESCRIPTION)
    set(FASTGARSIM_APP_MACRO_FILE ${_macro})

    set(_generated ${CMAKE_CURRENT_BINARY_DIR}/apps/${FASTGARSIM_APP_NAME}.cc)
    configure_file(${FASTGARSIM_CMAKE_DIR}/MacroApp.cc.in ${_generated} @ONLY)

    add_executable(${FASTGARSIM_APP_NAME} ${_generated})
    target_link_libraries(${FASTGARSIM_APP_NAME} PRIVATE ${ROOT_LIBRARIES} ${ARG_LINK})
    target_include_directories(${FASTGARSIM_APP_NAME} PRIVATE
                               ${FASTGARSIM_COMMON_INCLUDE_DIR}
                               ${_dir}
                               ${ARG_INCLUDE_DIRS})
    set_target_properties(${FASTGARSIM_APP_NAME} PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    install(TARGETS ${FASTGARSIM_APP_NAME} DESTINATION bin)
    list(APPEND _built ${FASTGARSIM_APP_NAME})
  endforeach()

  if(_built)
    string(REPLACE ";" ", " _list "${_built}")
    message(STATUS "  tools from ${ARG_DIRECTORY}/: ${_list}")
  endif()
endfunction()

#----------------------------------------------------------------------------
# fastgarsim_register_component(NAME <name>
#                               [BIN_DIRS <dir> ...]      executables
#                               [LIB_DIRS <dir> ...]      shared libraries
#                               [INCLUDE_DIRS <dir> ...]  headers for ROOT/ACLiC
#                               [LOAD_LIBS <name> ...]    libraries rootlogon loads
#                               [MACRO_DIRS <dir> ...])   ROOT macros
#
# Record what a sub-package contributes to the runtime environment. The lists
# are accumulated across all enabled sub-packages and written out by
# fastgarsim_write_env().
#----------------------------------------------------------------------------
function(fastgarsim_register_component)
  cmake_parse_arguments(ARG "" "NAME"
                        "BIN_DIRS;LIB_DIRS;INCLUDE_DIRS;LOAD_LIBS;MACRO_DIRS" ${ARGN})

  set_property(GLOBAL APPEND PROPERTY FASTGARSIM_COMPONENTS   ${ARG_NAME})
  set_property(GLOBAL APPEND PROPERTY FASTGARSIM_BIN_DIRS     ${ARG_BIN_DIRS})
  set_property(GLOBAL APPEND PROPERTY FASTGARSIM_LIB_DIRS     ${ARG_LIB_DIRS})
  set_property(GLOBAL APPEND PROPERTY FASTGARSIM_INCLUDE_DIRS ${ARG_INCLUDE_DIRS})
  set_property(GLOBAL APPEND PROPERTY FASTGARSIM_LOAD_LIBS    ${ARG_LOAD_LIBS})
  set_property(GLOBAL APPEND PROPERTY FASTGARSIM_MACRO_DIRS   ${ARG_MACRO_DIRS})
endfunction()

#----------------------------------------------------------------------------
# fastgarsim_write_env()
#
# Write setup.sh and rootlogon.C into CMAKE_CURRENT_BINARY_DIR, describing
# everything registered so far:
#
#   setup.sh     -- source it to use the build from any directory
#   rootlogon.C  -- picked up automatically by ROOT started in this directory
#
# Call it once, after every sub-package has been added.
#----------------------------------------------------------------------------
function(fastgarsim_write_env)
  get_property(_components GLOBAL PROPERTY FASTGARSIM_COMPONENTS)
  get_property(_bin_dirs   GLOBAL PROPERTY FASTGARSIM_BIN_DIRS)
  get_property(_lib_dirs   GLOBAL PROPERTY FASTGARSIM_LIB_DIRS)
  get_property(_inc_dirs   GLOBAL PROPERTY FASTGARSIM_INCLUDE_DIRS)
  get_property(_load_libs  GLOBAL PROPERTY FASTGARSIM_LOAD_LIBS)
  get_property(_macro_dirs GLOBAL PROPERTY FASTGARSIM_MACRO_DIRS)

  list(REMOVE_DUPLICATES _bin_dirs)
  list(REMOVE_DUPLICATES _lib_dirs)
  list(REMOVE_DUPLICATES _inc_dirs)
  list(REMOVE_DUPLICATES _load_libs)
  list(REMOVE_DUPLICATES _macro_dirs)

  # Shell-friendly, colon-separated
  string(REPLACE ";" ":" FASTGARSIM_ENV_BIN_PATH     "${_bin_dirs}")
  string(REPLACE ";" ":" FASTGARSIM_ENV_LIB_PATH     "${_lib_dirs}")
  string(REPLACE ";" ":" FASTGARSIM_ENV_INCLUDE_PATH "${_inc_dirs}")
  string(REPLACE ";" ":" FASTGARSIM_ENV_MACRO_PATH   "${_macro_dirs}")
  string(REPLACE ";" " " FASTGARSIM_ENV_COMPONENTS   "${_components}")

  # rootlogon.C bodies
  set(FASTGARSIM_ROOTLOGON_INCLUDES "")
  foreach(_dir ${_inc_dirs})
    string(APPEND FASTGARSIM_ROOTLOGON_INCLUDES
           "    gInterpreter->AddIncludePath(\"${_dir}\");\n")
    string(APPEND FASTGARSIM_ROOTLOGON_INCLUDES
           "    gSystem->AddIncludePath(\"-I${_dir}\");\n")
  endforeach()

  set(FASTGARSIM_ROOTLOGON_PATHS "")
  foreach(_dir ${_lib_dirs})
    string(APPEND FASTGARSIM_ROOTLOGON_PATHS
           "    gSystem->AddDynamicPath(\"${_dir}\");\n")
  endforeach()

  set(FASTGARSIM_ROOTLOGON_LOADS "")
  foreach(_lib ${_load_libs})
    string(APPEND FASTGARSIM_ROOTLOGON_LOADS
           "    if (gSystem->Load(\"lib${_lib}\") < 0)\n")
    string(APPEND FASTGARSIM_ROOTLOGON_LOADS
           "        Error(\"rootlogon\", \"could not load lib${_lib} -- did you run make?\");\n")
  endforeach()

  set(FASTGARSIM_ROOTLOGON_MACROPATH "")
  foreach(_dir ${_macro_dirs})
    string(APPEND FASTGARSIM_ROOTLOGON_MACROPATH
           "    gROOT->SetMacroPath(TString::Format(\"%s:${_dir}\", gROOT->GetMacroPath()));\n")
  endforeach()

  set(FASTGARSIM_ENV_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  configure_file(${FASTGARSIM_CMAKE_DIR}/setup.sh.in
                 ${CMAKE_CURRENT_BINARY_DIR}/setup.sh @ONLY)
  configure_file(${FASTGARSIM_CMAKE_DIR}/rootlogon.C.in
                 ${CMAKE_CURRENT_BINARY_DIR}/rootlogon.C @ONLY)
  configure_file(${FASTGARSIM_CMAKE_DIR}/rootrc.in
                 ${CMAKE_CURRENT_BINARY_DIR}/.rootrc @ONLY)
endfunction()

#----------------------------------------------------------------------------
# fastgarsim_print_summary()
#
# Configuration summary, so that a build with components switched off says so.
#----------------------------------------------------------------------------
function(fastgarsim_print_summary)
  get_property(_components GLOBAL PROPERTY FASTGARSIM_COMPONENTS)
  string(REPLACE ";" ", " _components "${_components}")

  message(STATUS "")
  message(STATUS "---- FastGArSim ------------------------------------------")
  message(STATUS "  build type   : ${CMAKE_BUILD_TYPE}")
  message(STATUS "  C++ standard : ${CMAKE_CXX_STANDARD}")
  message(STATUS "  components   : ${_components}")
  message(STATUS "  environment  : source ${CMAKE_CURRENT_BINARY_DIR}/setup.sh")
  message(STATUS "----------------------------------------------------------")
  message(STATUS "")
endfunction()
