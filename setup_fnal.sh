#!/bin/bash
#-----------------------------------------------------------------------------
# FastGArSim -- environment setup for the FNAL machines
#
# Source this, do not execute it:
#
#     source setup_fnal.sh
#
# It sets up the UPS products the build needs (CMake, Geant4, ROOT) and, if a
# build tree is present, the FastGArSim environment on top of them.
#
# The product versions can be overridden before sourcing, which is how to build
# against something other than the default stack:
#
#     export FASTGARSIM_ROOT_VERSION=v6_28_12
#     source setup_fnal.sh
#
# A build directory other than <source>/build can be given as an argument or in
# FASTGARSIM_BUILD_DIR:
#
#     source setup_fnal.sh /path/to/build
#
# Credentials for reading files from dCache are not set up here; run
# `setup_fnal_security` afterwards when they are needed.
#-----------------------------------------------------------------------------

# Resolve the directory holding this script, whichever shell sourced it
if [ -n "${BASH_SOURCE[0]}" ]; then
    _fastgarsim_this="${BASH_SOURCE[0]}"
else
    _fastgarsim_this="${(%):-%x}"
fi
FASTGARSIM_DIR="$(cd "$(dirname "${_fastgarsim_this}")" && pwd)"
export FASTGARSIM_DIR
unset _fastgarsim_this

#-----------------------------------------------------------------------------
# Product versions. Keep these in step with detector_simulation/jobs, whose
# grid jobs have to run against the same stack.
: "${FASTGARSIM_CMAKE_VERSION:=v3_27_4}"
: "${FASTGARSIM_GEANT4_VERSION:=v4_11_2_p02}"
: "${FASTGARSIM_GEANT4_QUAL:=e26:prof}"
: "${FASTGARSIM_ROOT_VERSION:=v6_28_12}"
: "${FASTGARSIM_ROOT_QUAL:=e26:p3915:prof}"

_fastgarsim_dune_setup=/cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

if [ ! -f "${_fastgarsim_dune_setup}" ]; then
    echo "[FastGArSim] ERROR: ${_fastgarsim_dune_setup} not found."
    echo "[FastGArSim]        This script is for the FNAL machines; elsewhere set up"
    echo "[FastGArSim]        CMake, Geant4 and ROOT yourself and source <build>/setup.sh."
    unset _fastgarsim_dune_setup
    return 1 2>/dev/null || exit 1
fi

echo "[FastGArSim] Setting up the DUNE environment"
source "${_fastgarsim_dune_setup}"
unset _fastgarsim_dune_setup

setup cmake  "${FASTGARSIM_CMAKE_VERSION}"                                    || return 1
setup geant4 "${FASTGARSIM_GEANT4_VERSION}" -q "${FASTGARSIM_GEANT4_QUAL}"    || return 1
setup root   "${FASTGARSIM_ROOT_VERSION}"   -q "${FASTGARSIM_ROOT_QUAL}"      || return 1

echo "[FastGArSim] cmake  ${FASTGARSIM_CMAKE_VERSION}"
echo "[FastGArSim] geant4 ${FASTGARSIM_GEANT4_VERSION} -q ${FASTGARSIM_GEANT4_QUAL}"
echo "[FastGArSim] root   ${FASTGARSIM_ROOT_VERSION} -q ${FASTGARSIM_ROOT_QUAL}"

#-----------------------------------------------------------------------------
# The build tree, if there is one. CMake writes setup.sh into the build
# directory; it puts the executables, libraries and headers on the search paths.
_fastgarsim_build="${1:-${FASTGARSIM_BUILD_DIR:-${FASTGARSIM_DIR}/build}}"

if [ -f "${_fastgarsim_build}/setup.sh" ]; then
    source "${_fastgarsim_build}/setup.sh"
else
    echo "[FastGArSim] No build found in ${_fastgarsim_build}; to build:"
    echo "[FastGArSim]     mkdir -p ${_fastgarsim_build} && cd ${_fastgarsim_build}"
    echo "[FastGArSim]     cmake ${FASTGARSIM_DIR} && make -j4"
    echo "[FastGArSim] then source this script again."
fi
unset _fastgarsim_build
