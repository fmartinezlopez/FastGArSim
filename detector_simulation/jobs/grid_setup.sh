#!/bin/bash

# we cannot rely on "whoami" in a grid job. We have no idea what the local username will be. 
# Use the GRID_USER environment variable instead (set automatically by jobsub).
USERNAME=${GRID_USER}

export WORKDIR=${_CONDOR_JOB_IWD} # if we use the RCDS the our tarball will be placed in $INPUT_TAR_DIR_LOCAL.
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo .`
fi

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup ifdhc
setup cmake  v3_27_4
setup geant4 v4_11_2_p02 -q e26:prof
setup root   v6_28_12    -q e26:p3915:prof
