#!/bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup cmake  v3_27_4
setup geant4 v4_11_2_p02 -q e26:prof
setup root   v6_28_12    -q e26:p3915:prof