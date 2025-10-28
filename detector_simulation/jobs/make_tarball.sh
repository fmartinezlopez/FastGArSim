#!/bin/bash

# Path to FastGArSim build directory
export FastGArSim=/exp/dune/app/users/jmartin4/FastGArSim/detector_simulation/build
# Path to grid setup script
export SetupScript=/exp/dune/app/users/jmartin4/FastGArSim/detector_simulation/jobs/grid_setup.sh

if [ -z "${FastGArSim}" ]; then
  echo "[ERROR]: FastGArSim is not set up, cannot tar up required binaries."
  exit 1
fi

mkdir tar_state
cp -r ${FastGArSim} ./tar_state/FastGArSim
cp ${SetupScript} ./tar_state

# Copy also any additinal files
for var in "$@"
do
    path=$(realpath "$var")
    cp ${path} ./tar_state
done

cd tar_state
tar -zcf FastGArSim.tar.gz ./*
cd ..
mv tar_state/FastGArSim.tar.gz .
rm -r tar_state
