#!/bin/bash

echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"

# Set the output location for copyback
OUTDIR=/pnfs/dune/scratch/users/jmartin4/FastGArSim

# Let's rename the output file so it's unique in case we send multiple jobs
OUTFILE=test_${CLUSTER}_${PROCESS}_$(date -u +%Y%m%d).root

# Make sure we see what we expect
pwd

ls -l $CONDOR_DIR_INPUT

# cd back to the top-level directory since we know that's writable
cd ${_CONDOR_JOB_IWD}

if [ -e ${INPUT_TAR_DIR_LOCAL}/grid_setup.sh ]; then
    . ${INPUT_TAR_DIR_LOCAL}/grid_setup.sh
else
  echo "Error, setup script not found. Exiting."
  exit 1
fi

# Set some other very useful environment variables for xrootd and IFDH
export IFDH_CP_MAXRETRIES=2
export XRD_CONNECTIONRETRY=32
export XRD_REQUESTTIMEOUT=14400
export XRD_REDIRECTLIMIT=255
export XRD_LOADBALANCERTTL=7200
export XRD_STREAMTIMEOUT=14400 # many vary for your job/file type

# Make sure the output directory exists
#ifdh ls $OUTDIR 0 # set recursion depth to 0 since we are only checking for the directory; we don't care about the full listing.
gfal-stat $OUTDIR

if [ $? -ne 0 ]; then
    # if ifdh ls failed, try to make the directory
    ifdh mkdir_p $OUTDIR || { echo "Error creating or checking $OUTDIR"; exit 2; }
fi

# Copy directory to current working directory and cd
cp -r ${INPUT_TAR_DIR_LOCAL}/FastGArSim ./FastGArSim
cd FastGArSim

# Set value of first event generated according to process number
first=$((PROCESS * 1000))
sed "s/first_event/$first/" macros/init_base.mac > macros/init.mac
echo "Starting processing at event $first."

# Run detector simulation
./GArSimulation -m macros/nu.mac

G4_RESULT=$?
if [ $G4_RESULT -ne 0 ]; then
    echo "Geant4 exited with abnormal status $G4_RESULT. See error outputs."
    exit $G4_RESULT
fi

# Run ntuple maker
echo "Running ntuple maker"
root -q -b 'EventToNtupleConverter.C("test.root", "out.root")'

NTUPLE_RESULT=$?
if [ $NTUPLE_RESULT -ne 0 ]; then
    echo "ROOT exited with abnormal status $NTUPLE_RESULT. See error outputs."
    exit $NTUPLE_RESULT
fi

if [ -f out.root ]; then

    mv out.root $OUTFILE
    
    # Copy output file back
    ifdh cp -D $OUTFILE $OUTDIR

    # Check the exit status to see if the copyback actually worked. Print a message if it did not
    IFDH_RESULT=$?
    if [ $IFDH_RESULT -ne 0 ]; then
	echo "Error during output copyback. See output logs."
	exit $IFDH_RESULT
    fi
fi

# If we got this far, we succeeded
echo "Completed successfully."
exit 0
