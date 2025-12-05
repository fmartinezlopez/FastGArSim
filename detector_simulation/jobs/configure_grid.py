#!/usr/bin/env python3
"""
Automated configuration tool for grid job submission scripts.

This script reads a configuration file and automatically generates/updates
the grid job submission scripts (grid_setup.sh, launch_job.sh, run_grid.sh, make_tarball.sh).
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Any


class GridJobConfigurator:
    """Configures grid job scripts based on user configuration."""

    def __init__(self, config_file: str, jobs_dir: str = None):
        """
        Initialize the configurator.

        Args:
            config_file: Path to JSON configuration file
            jobs_dir: Directory containing job scripts (defaults to script location)
        """
        self.config_file = Path(config_file)
        self.jobs_dir = Path(jobs_dir) if jobs_dir else Path(__file__).parent
        self.config = self._load_config()

    def _load_config(self) -> Dict[str, Any]:
        """Load and validate configuration file."""
        if not self.config_file.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_file}")

        with open(self.config_file, 'r') as f:
            config = json.load(f)

        # Validate required fields
        self._validate_config(config)
        return config

    def _validate_config(self, config: Dict[str, Any]):
        """Validate that required configuration fields are present."""
        required_fields = {
            'user': ['username', 'base_path'],
            'paths': ['build_dir'],
            'job_params': ['num_jobs', 'events_per_job'],
            'output': ['scratch_dir', 'output_prefix']
        }

        for section, fields in required_fields.items():
            if section not in config:
                raise ValueError(f"Missing required section: {section}")
            for field in fields:
                if field not in config[section]:
                    raise ValueError(f"Missing required field: {section}.{field}")

    def configure_all(self):
        """Configure all grid job scripts."""
        print("Configuring grid job scripts...")
        self.configure_grid_setup()
        self.configure_make_tarball()
        self.configure_run_grid()
        self.configure_launch_job()
        print("\n✓ All scripts configured successfully!")
        print(f"\nNext steps:")
        print(f"  1. Review the generated scripts in {self.jobs_dir}")
        print(f"  2. Create tarball: bash make_tarball.sh")
        print(f"  3. Submit jobs: bash launch_job.sh")

    def configure_grid_setup(self):
        """Generate grid_setup.sh with software versions."""
        software = self.config.get('software', {})

        script = f"""#!/bin/bash

# we cannot rely on "whoami" in a grid job. We have no idea what the local username will be.
# Use the GRID_USER environment variable instead (set automatically by jobsub).
USERNAME=${{GRID_USER}}

export WORKDIR=${{_CONDOR_JOB_IWD}} # if we use the RCDS the our tarball will be placed in $INPUT_TAR_DIR_LOCAL.
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo .`
fi

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup ifdhc
setup cmake  {software.get('cmake', 'v3_27_4')}
setup geant4 {software.get('geant4', 'v4_11_2_p02')} -q {software.get('geant4_qual', 'e26:prof')}
setup root   {software.get('root', 'v6_28_12')}    -q {software.get('root_qual', 'e26:p3915:prof')}
"""

        output_path = self.jobs_dir / 'grid_setup.sh'
        with open(output_path, 'w') as f:
            f.write(script)
        os.chmod(output_path, 0o755)
        print(f"✓ Generated {output_path}")

    def configure_make_tarball(self):
        """Generate make_tarball.sh with correct paths."""
        user = self.config['user']
        paths = self.config['paths']

        # Construct full paths
        build_path = f"{user['base_path']}/{paths['build_dir']}"
        setup_script = f"{user['base_path']}/{paths.get('jobs_subdir', 'detector_simulation/jobs')}/grid_setup.sh"

        script = f"""#!/bin/bash

# Path to FastGArSim build directory
export FastGArSim={build_path}
# Path to grid setup script
export SetupScript={setup_script}

if [ -z "${{FastGArSim}}" ]; then
  echo "[ERROR]: FastGArSim is not set up, cannot tar up required binaries."
  exit 1
fi

mkdir tar_state
cp -r ${{FastGArSim}} ./tar_state/FastGArSim
cp ${{SetupScript}} ./tar_state

# Copy also any additinal files
for var in "$@"
do
    path=$(realpath "$var")
    cp ${{path}} ./tar_state
done

cd tar_state
tar -zcf FastGArSim.tar.gz ./*
cd ..
mv tar_state/FastGArSim.tar.gz .
rm -r tar_state
"""

        output_path = self.jobs_dir / 'make_tarball.sh'
        with open(output_path, 'w') as f:
            f.write(script)
        os.chmod(output_path, 0o755)
        print(f"✓ Generated {output_path}")

    def configure_run_grid(self):
        """Generate run_grid.sh with output configuration."""
        user = self.config['user']
        output = self.config['output']
        job_params = self.config['job_params']

        outdir = f"/pnfs/dune/scratch/users/{user['username']}/{output['scratch_dir']}"
        output_prefix = output['output_prefix']
        events_per_job = job_params['events_per_job']

        script = f"""#!/bin/bash

echo "Running on $(hostname) at ${{GLIDEIN_Site}}. GLIDEIN_DUNESite = ${{GLIDEIN_DUNESite}}"

# Set the output location for copyback
OUTDIR={outdir}

# Let's rename the output file so it's unique in case we send multiple jobs
OUTFILE={output_prefix}_${{CLUSTER}}_${{PROCESS}}_$(date -u +%Y%m%d).root

# Make sure we see what we expect
pwd

ls -l $CONDOR_DIR_INPUT

# cd back to the top-level directory since we know that's writable
cd ${{_CONDOR_JOB_IWD}}

if [ -e ${{INPUT_TAR_DIR_LOCAL}}/grid_setup.sh ]; then
    . ${{INPUT_TAR_DIR_LOCAL}}/grid_setup.sh
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
    ifdh mkdir_p $OUTDIR || {{ echo "Error creating or checking $OUTDIR"; exit 2; }}
fi

# Copy directory to current working directory and cd
cp -r ${{INPUT_TAR_DIR_LOCAL}}/FastGArSim ./FastGArSim
cd FastGArSim

# Set value of first event generated according to process number
first=$((PROCESS * {events_per_job}))
sed "s/first_event/$first/" macros/init_base.mac > macros/init.mac
echo "Starting processing at event $first."

# Run detector simulation
./GArSimulation -m macros/run.mac

G4_RESULT=$?
if [ $G4_RESULT -ne 0 ]; then
    echo "Geant4 exited with abnormal status $G4_RESULT. See error outputs."
    exit $G4_RESULT
fi

# Run ntuple maker
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
\techo "Error during output copyback. See output logs."
\texit $IFDH_RESULT
    fi
fi

# If we got this far, we succeeded
echo "Completed successfully."
exit 0
"""

        output_path = self.jobs_dir / 'run_grid.sh'
        with open(output_path, 'w') as f:
            f.write(script)
        os.chmod(output_path, 0o755)
        print(f"✓ Generated {output_path}")

    def configure_launch_job(self):
        """Generate launch_job.sh with job submission parameters."""
        user = self.config['user']
        job_params = self.config['job_params']
        resources = self.config.get('resources', {})
        paths = self.config['paths']

        # Resource defaults
        num_jobs = job_params['num_jobs']
        memory = resources.get('memory', '2500MB')
        disk = resources.get('disk', '2GB')
        lifetime = resources.get('expected_lifetime', '4h')
        cpu = resources.get('cpu', 1)

        # Construct paths
        jobs_dir = f"{user['base_path']}/{paths.get('jobs_subdir', 'detector_simulation/jobs')}"
        tarball_path = f"dropbox://{jobs_dir}/FastGArSim.tar.gz"
        run_script = f"file://{jobs_dir}/run_grid.sh"

        # Singularity image
        singularity = resources.get('singularity_image',
                                    '/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest')

        script = f"""#!/bin/bash

jobsub_submit -G dune -N {num_jobs} --memory={memory} --disk={disk} --expected-lifetime={lifetime} --cpu={cpu} --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --tar_file_name={tarball_path} -l '+SingularityImage="{singularity}"'  --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' -e GFAL_PLUGIN_DIR=/usr/lib64/gfal2-plugins -e GFAL_CONFIG_DIR=/etc/gfal2.d -e UPS_OVERRIDE="-H Linux64bit+3.10-2.17" {run_script}

rm *.tbz2
"""

        output_path = self.jobs_dir / 'launch_job.sh'
        with open(output_path, 'w') as f:
            f.write(script)
        os.chmod(output_path, 0o755)
        print(f"✓ Generated {output_path}")

    def create_template_config(self, output_path: str):
        """Create a template configuration file."""
        template = {
            "_comment": "Configuration file for grid job submission",
            "user": {
                "username": "your_username",
                "base_path": "/exp/dune/app/users/your_username/path/to/FastGArSim"
            },
            "paths": {
                "build_dir": "detector_simulation/build",
                "jobs_subdir": "detector_simulation/jobs"
            },
            "job_params": {
                "num_jobs": 10,
                "events_per_job": 1000
            },
            "output": {
                "scratch_dir": "geant_test",
                "output_prefix": "test"
            },
            "resources": {
                "memory": "2500MB",
                "disk": "2GB",
                "expected_lifetime": "4h",
                "cpu": 1,
                "singularity_image": "/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest"
            },
            "software": {
                "cmake": "v3_27_4",
                "geant4": "v4_11_2_p02",
                "geant4_qual": "e26:prof",
                "root": "v6_28_12",
                "root_qual": "e26:p3915:prof"
            }
        }

        with open(output_path, 'w') as f:
            json.dump(template, f, indent=2)
        print(f"✓ Created template configuration: {output_path}")
        print(f"\nEdit this file with your settings, then run:")
        print(f"  python configure_grid.py -c {output_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Configure grid job submission scripts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create a template configuration file
  python configure_grid.py --create-template my_config.json

  # Configure scripts using a config file
  python configure_grid.py -c my_config.json

  # Specify a different jobs directory
  python configure_grid.py -c my_config.json -d /path/to/jobs
        """
    )

    parser.add_argument('-c', '--config',
                       help='Path to configuration JSON file')
    parser.add_argument('-d', '--jobs-dir',
                       help='Directory containing job scripts (default: script location)')
    parser.add_argument('--create-template',
                       metavar='OUTPUT',
                       help='Create a template configuration file')

    args = parser.parse_args()

    try:
        if args.create_template:
            configurator = GridJobConfigurator.__new__(GridJobConfigurator)
            configurator.create_template_config(args.create_template)
        elif args.config:
            configurator = GridJobConfigurator(args.config, args.jobs_dir)
            configurator.configure_all()
        else:
            parser.print_help()
            print("\nError: Either --config or --create-template is required")
            sys.exit(1)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
