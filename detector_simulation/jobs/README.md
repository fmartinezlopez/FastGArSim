## Grid job submission

### Automated Configuration (Recommended)

Use the Python configuration script to automatically generate all job scripts:

1. **Create or edit a configuration file:**
   ```bash
   # Create a template (if you don't have one)
   python configure_grid.py --create-template my_config.json

   # Or copy and edit the example
   cp grid_config_example.json my_config.json
   ```

2. **Edit your configuration file** (`my_config.json`) with your settings:
   - Username and paths
   - Number of jobs and events per job
   - Output directories and file naming
   - Resource requirements (memory, disk, CPU, time)
   - Software versions

3. **Generate all scripts:**
   ```bash
   python configure_grid.py -c my_config.json
   ```

   This will automatically configure:
   - `grid_setup.sh` - Software environment setup
   - `make_tarball.sh` - Build directory paths
   - `run_grid.sh` - Output configuration and processing
   - `launch_job.sh` - Job submission parameters

4. **Create the tarball:**
   ```bash
   bash make_tarball.sh
   ```

5. **From an Alma9 session**, set up the DUNE environment and load `fife-utils`:
   ```bash
   source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

   spack load fife-utils@3.7.0 # ifdh
   export IFDH_TOKEN_ENABLE=1
   export IFDH_PROXY_ENABLE=0
   ```

6. **Submit jobs:**
   ```bash
   bash launch_job.sh
   ```

### Manual Configuration (Legacy)

Modify `make_tarball.sh` with your complete path to the build directory and the `grid_setup.sh` script. Create tarball running:
```bash
bash make_tarball.sh
```

Depending on what you are doing `run_grid.sh` has to be modified accordingly. At minimum change the output directory and name, given by `OUTDIR` and `OUTFILE`.

The default script and macros process an input file in batches of `1000` events, using the jobid in `PROCESS` to shift the first event index in the job.

From an Alma9 session, set up the DUNE environment and load `fife-utils`:
```bash
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

spack load fife-utils@3.7.0 # ifdh
export IFDH_TOKEN_ENABLE=1
export IFDH_PROXY_ENABLE=0
```

Then, launch the grid job by executing:
```bash
bash launch_job.sh
```
Modify the number of jobs submitted and the other job submision parameters as needed.