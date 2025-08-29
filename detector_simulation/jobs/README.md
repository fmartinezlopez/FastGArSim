## Grid job submission

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