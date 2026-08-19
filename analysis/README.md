# GArAnalysis

Compiled analysis framework for FastGArSim, plus the ROOT macros built on top of it.

Everything that used to be duplicated in each analysis macro — opening files, locating trees, wiring up dozens of branches, driving the event loop, reporting progress and writing the output — now lives in a shared library. A macro derives from `ana::AnalysisBase` and implements a single method, `Run()`, which is called once per event with all input variables already loaded.

## Layout

```
analysis/
├── include/
│   ├── AnalysisEvent.hh    # GenieEvent, SimEvent, GeometryInfo -- the input trees
│   ├── AnalysisBase.hh     # AnalysisConfig and the AnalysisBase driver
│   ├── AnalysisInput.hh    # Wildcards, file lists, /pnfs -> XRootD
│   ├── AnalysisMath.hh     # Truncated mean, logarithmic binning
│   ├── Clustering.hh       # DBSCAN clustering of 3D points
│   └── PlotStyle.hh        # Shared plotting style and drawing helpers
├── src/                    # Implementations
├── macros/
│   ├── ExampleAnalysis.C   # Worked example -- copy this to start a new analysis
│   └── TruncatedDEDX.C     # Truncated mean dE/dx vs momentum
└── CMakeLists.txt
```

## Build

Built as part of the top-level FastGArSim build, or on its own — it needs ROOT but not Geant4:

```bash
cd analysis
mkdir build && cd build
cmake .. && make -j4
```

Build options and environment setup are in the [top-level README](../README.md#building).

This produces `libGArAnalysis` together with two generated helpers, in the build directory:

- `rootlogon.C` — picked up automatically when ROOT is started from the build
  directory; it adds the headers to the include path, loads the library and puts the
  macros on ROOT's macro search path.
- `setup.sh` — source it to get the same environment from any directory.

## Running a macro

From the build directory:

```bash
root -l 'macros/ExampleAnalysis.C("ntuple.root", "example_out.root")'
```

From anywhere else:

```bash
source /path/to/build/setup.sh
ExampleAnalysis ntuple.root example_out.root
```

The input is the flat ntuple produced by `EventToNtupleConverter.C` (see the [simulation README](../detector_simulation/README.md)), which contains the `AnaTree` and `GeoTree` trees. A GENIE `gst` file is optional and only needed if the analysis uses truth-level information.

Macros are ordinary ROOT macros, so they can also be compiled with ACLiC (`root -l 'ExampleAnalysis.C+("ntuple.root", "out.root")'`) once the environment is set up.

## Writing an analysis

```cpp
R__LOAD_LIBRARY(libGArAnalysis)

#include "AnalysisBase.hh"

class MyAnalysis : public ana::AnalysisBase {
protected:
    // Optional: book output branches and histograms
    void BeginJob() override {
        Output()->Branch("nECalHits", &fNECalHits);
    }

    // Required: called once per event
    void Run() override {
        fNECalHits = static_cast<Int_t>(sim.NECalHits());
        Fill();
    }

    // Optional: summaries, plots
    void EndJob() override {}

private:
    Int_t fNECalHits = 0;
};

void my_analysis(const char* simFiles, const char* outFile) {
    ana::AnalysisConfig config;
    config.simFiles = simFiles;
    config.outputFile = outFile;

    MyAnalysis analysis;
    analysis.Execute(config);
}
```

The name of the entry function must match the file name, and must differ from the name of the analysis class.

### What the base class gives you

| Member | Contents |
| --- | --- |
| `sim` | `SimEvent` — the `AnaTree` particle and hit collections, reloaded every event |
| `genie` | `GenieEvent` — the GENIE `gst` record for this event, when a gst file was supplied |
| `geo` | `GeometryInfo` — the `GeoTree` detector configuration, loaded once |

| Helper | Purpose |
| --- | --- |
| `Output()` | Output `TTree`; `nullptr` if no output file was configured |
| `OutputFile()` | Output `TFile`; histograms booked in `BeginJob()` are attached to it and written automatically |
| `Fill()` | Fill one entry of the output tree |
| `HasGenie()`, `HasGeometry()` | Whether the optional inputs are present |
| `CurrentEntry()`, `NEvents()` | Position in, and length of, the event loop |
| `Abort()` | Stop the event loop early from inside `Run()` |

`SimEvent` exposes the raw branch vectors (`ecalHitEdep`, `tpcHitX`, …) plus convenience accessors: `NParticles()`, `NECalHits()`, `StartMomentum(i)`, `ECalHitPosition(i)`, `IndexOfTrack(id)`, `ECalHitsOfTrack(id)` and `ECalEdepOfTrack(id)`. The per-track lookups are backed by maps built once per event, so they are cheap to call repeatedly.

Branches that a given file does not provide are reported once at start-up and left as `nullptr` (vectors) or at their default value (scalars), so the same reader works with gun samples, older files and files without a geometry tree. Check the pointer before dereferencing an optional branch.

### Configuration

`ana::AnalysisConfig` replaces the long list of positional `const char*` arguments:

| Field | Default | Meaning |
| --- | --- | --- |
| `simFiles` | — | FastGArSim flat ntuples (required); path, wildcard or comma-separated list |
| `outputFile` | `""` | Output file; empty means no output file is created |
| `genieFiles` | `""` | GENIE gst files; empty means no truth information |
| `xrootdForPnfs` | `kTRUE` | Stream resolved `/pnfs` paths over XRootD |
| `simTreeName` | `"AnaTree"` | Name of the analysis tree |
| `geoTreeName` | `"GeoTree"` | Name of the geometry tree |
| `genieTreeName` | `"gst"` | Name of the GENIE tree |
| `outputTreeName` | `"AnaOutput"` | Name of the output tree |
| `firstEvent` | `0` | First entry to process |
| `maxEvents` | `-1` | Number of entries to process; `-1` means all |
| `nReports` | `10` | Number of progress messages over the whole job |
| `printGeometry` | `kTRUE` | Print the detector configuration at start-up |

For quick jobs there is a shorthand:

```cpp
MyAnalysis().Execute("ntuple.root", "out.root");
```

### Input files

`simFiles` and `genieFiles` are input *specifications*, not single paths. Each is a comma-separated list whose entries may be plain paths or shell wildcard patterns; `~` and `$VARIABLES` are expanded first. The resulting files are chained, so one job can span a whole production.

```cpp
config.simFiles = "run1.root";                       // one file
config.simFiles = "gun_*.root";                      // wildcard
config.simFiles = "gun_a.root, gun_b.root";          // list
config.simFiles = "set1/*.root, set2/*.root";        // several patterns
```

Files matched by more than one pattern are read once, so an overlapping pattern cannot silently double your statistics. A pattern that matches nothing, and a named file that does not exist, are both reported.

The **geometry is taken from the first file** of the chain; if the files in a chain were produced with different detector configurations, only the first one is described by `geo`.

### Reading from dCache over XRootD

Paths under `/pnfs` are rewritten as XRootD URLs and streamed, rather than being read through the NFS mount:

```
/pnfs/dune/scratch/users/me/f.root
    -> root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/scratch/users/me/f.root
```

A named `/pnfs` file therefore does **not** need the mount to be visible. A `/pnfs` *wildcard* does, because a remote area cannot be listed — the pattern is expanded against the mount first, and only the resolved paths are converted. Entries that are already URLs are passed through untouched.

The door and dCache prefix can be changed for a different site without recompiling:

```bash
export FASTGAR_XROOTD_DOOR=root://someotherdoor.fnal.gov:1094
export FASTGAR_PNFS_PREFIX=/pnfs/fnal.gov/usr/
```

Set `config.xrootdForPnfs = kFALSE` to read through the mount instead. This needs a ROOT built with XRootD support, which the CVMFS DUNE builds have.

#### Authentication

Reading from dCache needs valid credentials. Without them the open fails with

```
Error in <TNetXNGFile::Open>: [FATAL] Auth failed: No protocols left to try
```

which means XRootD ran out of authentication mechanisms — the credentials are missing or expired, and the files themselves are fine. On a FNAL/DUNE node:

```bash
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup_fnal_security
```

or individually:

```bash
htgettoken -a htvaultprod.fnal.gov -i dune
voms-proxy-init -rfc -noregen -voms dune:/dune/Role=Analysis
```

Credentials expire, so this has to be repeated in a new session or after a few hours. Check with `voms-proxy-info -all` and `httokendecode`, and test the door on its own — which separates a credential problem from a ROOT one:

```bash
xrdfs fndca1.fnal.gov ls /pnfs/fnal.gov/usr/dune/scratch/users/$USER
```

The analysis prints this guidance itself when an input fails to open.

Chaining several files while also reading GENIE is only meaningful if `eventID` is unique across the whole chain, since the truth record is looked up by it. The job warns when it detects this situation; check your production before trusting the truth-level quantities in that case.

## Utilities

- `ana::DBSCAN3D(points, eps, minPts)` — density-based clustering of 3D points, with a
  uniform grid for fast neighbour lookup. Returns `ana::Cluster` objects that carry the
  indices of their members and can compute a centroid (optionally energy-weighted),
  a total weight and a radius.
- `ana::TruncatedMean(values, fraction, &nKept)` — mean of the lowest `fraction` of
  `values`, the standard way to tame the Landau tail of ionisation samples. Sorts
  `values` in place.
- `ana::LogSpacedBins(nBins, min, max)` — logarithmically spaced bin edges, for the
  `TH1`/`TH2` constructors that take an explicit edge array.
- `ana::SetPlotStyle()` — the standard plotting style.
- `ana::DrawComparison(...)` — draw two area-normalised histograms on one canvas and
  save it.

## Macros

### `ExampleAnalysis.C`

Worked example, and the starting point for a new analysis. Clusters the ECal barrel hits of each layer with DBSCAN and writes one output entry per cluster.

### `TruncatedDEDX.C`

Truncated mean dE/dx of primary particles in the gas TPC versus their initial momentum, for particle-gun samples.

```bash
TruncatedDEDX gun_ntuple.root dedx.root
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `inputFilesG4` | — | Flat ntuples from a gun sample; path, wildcard or comma-separated list |
| `outputFileName` | — | Output file |
| `pMin` | `-1` | Lower edge of the momentum axis, GeV/c; negative takes it from the data |
| `pMax` | `-1` | Upper edge of the momentum axis, GeV/c; negative takes it from the data |
| `dedxMin` | `-1` | Lower edge of the dE/dx axis, keV/cm; negative uses 0 |
| `dedxMax` | `-1` | Upper edge of the dE/dx axis, keV/cm; negative takes it from the data |
| `truncation` | `0.6` | Fraction of the dE/dx samples to keep (lowest first) |
| `excludeSecondaries` | `kTRUE` | Drop deposits flagged `tpcHitIsSec` |
| `requireExit` | `kFALSE` | Keep only tracks that leave the TPC (punch-through) |
| `minHits` | `10` | Minimum samples required per particle |
| `maxEvents` | `-1` | Number of events to process |

The ranges apply to the binning of the plots, not to the event selection: entries outside them land in the under/overflow bins and are excluded from the profile, but the output tree always holds every primary that passed `minHits`.

Primaries are identified as the particles with `motherID == 0`. Each of their TPC hits gives one dE/dx sample, `tpcHitEdep / tpcHitStepSize`; the lowest `truncation` of those samples are averaged. `tpcHitIsSec` marks deposits made by unstored secondaries — mostly delta rays — that the converter folds back onto the parent's track ID, so excluding them is what isolates the parent's own ionisation.

`requireExit` keeps only the tracks that leave the TPC gas volume, so that the sample is restricted to particles whose range exceeds the detector and excludes those that stop inside and deposit their whole remaining energy there. A track counts as exiting when it stops outside the gas cylinder — `r > gar_tpc_radius` or `|z| > gar_tpc_length / 2`, taken from the geometry tree and applied to the track's end point, which the simulation records wherever the particle finally stops in the detector.

The macro converts to the units conventionally used for dE/dx plots — **GeV/c** and **keV/cm** — as the values are computed, so the output tree, the plots and the `pMin`/`pMax`/`dedxMin`/`dedxMax` arguments are all in those units.

The output tree carries `momentum` (GeV/c), `dedx` (keV/cm), `pdgCode`, `nHits`, `nHitsKept` and `trackLength` (cm) per primary. The macro also writes a 2D histogram of dE/dx against momentum with logarithmic momentum bins, overlays the profile of the mean, and saves `truncated_dedx_vs_momentum.png`. Left to themselves both axes are ranged from the data, so the plot works whatever the gas pressure and gun momentum range happen to be.

The output file also contains a `dedx_slices/` folder holding the 1D dE/dx distribution of each momentum bin, as `hDEDX_pbin01` … `hDEDX_pbin60`. These are Y projections of the 2D histogram, so slice `N` is exactly bin `N` of the plot; each title carries the momentum range of its bin. One histogram is written per bin, including the empty ones at the edges of the range, so the numbering always lines up with the 2D histogram.

```cpp
TFile f("dedx.root");
f.Get<TH1D>("dedx_slices/hDEDX_pbin20")->Draw();
```
