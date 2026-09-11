# GArAnalysis

Analysis framework for FastGArSim: one program, `GArAnalysis`, and the ROOT macros it runs.

An analysis is a ROOT macro holding a class that derives from `ana::AnalysisBase` and implements one method, `Run()`, called once per event with every input variable already loaded. The base class owns the rest of the job — opening the files, locating the trees, wiring up the branches, driving the event loop, reporting progress and writing the output — and lives in `libGArAnalysis`.

`GArAnalysis` takes the macro as an argument and compiles it when the job starts, so an analysis is a file and nothing else, and every run uses the macro as it stands at that moment. What to run it over goes on the command line; what the analysis is to make of it goes in a job macro, the same kind of file that configures `GArReconstruction`:

```bash
GArAnalysis -a TruncatedDEDX.C -i 'gun_*.root' -o dedx.root -m TruncatedDEDX.mac
```

The input is read as it was written — the simulation's `root::Event` objects and the reconstruction's product collections, straight out of the file. The reconstruction is modular, so which of its products a file holds varies; an analysis names the ones it needs and is told at start-up if they are not there.

## Layout

```
analysis/
├── GArAnalysis.cc          # The one executable: compiles a macro and runs it
├── include/
│   ├── AnalysisEvent.hh    # GenieEvent, SimEvent, GeometryInfo -- the input trees
│   ├── ProductStore.hh     # Reconstruction products, asked for by name
│   ├── AnalysisBase.hh     # AnalysisConfig, the AnalysisBase driver, ANA_ANALYSIS
│   ├── ParameterSet.hh     # The analysis's own parameters, read in Configure()
│   ├── AnalysisMacro.hh    # Job macro parsing and the run-time compilation
│   ├── AnalysisInput.hh    # Wildcards, file lists, /pnfs -> XRootD
│   ├── AnalysisMath.hh     # Truncated mean, logarithmic binning
│   ├── Clustering.hh       # DBSCAN clustering of 3D points
│   └── PlotStyle.hh        # Shared plotting style and drawing helpers
├── src/                    # Implementations
├── macros/                 # One .C per analysis, with its .mac beside it
│   ├── ExampleAnalysis.C   # Worked example -- copy this to start a new analysis
│   ├── TruncatedDEDX.C     # Truncated mean dE/dx vs momentum
│   └── TrackingEfficiency.C # TPC tracking efficiency of primaries
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

This produces the `GArAnalysis` executable and `libGArAnalysis`, together with two generated helpers, in the build directory:

- `rootlogon.C` — picked up automatically when ROOT is started from the build
  directory; it adds the headers to the include path, loads the library and puts the
  macros on ROOT's macro search path.
- `setup.sh` — source it to get the same environment from any directory.

## Running an analysis

```bash
GArAnalysis -a ExampleAnalysis.C -i tpc_reco.root -o example_out.root
```

From the build directory that works as it stands; from anywhere else, `source /path/to/build/setup.sh` first, or give the macros by their full path.

| Option | Meaning |
| --- | --- |
| `-a <file.C>` | Analysis macro to run. May also be given as the first argument, without the flag, or by the job macro |
| `-i <spec>` | Input file(s): a path, a wildcard pattern or a comma-separated list of either |
| `-o <file>` | Output ROOT file; leave it out for an analysis that only prints or draws |
| `-m <file.mac>` | Job macro holding the analysis parameters |
| `-g <spec>` | GENIE gst file(s) with the truth record |
| `-n <N>` | Events to process; `-1` (the default) is all of them |
| `-h`, `--help` | Show the usage message |

Both macros are looked for in the working directory, then in `FASTGARSIM_MACRO_PATH` (exported by `setup.sh`) and then next to the executable, so `-a TruncatedDEDX.C` works from anywhere once the environment is set up. Where the command line and the job macro say the same thing, the command line wins, which is what makes one job macro usable across a production:

```bash
for file in gun_*.root; do
    GArAnalysis -m TruncatedDEDX.mac -i "$file" -o "dedx_${file}"
done
```

The input is a FastGArSim ROOT file: either simulation output (`Events` and `Geometry` trees) or reconstruction output, which by default carries those trees along with its own `Reco` tree. A GENIE `gst` file is optional, and needed only by an analysis that uses truth-level information.

Anything the analysis draws is written to file rather than to a window, since a batch job has no display; set `FASTGARSIM_NO_BATCH` to override that.

### The compilation step

The macro is compiled with ACLiC when the job starts, which takes a few seconds and is reported as it happens:

```
 Compiling macros/TruncatedDEDX.C ...
   compiled in 5.3 s
```

It is compiled **afresh every run**, into a private build directory under `$TMPDIR` that is removed when the job ends. Nothing is cached and nothing is written next to the macro, so an edit between two runs always takes effect, and an analysis can live in a read-only directory — an installed build, or a shared production area.

Because ACLiC builds a ROOT dictionary for it, a macro is subject to the usual rules for compiled ROOT code. The one that comes up in practice: a `std::vector` member cannot be built on a type nested inside the analysis class, because the generated dictionary has no access to it. Define such a type at file scope instead, as `TrackingEfficiency.C` does with `EffPlot`.

A macro can also be compiled by hand in an interactive session, which is the quickest way to check one over:

```bash
root -l          # from the build directory, so rootlogon.C loads the library
root [0] .L macros/TruncatedDEDX.C+
root [1] ana::AnalysisConfig cfg; ana::ParseJobMacro("macros/TruncatedDEDX.mac", cfg);
root [2] cfg.inputFiles = "gun.root"; cfg.outputFile = "dedx.root";
root [3] ana_MakeAnalysis()->Execute(cfg);
```

## The job macro

Everything an analysis takes beyond its input and output is set in a job macro, read the same way as the reconstruction's: one command per line, `#` starts a comment, blank lines are ignored.

```
# TruncatedDEDX.mac
/ana/global/analysis      TruncatedDEDX.C
/ana/global/outputTree    DEDXTree
/ana/global/maxEvents     1000

/ana/truncation           0.6
/ana/excludeSecondaries   true
/ana/minHits              10
```

Commands under `/ana/global/` configure the job itself and go into `ana::AnalysisConfig`; their names are fixed:

| Command | Field | Meaning |
| --- | --- | --- |
| `/ana/global/analysis` | — | The analysis macro to run, so that `-m` alone describes a whole job |
| `/ana/global/input` | `inputFiles` | Same as `-i` |
| `/ana/global/output` | `outputFile` | Same as `-o` |
| `/ana/global/genie` | `genieFiles` | Same as `-g` |
| `/ana/global/outputTree` | `outputTreeName` | Name of the tree in the output file |
| `/ana/global/maxEvents` | `maxEvents` | Same as `-n` |
| `/ana/global/firstEvent` | `firstEvent` | First entry to process |
| `/ana/global/nReports` | `nReports` | Progress messages over the whole job |
| `/ana/global/requireSim` | `requireSim` | Stop if the input has no simulation tree |
| `/ana/global/printProducts` | `printProducts` | List the reconstruction products at start-up |
| `/ana/global/printGeometry` | `printGeometry` | Print the detector configuration at start-up |
| `/ana/global/xrootdForPnfs` | `xrootdForPnfs` | Stream resolved `/pnfs` paths over XRootD |
| `/ana/global/simTreeName` | `simTreeName` | Name of the simulation tree |
| `/ana/global/recoTreeName` | `recoTreeName` | Name of the reconstruction tree |
| `/ana/global/geoTreeName` | `geoTreeName` | Name of the geometry tree |
| `/ana/global/genieTreeName` | `genieTreeName` | Name of the GENIE tree |
| `/ana/global/simBranchName` | `simBranchName` | Branch of the simulation tree holding the event |

Everything else under `/ana/` belongs to the analysis, which reads it in `Configure()`. Those names are whatever that analysis chose to read, and are documented in its own header comment and in its `.mac`.

Nothing checks a parameter name at compile time, so both kinds of mistake are caught when the job starts, before anything is opened:

- a value that cannot be read as the type the analysis asked for (`/ana/minHits ten`) stops the job;
- a name no analysis ever asked for (`/ana/truncatoin 0.6`) is reported as ignored, which is what catches a misspelling. It is a warning rather than an error, since one macro may be shared between several analyses;
- a `/ana/global/` name that does not exist stops the job and the known names are listed.

A parameter the macro does not mention keeps the default its member is declared with, so a job macro only has to say what is to differ — and `-m` can be left out entirely for a run at the defaults.

## Writing an analysis

```cpp
#include "AnalysisBase.hh"

class MyAnalysis : public ana::AnalysisBase {
protected:
    // Optional: read this analysis's parameters from the job macro
    void Configure(const ana::ParameterSet& params) override {
        params.Get("minHits", fMinHits);     // /ana/minHits
    }

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
    Int_t fMinHits = 10;
    Int_t fNECalHits = 0;
};

ANA_ANALYSIS(MyAnalysis)
```

Save that as `MyAnalysis.C` in `macros/` and it is ready to run:

```bash
GArAnalysis -a MyAnalysis.C -i tpc_reco.root -o my_out.root
```

`ANA_ANALYSIS()` is the whole of the plumbing: it names the class `GArAnalysis` is to instantiate. That class has to be default-constructible, since its parameters arrive through `Configure()` rather than through a constructor. Which files to read and where to write them are the job's business, and reach the analysis through the base class.

### Parameters

`Configure()` is called once, before any file is opened, with whatever the job macro set under `/ana/`:

```cpp
void Configure(const ana::ParameterSet& params) override {
    params.Get("minHits",   fMinHits);      // /ana/minHits
    params.Get("minLength", fMinLength);    // /ana/minLength
    params.Get("useTruth",  fUseTruth);     // /ana/useTruth

    if (fMinHits < 1) {
        std::cerr << "MyAnalysis: /ana/minHits has to be at least 1" << std::endl;
        Abort();                            // ends the job before it opens anything
    }
}
```

`Get()` converts the text to whatever type the member has — `std::string`, `Bool_t`, `Int_t`, `Long64_t`, `Float_t` or `Double_t` — and **leaves the member alone when the macro did not set that name**, so the value it is declared with is its default and there is no second place where the defaults have to be kept in step. A value that cannot be read as that type stops the job and says so.

Calling `Get()` is also what declares the name: anything in the macro that no `Get()` ever asked for is reported at start-up, which is what catches a misspelling. For a parameter used in an expression rather than stored in a member there are `GetInt()`, `GetDouble()`, `GetBool()`, `GetLong()` and `GetString()`, each taking the fallback to use.

Document the names the analysis reads in its header comment, and ship a `.mac` beside it listing them with their defaults — as the three macros here do.

### What the base class gives you

| Member | Contents |
| --- | --- |
| `sim` | `SimEvent` — the `Events` tree, reloaded every event |
| `reco` | `ProductStore` — the reconstruction products this file holds, if any |
| `genie` | `GenieEvent` — the GENIE `gst` record for this event, when a gst file was supplied |
| `geo` | `GeometryInfo` — the `Geometry` detector configuration, loaded once |

| Helper | Purpose |
| --- | --- |
| `Params()` | The parameters this job was configured with, outside `Configure()` |
| `Output()` | Output `TTree`; `nullptr` if no output file was configured |
| `OutputFile()` | Output `TFile`; histograms booked in `BeginJob()` are attached to it and written automatically |
| `Fill()` | Fill one entry of the output tree |
| `Require<T>(name)`, `Optional<T>(name)` | Ask for a reconstruction product; call from `BeginJob()` |
| `HasGenie()`, `HasGeometry()`, `HasReco()` | Whether the optional inputs are present |
| `Schema()` | What the first input file says it holds, and how it was made |
| `CurrentEntry()`, `NEvents()` | Position in, and length of, the event loop |
| `Abort()` | Stop the event loop early from inside `Run()`, or turn the job down from inside `Configure()` |

### Reading the simulation

`sim.event` is the `root::Event` itself, so anything the simulation wrote is reachable. On top of it `SimEvent` offers the two views an analysis usually wants:

- **particles by index** — `NParticles()`, `Particle(i)`, `TrackID(i)`, `PdgCode(i)`, `MotherID(i)`, `CreatorProcess(i)`, `StartPosition(i)`, `EndPosition(i)`, `StartMomentum(i)`, `EndMomentum(i)`, `HasTrajectory(i)`
- **every hit of a sub-detector in one flat list**, regardless of which particle produced it — `NECalHits()`, `ECalHit(k)`, `ECalHitPosition(k)`, `ECalHitTrackID(k)`, `ECalHitIsSecondary(k)`, and the same for `TPC` and `MuID`

with `IndexOfTrack(id)`, `ECalHitsOfTrack(id)` and `ECalEdepOfTrack(id)` connecting the two. The flat list and the per-track maps are built lazily, once per event and only if something asks for them, so an analysis that only walks particles pays nothing for them.

### Reading reconstruction products

Which products a file holds depends on which reconstruction modules were run, so they are asked for by name in `BeginJob()`:

```cpp
class MyAnalysis : public ana::AnalysisBase {
protected:
    void BeginJob() override {
        fClusters  = Require<std::vector<digi::TPCCluster>>("TPCClusters");
        fWaveforms = Optional<std::vector<digi::TPCWaveform>>("TPCWaveforms");
    }

    void Run() override {
        for (const digi::TPCCluster& cluster : *fClusters) { /* ... */ }

        if (fWaveforms) { /* only when the file has them */ }
    }

private:
    ana::Handle<std::vector<digi::TPCCluster>> fClusters;
    ana::Handle<std::vector<digi::TPCWaveform>> fWaveforms;
};
```

`Require()` ends the job before the first event if the product is absent, or present with a different type, saying what the file does hold and which module produced it. `Optional()` gives back a handle that is simply never valid, so one `if` is all the analysis needs to run on files with and without it. A handle also goes invalid on its own if the chain moves on to a file lacking that branch, so a mixed set of inputs gives "not there" rather than the previous file's contents.

`DumpSchema file.root` prints the same information from the shell, which is the quickest way to find out what a file offers before writing anything.

### Missing inputs

A file without a `Geometry` tree leaves `geo` at its default values, with a warning; a file without a `Reco` tree leaves every product handle invalid. Either is reported once at start-up. Passing a file with no `Events` tree is an error unless `AnalysisConfig::requireSim` is set false, which is for analyses that run on reconstruction products alone.

### Configuration

`ana::AnalysisConfig` describes the job itself. `GArAnalysis` fills it in from the command line and the `/ana/global/` commands of the job macro, which map straight onto its fields; a macro compiled by hand in ROOT builds one itself. An analysis can read it through `Config()`, though it rarely needs to.

| Field | Default | Meaning |
| --- | --- | --- |
| `inputFiles` | — | FastGArSim simulation or reconstruction files (required); path, wildcard or comma-separated list |
| `outputFile` | `""` | Output file; empty means no output file is created |
| `genieFiles` | `""` | GENIE gst files; empty means no truth information |
| `xrootdForPnfs` | `kTRUE` | Stream resolved `/pnfs` paths over XRootD |
| `simTreeName` | `"Events"` | Name of the simulation tree |
| `recoTreeName` | `"Reco"` | Name of the reconstruction tree |
| `geoTreeName` | `"Geometry"` | Name of the geometry tree |
| `simBranchName` | `"Event"` | Branch of the simulation tree holding the event |
| `genieTreeName` | `"gst"` | Name of the GENIE tree |
| `outputTreeName` | `"AnaOutput"` | Name of the output tree |
| `requireSim` | `kTRUE` | Stop if the input has no simulation tree |
| `printProducts` | `kTRUE` | List the reconstruction products at start-up |
| `firstEvent` | `0` | First entry to process |
| `maxEvents` | `-1` | Number of entries to process; `-1` means all |
| `nReports` | `10` | Number of progress messages over the whole job |
| `printGeometry` | `kTRUE` | Print the detector configuration at start-up |
| `params` | empty | The analysis's own parameters, everything the job macro set that is not one of the fields above |

For quick jobs there is a shorthand:

```cpp
MyAnalysis().Execute("tpc_reco.root", "out.root");
```

### Input files

`inputFiles` and `genieFiles` are input *specifications*, not single paths. Each is a comma-separated list whose entries may be plain paths or shell wildcard patterns; `~` and `$VARIABLES` are expanded first. The resulting files are chained, so one job can span a whole production.

```bash
GArAnalysis -a MyAnalysis.C -i run1.root                    # one file
GArAnalysis -a MyAnalysis.C -i 'gun_*.root'                 # wildcard
GArAnalysis -a MyAnalysis.C -i 'gun_a.root, gun_b.root'     # list
GArAnalysis -a MyAnalysis.C -i 'set1/*.root, set2/*.root'   # several patterns
```

Quote a pattern so that the shell hands it over whole; an unquoted one the shell expands itself becomes several arguments, and `GArAnalysis` takes only one. The same specifications work as `/ana/global/input` in a job macro.

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

Set `/ana/global/xrootdForPnfs false` to read through the mount instead. This needs a ROOT built with XRootD support, which the CVMFS DUNE builds have.

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

Each analysis ships with a `.mac` of the same name, listing its parameters with their defaults. Run one with `-m` to use it, or without to run at the defaults throughout.

### `ExampleAnalysis.C`

Worked example, and the starting point for a new analysis. Clusters the ECal barrel hits of each layer with DBSCAN and writes one output entry per cluster.

```bash
GArAnalysis -a ExampleAnalysis.C -i tpc_reco.root -o example_out.root -m ExampleAnalysis.mac
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `/ana/eps` | `0.5` | Cluster radius handed to DBSCAN, cm |
| `/ana/minPts` | `5` | Neighbours a hit needs to be a core hit |

### `TrackingEfficiency.C`

Tracking efficiency of primary particles in the gas TPC, for particle-gun samples. The input has to have been through the TPC reconstruction, since it reads `TPCHits` and `TPCClusters`:

```bash
GArReconstruction -i gun.root -m tpc_reco.mac -o tpc_reco.root
GArAnalysis -a TrackingEfficiency.C -i tpc_reco.root -o tracking_eff.root \
            -m TrackingEfficiency.mac
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `/ana/minHits` | `20` | Hits the candidate needs to count as a findable track |
| `/ana/minLength` | `10.` | 3D span the candidate has to cover, cm |
| `/ana/minPurity` | `0.5` | Fraction of the candidate's hits that must belong to the primary |
| `/ana/minTrueLength` | `5.` | True path in the gas needed to enter the denominator, cm |
| `/ana/minHitFraction` | `0.5` | Share of a hit's true energy the primary needs to own it |
| `/ana/minClusterPurity` | `0.5` | Hit purity a cluster needs to join the candidate |
| `/ana/minClusterHits` | `3` | Matched hits a cluster needs to join the candidate |
| `/ana/excludeSecondaries` | `true` | Leave delta-ray deposits out of the true path length |

There is no pattern recognition in the chain yet, so what this measures is the question that comes before it: does the reconstruction leave behind enough of a primary for a track to be found at all?

Each primary (`motherID == 0`) is followed through in four steps.

**Denominator.** Primaries whose true path through the gas — `Sum(stepSize)` over their own TPC deposits — reaches `minTrueLength`. Neutrals and particles that barely clipped the gas drop out here on their own, so what remains measures the reconstruction rather than the acceptance.

**Matching.** A reconstructed hit belongs to the primary when the primary is among its `trackIDs` with a `trackFraction` of at least `minHitFraction`, i.e. when the primary's ionisation dominated the pulse. Delta rays count towards the primary here — the digitisation books their charge under the parent that made them, and no detector could separate the two — while the true path length above deliberately excludes them, so that it is the length of the particle's own track.

**Candidate.** Every cluster the primary dominates, meaning at least `minClusterHits` matched hits and a hit purity of at least `minClusterPurity`, joins the candidate. Several usually do: with the default clustering a found track is spread over four or five clusters, so restricting the candidate to the best one would badly understate what is there. The number of clusters it took is written out and reported at the end of the job, which is how fragmentation stays visible rather than being hidden by the merge.

**Numerator.** The candidate is a findable track when it holds at least `minHits` hits, spans at least `minLength` in 3D and is at least `minPurity` pure.

The 3D span is the distance between the two most widely separated hits of the candidate. For a straight track that is its length; for one that curls up in the field it is the chord across the spiral, which is the right thing to cut on, since a tight curler covers little ground however many hits it leaves. The extent along the candidate's principal axis is written out next to it as `candPcaLength`.

Two efficiencies are printed. The first is over the whole denominator and folds in the primaries whose true path is shorter than the `minLength` their candidate would have to span — those cannot be found however well they were reconstructed. The second is over just the ones long enough to be in reach, and is the one to quote for the reconstruction itself.

The output tree carries one entry per primary in the denominator, holding the raw quantities behind every cut — `nMatchedHits`, `nCandClusters`, `nCandHits`, `candLength`, `candPurity`, `candCompleteness`, the leading cluster's numbers on their own as `leadHits`/`leadLength`/`leadPurity`, and the truth (`momentum`, `momentumT`, `cosTheta`, `trueLength`, `trueEdep`) — so the working point can be moved afterwards without running the job again.

Four `TEfficiency` objects are written, with Clopper-Pearson intervals: `effP`, `effPt`, `effCosTheta` and `effLength`. The magnetic field runs along z, so the helix radius is set by the transverse momentum, which is why `effPt` is there beside `effP`. Their axes are binned over the range the sample turns out to cover, so the plots work whatever momenta the gun was thrown at. The same four are drawn to `tracking_efficiency.png`.

Note that on a single-particle gun sample almost every hit belongs to the one primary, so the purity cuts do nothing and `candPurity` sits at 1. They start to matter on samples where several particles cross the same part of the detector.

### `TruncatedDEDX.C`

Truncated mean dE/dx of primary particles in the gas TPC versus their initial momentum, for particle-gun samples.

```bash
GArAnalysis -a TruncatedDEDX.C -i gun.root -o dedx.root -m TruncatedDEDX.mac
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `/ana/truncation` | `0.6` | Fraction of the dE/dx samples to keep (lowest first) |
| `/ana/excludeSecondaries` | `true` | Drop deposits booked as a particle's secondaries |
| `/ana/requireExit` | `false` | Keep only tracks that leave the TPC (punch-through) |
| `/ana/minHits` | `10` | Minimum samples required per particle |
| `/ana/pMin` | `-1` | Lower edge of the momentum axis, GeV/c; negative takes it from the data |
| `/ana/pMax` | `-1` | Upper edge of the momentum axis, GeV/c; negative takes it from the data |
| `/ana/dedxMin` | `-1` | Lower edge of the dE/dx axis, keV/cm; negative uses 0 |
| `/ana/dedxMax` | `-1` | Upper edge of the dE/dx axis, keV/cm; negative takes it from the data |

The ranges apply to the binning of the plots, not to the event selection: entries outside them land in the under/overflow bins and are excluded from the profile, but the output tree always holds every primary that passed `minHits`.

Primaries are identified as the particles with `motherID == 0`. Each of their TPC hits gives one dE/dx sample, `energyDeposit / stepSize`; the lowest `truncation` of those samples are averaged. `SimEvent::TPCHitIsSecondary()` marks the deposits made by unstored secondaries — mostly delta rays — that the simulation folds back onto the parent (its `sec_tpcHits`), so excluding them is what isolates the parent's own ionisation.

`requireExit` keeps only the tracks that leave the TPC gas volume, so that the sample is restricted to particles whose range exceeds the detector and excludes those that stop inside and deposit their whole remaining energy there. A track counts as exiting when it stops outside the gas cylinder — `r > gar_tpc_radius` or `|z| > gar_tpc_length / 2`, taken from the geometry tree and applied to the track's end point, which the simulation records wherever the particle finally stops in the detector.

The macro converts to the units conventionally used for dE/dx plots — **GeV/c** and **keV/cm** — as the values are computed, so the output tree, the plots and the `pMin`/`pMax`/`dedxMin`/`dedxMax` parameters are all in those units.

The output tree carries `momentum` (GeV/c), `dedx` (keV/cm), `pdgCode`, `nHits`, `nHitsKept` and `trackLength` (cm) per primary. The macro also writes a 2D histogram of dE/dx against momentum with logarithmic momentum bins, overlays the profile of the mean, and saves `truncated_dedx_vs_momentum.png`. Left to themselves both axes are ranged from the data, so the plot works whatever the gas pressure and gun momentum range happen to be.

The output file also contains a `dedx_slices/` folder holding the 1D dE/dx distribution of each momentum bin, as `hDEDX_pbin01` … `hDEDX_pbin60`. These are Y projections of the 2D histogram, so slice `N` is exactly bin `N` of the plot; each title carries the momentum range of its bin. One histogram is written per bin, including the empty ones at the edges of the range, so the numbering always lines up with the 2D histogram.

```cpp
TFile f("dedx.root");
f.Get<TH1D>("dedx_slices/hDEDX_pbin20")->Draw();
```
