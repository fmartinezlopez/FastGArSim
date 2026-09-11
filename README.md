# FastGArSim

A Geant4-based simulation, reconstruction and analysis framework for a high-pressure gaseous argon detector, targeting neutrino interaction studies.

The detector is a cylindrical high-pressure gaseous argon TPC in a magnetic field, surrounded by an electromagnetic calorimeter with separate barrel and endcap regions and a muon identification system. Geometry, physics models and event generation are all configurable at run time through macro commands.

## Modules

| Stage | Package | Produces |
| --- | --- | --- |
| 1. Detector simulation | [detector_simulation/](detector_simulation/) | `Events` and `Geometry` trees |
| 2. Reconstruction | [reconstruction/](reconstruction/) | a copy of its input with a `Reco` tree added, so the result holds `Events`, `Geometry` and `Reco` together |
| 3. Analysis | [analysis/](analysis/) | Whatever your macro writes |

The classes written to and read from these files are defined once, in [common/include/](common/include/), so that all three packages agree on them. The analysis reads them as they were written; there is no flat ntuple stage in between. [common/utils/MakeNtuple.C](common/utils/MakeNtuple.C) still produces one for reading outside the framework, from uproot or bare ROOT.

### What a file says about itself

The reconstruction is modular, so which products a file holds depends on which modules were run. Rather than leaving that to be remembered, every file it writes carries two extra trees describing itself:

- `Schema` — one entry per branch: the tree it lives in, its C++ type, and the module that produced it with the parameters it was given
- `Provenance` — one entry per processing pass: when it ran, on what input, and the configuration macro verbatim

`DumpSchema file.root` prints both. Analyses use the same record to fail at start-up, naming what the file does hold, rather than part-way through an event loop. Simulation files and files written before this existed have no `Schema` tree; readers fall back to walking the trees, so nothing has to be regenerated.

## Repository structure

```
FastGArSim/
├── CMakeLists.txt              # Top-level build, selects the sub-packages
├── setup_fnal.sh               # Environment setup for the FNAL machines
├── cmake/                      # Shared CMake settings, helpers and templates
├── common/                     # Shared data types and file-level utilities
│   ├── include/
│   │   ├── SimDataTypes.hh     # Simulation output types
│   │   ├── DigiDataTypes.hh    # Digitization types
│   │   ├── ProductSchema.hh    # What a file holds, and how it was made
│   │   └── TreeFlattener.hh    # Object trees -> flat vector ntuples
│   ├── src/                    # Their implementations (libGArCommon)
│   └── utils/                  # MakeNtuple, DumpSchema
├── detector_simulation/        # Geant4 simulation
│   ├── src/, include/          # Sources and headers
│   ├── macros/                 # Run and configuration macros
│   ├── utils/                  # Ntuple converter, geometry viewer
│   └── jobs/                   # Grid job submission tools
├── reconstruction/             # Modular reconstruction framework
│   ├── src/, include/          # Modules and framework
│   ├── macros/                 # Module configuration macros
│   └── utils/                  # Macros for reco output
└── analysis/                   # Compiled analysis framework
    ├── src/, include/          # Event readers, event loop, plotting, and utils
    └── macros/                 # Analyses built on the framework
```

## Requirements

- ROOT 6
- Geant4 11 (10.7 or later works), for the simulation only
- CMake >= 3.16
- A C++17 compiler

The reconstruction and analysis packages need only ROOT, so they build on a machine without Geant4.

## Building

Everything, in one go:

```bash
mkdir build && cd build
cmake ..
make -j4
```

### Selecting what to build

List the sub-packages you want. Names are `sim`, `reco`, `analysis` and `all`:

```bash
cmake -DCOMPONENTS="reco,analysis" ..   # skip the Geant4 simulation
cmake -DCOMPONENTS="sim" ..             # simulation only
```

Or switch them off one at a time:

```bash
cmake -DBUILD_SIMULATION=OFF ..
```

| Option | Default | Effect |
| --- | --- | --- |
| `COMPONENTS` | *(empty)* | Comma-separated list of sub-packages; overrides the `BUILD_*` options |
| `BUILD_SIMULATION` | `ON` | Build `GArSimulation` |
| `BUILD_RECONSTRUCTION` | `ON` | Build `GArReconstruction` |
| `BUILD_ANALYSIS` | `ON` | Build `GArAnalysis` and `libGArAnalysis` |
| `WITH_GEANT4_UIVIS` | `ON` | Build the simulation with the Geant4 UI and visualisation drivers |
| `CMAKE_BUILD_TYPE` | `RelWithDebInfo` | Standard CMake build types |

Each sub-package also still configures on its own, from its own directory — useful when working on one of them:

```bash
cd analysis && mkdir build && cd build && cmake .. && make -j4
```

`common/` is pulled in automatically either way, so a standalone build gets the same data type dictionaries as the full one.

### Build layout

The sub-packages keep their own subdirectory of the build tree, so the layout is the same whether you built everything or just one package:

```
build/
├── setup.sh                    # Environment for this build (generated)
├── rootlogon.C                 # Same, for ROOT sessions (generated)
├── common/                     # libSimDataDict, libDigiDataDict + .pcm/.rootmap,
│                               # libGArCommon, MakeNtuple, DumpSchema
├── detector_simulation/        # GArSimulation, macros/, GeoVis.C
├── reconstruction/             # GArReconstruction, libRecoDataDict, macros/
└── analysis/                   # GArAnalysis, libGArAnalysis, macros/
```

## Setting up the environment

CMake writes a `setup.sh` into the build directory describing what was built. Source it to use the build from any directory — it puts the executables on `PATH`, the shared libraries on the dynamic loader path and the headers on `ROOT_INCLUDE_PATH`:

```bash
source build/setup.sh
GArSimulation -m macros/gun.mac
```

It also writes a `rootlogon.C`, which ROOT runs automatically when started from the build directory, loading the dictionaries and `libGArAnalysis`. To get it in any directory, run `root -l "$FASTGARSIM_ROOTLOGON"`, or point `Rint.Logon` in your `~/.rootrc` at it.

### Command line tools

The utility macros are also built as executables, installed into `bin` and put on `PATH` by `setup.sh`. They take the macro's arguments positionally, so there is no ROOT invocation to quote and no environment to set up:

```bash
MakeNtuple reco.root ntuple.root
ECalDigiAnalysis reco.root
```

These tools are discovered, not listed. Each package scans one directory and builds an executable per macro it finds, named after the macro with the `.C` dropped:

| Directory scanned | Tools currently built |
| --- | --- |
| `common/utils/` | `MakeNtuple`, `DumpSchema` |
| `detector_simulation/utils/` | `GeoVis` |
| `reconstruction/utils/` | `ECalDigiAnalysis`, `RecoExample`, `TPCRecoAnalysis` |

Dropping a new macro into one of those directories is all it takes to get a tool for it — no CMake edits. The macro's entry function has to share its file name, which is the convention ROOT already requires for `.x Macro.C`.

Each tool takes `-h`/`--help`, which reports the argument count it accepts and the description taken from the macro's own header comment. Optional arguments default to whatever the macro declares, and anything that draws runs in batch mode and writes its canvases to file (set `FASTGARSIM_NO_BATCH` to override).

A macro that only makes sense interactively, or that cannot be compiled, is named in the `EXCLUDE` list of its package's `fastgarsim_add_macro_apps()` call.

The analyses in `analysis/macros/` work the other way round. They are not built at all: `GArAnalysis` takes one as an argument and compiles it when the job starts, the same way `GArSimulation` and `GArReconstruction` take a macro describing what to do:

```bash
GArAnalysis -a TruncatedDEDX.C -i 'gun_*.root' -o dedx.root -m TruncatedDEDX.mac
```

so writing an analysis needs no rebuild of anything, and a run always uses the macro as it stands at that moment. See the [analysis README](analysis/README.md#running-an-analysis).

### On the FNAL machines

[setup_fnal.sh](setup_fnal.sh) sets up the UPS products (CMake, Geant4, ROOT) from CVMFS and then the build on top of them. Source it, don't run it:

```bash
source setup_fnal.sh
```

If there is no build yet it says so and tells you how to make one; source it again afterwards. A build directory other than `<source>/build` can be given as an argument, and the product versions can be overridden beforehand:

```bash
export FASTGARSIM_ROOT_VERSION=v6_28_12
source setup_fnal.sh /path/to/build
```

Credentials for reading files from dCache are not set up — run `setup_fnal_security` afterwards when you need them. See the [analysis README](analysis/README.md#reading-from-dcache-over-xrootd) for what that involves.

## Quick start

The whole chain, from a build of everything:

```bash
source build/setup.sh

# 1. Simulate. macros/nu.mac reads GENIE events, macros/gun.mac fires a particle gun
GArSimulation -m macros/gun.mac

# 2. Reconstruct. The output is a copy of the input with a Reco tree added, so
#    tpc_reco.root holds Events, Geometry and Reco together
GArReconstruction -i output.root -m macros/tpc_reco.mac -o tpc_reco.root

# 3. See what came out
DumpSchema tpc_reco.root

# 4. Analyse. The analysis macro is compiled when the job starts; the simulation
#    objects and the reconstruction products are both read straight from the file
GArAnalysis -a ExampleAnalysis.C -i tpc_reco.root -o example_out.root

# 5. Optional: flatten it for uproot or bare ROOT
MakeNtuple tpc_reco.root ntuple.root
```

Step 2 is optional: an analysis that only uses truth-level information runs on `output.root` directly.

Run the simulation with `-v` instead of `-m` for the interactive viewer; that is what creates the Geant4 visualisation manager.

`GArSimulation`, `GArReconstruction` and `GArAnalysis` find their macros through `FASTGARSIM_MACRO_PATH` (exported by `setup.sh`) and through their own location, so they can be started from any directory. Macros still resolve against the working directory first, so a local `macros/` overrides the built-in one.

## Where to go next

| Document | Covers |
| --- | --- |
| [detector_simulation/README.md](detector_simulation/README.md) | Detector geometry, physics and generator commands, the output format, the ntuple branches |
| [detector_simulation/jobs/README.md](detector_simulation/jobs/README.md) | Grid production: tarballs, job submission, the configuration script |
| [reconstruction/README.md](reconstruction/README.md) | Reconstruction modules and their parameters, output branches, writing a new module |
| [analysis/README.md](analysis/README.md) | Writing an analysis, input specifications, reading from dCache, the shipped macros |

## Citation

See [CITATION.cff](CITATION.cff).

## Authors

- Francisco Martinez Lopez — [frmart@iu.edu](mailto:frmart@iu.edu)
- Jude Martin — [j.martin24@imperial.ac.uk](mailto:j.martin24@imperial.ac.uk)
