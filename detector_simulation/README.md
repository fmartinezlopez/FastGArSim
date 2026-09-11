# GArSimulation

A Geant4 simulation featuring a high-pressure gaseous argon TPC with an electromagnetic calorimeter and muon tagger for neutrino detector studies.

## Overview

This simulation models a detector with the following components:
- Cylindrical high-pressure gaseous argon Time Projection Chamber (TPC)
- Electromagnetic calorimeter (ECal) with separate barrel and endcap regions
  - Barrel and endcap ECal with configurable high-gain (HG) and low-gain (LG) layers
  - Alternating layers of absorber and scintillator
- Muon identification system (MuID) barrel with layers of absorber and scintillator
- Configurable magnetic field

## Features

- Fully configurable detector geometry through macro commands
- Customizable physics models
- Multiple event generator options:
  - Particle gun with extensive configuration options
  - GENIE neutrino event generator interface
  - NuWro neutrino event generator interface
- Energy deposition recording in all detector components
- Particle trajectory tracking
- Multi-threaded event processing support
- Grid job submission tools for large-scale production

## Building

Built as part of the top-level FastGArSim build, or on its own:

```bash
mkdir build && cd build
cmake .. && make -j4
```

Requirements, build options and environment setup (including the FNAL machines) are in the [top-level README](../README.md#building).

## Running the Simulation

### Batch Mode

Run the simulation with a macro file:

```bash
# Using neutrino events from GENIE
./GArSimulation -m macros/nu.mac

# Using particle gun
./GArSimulation -m macros/gun.mac
```

### Interactive Mode

Run the simulation with visualization:

```bash
./GArSimulation -v
```

`-v` is what constructs the Geant4 visualization manager, and it also executes [macros/vis.mac](macros/vis.mac) for you. The `/vis/` commands do not exist without it, so `-m macros/vis.mac` on its own fails with `COMMAND NOT FOUND </vis/open ...>`.

### Where the macros are looked up

The run macros refer to each other by paths relative to the application directory (`/control/execute macros/init.mac`). `GArSimulation` therefore sets Geant4's macro search path at start-up to the working directory, then `FASTGARSIM_MACRO_PATH` (exported by the generated `setup.sh`), then the directory holding the executable and its `macros/` subdirectory. The program can consequently be started from anywhere, including from a relocated or tarballed build as used by the grid jobs. Because the working directory is searched first, a local `macros/` still takes precedence.

## Configuring the Simulation

### Generator Selection

Choose between different event generators using `/generator/select`:

**Neutrino Events (GENIE):**
```
/generator/select genie
/generator/genieFile /path/to/genie.gst.root
/generator/initialEvent 0
```

**Neutrino Events (NuWro):**
```
/generator/select nuwro
/generator/nuwroFile /path/to/nuwro.root
/generator/initialEvent 0
```

**Particle Gun:**
```
/generator/select particle
```

### Particle Gun Configuration

The particle gun offers extensive configuration options:

```
# Particle type and momentum
/particle/particleType pi0
/particle/momentum 2500.0 MeV
/particle/momentumSpread 2500.0 MeV
/particle/momentumDistribution uniform

# Position configuration
/particle/position 0 0 0 cm
/particle/positionSpread 260 260 250 cm
/particle/positionDistribution uniform
/particle/positionRMax 260 cm  # Maximum radial position

# Angular configuration
/particle/angleXZ 90.0 deg
/particle/angleXZSpread 0.0 deg
/particle/angleXZDistribution isotropic
/particle/angleXY 180.0 deg
/particle/angleXYSpread 180.0 deg
/particle/angleXYDistribution uniform
```

### Detector Geometry

Configure the detector geometry parameters:

```
# TPC configuration
/detector/TPCRadius 260 cm
/detector/TPCLength 500 cm
/detector/GasPressure 10.0 bar
/detector/TPCMaxStep 1.0 mm
/detector/BField 0.5 tesla

# ECal absorber and scintillator materials
/detector/ECalAbsorberMaterial G4_Pb
/detector/ECalScintillatorMaterial G4_PLASTIC_SC_VINYLTOLUENE

# ECal layer thicknesses
/detector/ECalHGAbsorberThickness 0.7 mm
/detector/ECalLGAbsorberThickness 1.4 mm
/detector/ECalHGScintillatorThickness 5 mm
/detector/ECalLGScintillatorThickness 10 mm

# ECal layer counts
/detector/ECalBarrelHGLayers 8
/detector/ECalBarrelLGLayers 34
/detector/ECalEndcapHGLayers 6
/detector/ECalEndcapLGLayers 36

# MuID configuration
/detector/MuIDLayers 3
```

### Physics Models

Select physics models and production cuts:

```
/physics/emModel emstandard_opt4
/physics/hadronicModel FTFP_BERT
/physics/cutValue 1.0 mm
```

### Output Configuration

```
/run/OutputFileName output_name
/analysis/TPCEnergyCut 0.0 MeV     # Energy threshold for recording TPC gas hits
/analysis/CaloEnergyCut 0.001 MeV  # Energy threshold for recording ECal and MuID hits
```

The two thresholds are deliberately different.

### Stepping in the TPC Gas

`/detector/TPCMaxStep` sets the maximum step length in the gas volume, and therefore the granularity of the energy deposits handed to the drift simulation. The gas is thin enough that Geant4's own step limits are of order metres for a GeV-scale track, so without an explicit limit a particle crosses the whole TPC in a handful of steps. The total energy
loss is still correct, but the deposits are far too sparse to seed drift.

A limit of 1--2 mm is a sensible default: pads are a few mm across and transverse diffusion over the full drift is itself of that order, so finer stepping costs CPU and output size without buying resolution.

### Multi-threading

```
/run/numberOfThreads 1  # Set number of threads (use 1 for single-threaded)
```

## Example Macros

The [macros](macros/) directory contains several example configurations:

- [init.mac](macros/init.mac) - Initialization macro for GENIE neutrino events
- [nu.mac](macros/nu.mac) - Run neutrino events from GENIE file
- [init_gun.mac](macros/init_gun.mac) - Initialization macro for particle gun
- [gun.mac](macros/gun.mac) - Run particle gun events
- [vis.mac](macros/vis.mac) - Visualization configuration

## Grid Job Submission

For large-scale production runs, see the [jobs](jobs/) directory for grid submission scripts and configuration tools. The automated configuration script simplifies setting up batch jobs on computing grids.

## Analysis

ROOT-based analysis macros are available in the `../analysis/` directory for processing simulation output files.

## Output Format

The simulation writes a ROOT file containing two TTrees:

- **`Events`** — one entry per simulated event, stored as a `root::Event` object with the following structure:
  - `eventID` — event index
  - `particles` — vector of `root::Particle` objects, one per tracked particle. Each particle contains:
    - `trackID`, `pdgCode`, `motherID` — particle identity and parentage
    - `creatorProcess`, `endProcess` — Geant4 process names at creation and termination
    - `trajectory` — step-by-step positions, momenta, and volume names
    - `tpcHits` / `ecalHits` / `muidHits` — direct energy deposits (position, energy, step length / layer / segment / detID)
    - `sec_tpcHits` / `sec_ecalHits` / `sec_muidHits` — energy deposits from secondary particles that were not individually recorded, accumulated onto the nearest recorded ancestor

  Only particles satisfying at least one of the following criteria are stored: primary particles, decay or conversion daughters of primaries, particles with TPC path length > 1 cm, or particles originating in the TPC and stopping in the ECal. Energy deposits from all other secondaries are folded into their nearest recorded ancestor.

- **`Geometry`** — one entry storing the detector geometry parameters used in the run (TPC radius/length/pressure/field, ECal and MuID layer counts and thicknesses, etc.).

The data types are defined in [common/include/SimDataTypes.hh](../common/include/SimDataTypes.hh) and shared between the simulation and the reconstruction.

### Flat ntuples

The `Events` tree holds objects, which the analysis framework reads directly — see the [analysis README](../analysis/README.md). For reading the output outside that framework, from uproot or a bare ROOT session, [common/utils/MakeNtuple.C](../common/utils/MakeNtuple.C) converts any FastGArSim file into flat `std::vector` branches:

```bash
MakeNtuple simulation.root ntuple.root
```

Nothing in it is specific to the simulation: the columns are worked out from the ROOT dictionaries of whatever the file contains, so the same tool flattens reconstruction output with any set of modules. It writes one flat tree per input tree, keeping the names, and a `Schema` tree recording which column came from which branch. See the [main README](../README.md) for the naming rule and the options.
