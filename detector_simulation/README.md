# GArSimulation

A Geant4 simulation featuring a high-pressure gaseous argon TPC with an electromagnetic calorimeter and muon tagger for neutrino detector studies.

## Overview

This simulation models a detector with the following components:
- Cylindrical high-pressure gaseous argon Time Projection Chamber (TPC)
- Electromagnetic calorimeter (ECal) with separate barrel and endcap regions
  - Barrel and endcap ECal with configurable high-gain (HG) and low-gain (LG) layers
  - Alternating layers of lead absorber and plastic scintillator
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

## Requirements

- Geant4 (10.7 or later recommended)
- CMake (3.16 or later)
- C++ compiler with C++17 support
- ROOT

## Building the Project

```bash
# Create a build directory
mkdir build
cd build

# Configure with CMake
cmake ..

# Build the application
make -j4
```

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
/analysis/EnergyCut 0.001 MeV  # Energy threshold for recording hits
```

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

### Ntuple Maker

The object-based `Events` TTree is not directly suited for event-loop analyses. The macro [utils/EventToNtupleConverter.C](utils/EventToNtupleConverter.C) converts it into a flat ntuple that is easier to work with.

**Usage (from ROOT):**

```cpp
// Load the simulation dictionary first
gSystem->Load("libROOTDataDict");
.x utils/EventToNtupleConverter.C("input.root", "output_ntuple.root")
```

The converter produces a file with two TTrees:

- **`AnaTree`** — one entry per event, with all data stored as flat `std::vector` branches (indexed by particle or hit):
  - Particle identity: `eventID`, `trackID`, `pdgCode`, `motherID`, `creatorProcess`, `endProcess`
  - Trajectory endpoints: `startX/Y/Z`, `endX/Y/Z` [cm] and `startPX/Y/Z`, `endPX/Y/Z` [MeV/c]
  - TPC hits: `tpcHitTrackID`, `tpcHitIsSec`, `tpcHitX/Y/Z`, `tpcHitEdep`, `tpcHitStepSize`
  - ECal hits: `ecalHitTrackID`, `ecalHitIsSec`, `ecalHitX/Y/Z`, `ecalHitTime`, `ecalHitEdep`, `ecalHitSegment`, `ecalHitLayer`, `ecalHitDetID`
  - MuID hits: `muidHitTrackID`, `muidHitIsSec`, `muidHitX/Y/Z`, `muidHitTime`, `muidHitEdep`, `muidHitSegment`, `muidHitLayer`, `muidHitDetID`

  The `*IsSec` flag distinguishes direct hits from hits accumulated from unrecorded secondaries. All hit vectors include both primary and secondary contributions tagged accordingly.

- **`GeoTree`** — a copy of the `Geometry` tree from the simulation file, renamed for consistency.

This flat ntuple is the expected input format for the analysis macros in [analysis/](../analysis/).

## Authors

- Francisco Martinez Lopez — [frmart@iu.edu](mailto:frmart@iu.edu)
- Jude Martin — [j.martin24@imperial.ac.uk](mailto:j.martin24@imperial.ac.uk)
