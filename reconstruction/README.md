# GArReconstruction

Modular reconstruction framework for the FastGArSim detector simulation.

## Overview

This directory contains a modular reconstruction framework for processing simulation output from `detector_simulation`. Reconstruction algorithms are implemented as independent modules that can be configured and composed via macro files.

### Key Features

- **Modular architecture**: Reconstruction algorithms are independent, self-registering modules
- **Macro configuration**: Specify which modules to run and their parameters via text files
- **Extensible**: Add new modules without modifying the core framework
- **ROOT-based I/O**: No Geant4 dependency, only ROOT required

### Available Modules

- **ECalDigiModule**: ECal digitization — maps raw Geant4 ECal hits to readout channels (HG tiles and LG strips), applies SiPM photon statistics, waveform digitization, and time-based clustering

## Building

Built as part of the top-level FastGArSim build, or on its own — it needs ROOT but not
Geant4:

```bash
cd reconstruction
mkdir build && cd build
cmake .. && make -j4
```

Build options and environment setup are in the
[top-level README](../README.md#building).

## Usage

### Basic Reconstruction

```bash
GArReconstruction -i <input_simulation.root> -m <config.mac> -o <output_reconstruction.root>
```

### Command Line Options

- `-i <file>` : Input ROOT file from detector simulation (required)
- `-m <file>` : Macro file for reconstruction configuration (required)
- `-o <file>` : Output ROOT file (default: `reconstruction_output.root`)
- `-h, --help` : Show help message

## Configuration Macros

### Macro Syntax

```
# Set global options
/reco/global/inputTreeName Events

# Add a module
/reco/addModule <instanceName> <ModuleType>

# Set module parameters
/reco/<instanceName>/<parameter> <value>
```

### ECalDigiModule Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tileSize` | double | 0.5 | HG tile readout cell size [cm] |
| `stripSize` | double | 0.5 | LG strip width [cm] |
| `timeRes` | double | 0.1 | Time bin width / resolution [ns] |
| `propVelocity` | double | 10.0 | Light propagation velocity in strip [cm/ns] |
| `lightYield` | double | 1000.0 | Scintillator light yield [PE/MeV] |
| `effectivePE` | double | 5000.0 | SiPM effective number of pixels |
| `sipmNoise` | double | 1.0 | SiPM dark noise [PE RMS] |
| `sipmGain` | double | 10.0 | SiPM gain (PE → ADC counts) |
| `attLength` | double | 150.0 | Scintillator attenuation length [cm] |
| `adcRange` | int | 4096 | ADC dynamic range |
| `adcThreshold` | int | 10 | ADC threshold for cluster finding |
| `minClusterBins` | int | 2 | Minimum time bins per cluster |
| `maxPropTime` | double | 50.0 | Maximum propagation time for strip coincidence [ns] |

## Output Structure

The reconstruction creates a ROOT file with a `RecoTree` containing one entry per event.

### ECalDigiModule Output Branches

- **`ECalTileDigiHits`** (`std::vector<digi::TileDigiHit>`) — digitized hits in HG (fine-granularity tile) layers:
  - `segment`, `layer`, `row`, `col` — readout channel address
  - `adcSum` — integrated ADC counts
  - `time` — ADC-weighted hit time [ns]
  - `trueEnergy` — true deposited energy [MeV]
  - `trackIDs`, `trackFractions` — contributing track IDs and their energy fractions

- **`ECalStripDigiHits`** (`std::vector<digi::StripDigiHit>`) — digitized hits in LG (strip) layers:
  - `segment`, `layer`, `row` — readout channel address
  - `stripLength` — length of the strip [cm]
  - `adcLeft`, `adcRight` — ADC counts at each strip end
  - `adcCombined` — geometric mean of `adcLeft` and `adcRight`
  - `timeLeft`, `timeRight` — arrival times at each strip end [ns]
  - `recoTime` — reconstructed hit time (corrected for propagation) [ns]
  - `recoPosition` — reconstructed position along the strip from time difference [cm]
  - `trueEnergy` — true deposited energy [MeV]
  - `trackIDs`, `trackFractions` — contributing track IDs and their energy fractions

Data type definitions are in [common/include/DigiDataTypes.hh](../common/include/DigiDataTypes.hh).

## Analysis

Example ROOT macro for plotting ECalDigiModule output. It needs `libDigiDataDict` and
the `common/` headers, which is what sourcing the build environment provides:

```bash
ECalDigiAnalysis reco_output.root
```

## Directory Structure

```
reconstruction/
├── CMakeLists.txt              # Build configuration
├── GArReconstruction.cc        # Main executable
├── include/                    # Header files
│   ├── RecoDataTypes.hh        # Reconstruction output types
│   ├── RecoDataTypesLinkDef.h  # ROOT dictionary linkdef
│   ├── RecoManager.hh          # I/O and module orchestration
│   ├── RecoModule.hh           # Base class for reconstruction modules
│   ├── ModuleFactory.hh        # Self-registration factory
│   ├── MacroParser.hh          # Configuration macro parser
│   ├── RecoStore.hh            # Inter-module data store
│   └── ECalDigiModule.hh       # ECal digitization module
├── src/                        # Source files
│   ├── RecoManager.cc
│   ├── RecoModule.cc
│   ├── ModuleFactory.cc
│   ├── MacroParser.cc
│   ├── RecoStore.cc
│   └── ECalDigiModule.cc
├── macros/
│   └── ecal_digi.mac             # ECalDigiModule configuration
└── utils/
    ├── ECalDigiAnalysis.C      # Example macro for ECal digi output
    └── RecoExample.C           # Generic reconstruction output example
```

## Architecture

### Core Components

1. **RecoManager** — handles file I/O and module lifecycle (`Initialize → Execute → Finalize`)
2. **RecoModule** — base class providing a common interface, parameter access, and handles to input/output trees and the RecoStore
3. **ModuleFactory** — singleton registry; modules self-register at program startup via a static initializer
4. **RecoStore** — key-value store for passing objects between modules within an event
5. **MacroParser** — reads the configuration macro and instantiates modules

### Data Flow

```
Input ROOT File → RecoManager → Module::Initialize()
                                ↓
                 For each event: Module::Execute()
                                 → Fill output tree
                                ↓
                                Module::Finalize()
                                ↓
                             Output ROOT File
```

## Adding a New Reconstruction Module

1. **Create header** in [include/](include/) — inherit from `RecoModule`, override `GetInputSpec`, `GetOutputSpec`, `Initialize`, `Execute`, `Finalize`
2. **Create implementation** in [src/](src/) — register with `ModuleFactory` via a static initializer (see `ECalDigiModule.cc` for the pattern)
3. CMakeLists.txt picks up new `.cc` files automatically via `file(GLOB sources src/*.cc)` — no changes needed
4. Add the module to your macro with `/reco/addModule`
