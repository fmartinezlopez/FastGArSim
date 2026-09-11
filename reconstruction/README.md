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
- **TPCDigiModule**: TPC readout simulation — converts energy deposits to ionization electrons, drifts them to the pad plane, collects them with a parametrised pad response, and digitizes the resulting pad waveforms
- **TPCHitFinderModule**: TPC hit finding — finds charge pulses on the digitized waveforms and turns each into a 3D hit with a calibrated charge
- **TPCClusterModule**: TPC clustering — groups TPC hits with DBSCAN and summarizes each group's charge, centroid and principal axis

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

The output is a copy of the input with a `Reco` tree added to it, so `<output>` holds `Events`, `Geometry` and `Reco` together — see [Output Structure](#output-structure). Leave `-o` out and it is named after the input, `sim.root` becoming `sim_reco.root`.

### Where the macros are looked up

The configuration macros are copied next to the executable at build time, and `GArReconstruction` searches for the `-m` argument in

1. the working directory, so a local `macros/` still wins;
2. `FASTGARSIM_MACRO_PATH`, exported by the generated `setup.sh`;
3. the directory holding the executable, and its `macros/` subdirectory.

So `-m macros/ecal_digi.mac` and `-m ecal_digi.mac` both work from any directory, including from a relocated or installed build. If nothing matches, the error lists the directories that were searched. The same resolution is used by `GArSimulation`; it lives in [common/include/MacroPath.hh](../common/include/MacroPath.hh).

```bash
# all equivalent, from anywhere
GArReconstruction -i sim.root -m macros/ecal_digi.mac -o out.root
GArReconstruction -i sim.root -m ecal_digi.mac        -o out.root
```

### Command Line Options

- `-i <file>` : Input ROOT file from detector simulation (required)
- `-m <file>` : Macro file for reconstruction configuration (required)
- `-o <file>` : Output ROOT file (default: the input name with `_reco` before `.root`)
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

### TPCDigiModule Parameters

The pad plane radius and the drift length come from the `Geometry` tree; everything else is set here.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `padPitch` | double | 0.6 | Square pad pitch [cm] |
| `doubleSided` | bool | true | Cathode at `z = 0` and an anode at each end; `false` drifts the full length to one anode |
| `readoutSide` | int | 1 | Single-sided only: `+1` reads out at the `+z` end, `-1` at `-z` |
| `padResponseWidth` | double | 2.5 | Pad response `sigma` is `padPitch / padResponseWidth` |
| `padResponseShape` | double | 4 | Exponent `p` of the generalized Gaussian |
| `padResponseRange` | int | 1 | Pads either side of the arrival point that share the charge |
| `padResponseNormalize` | bool | true | Normalize the response over those pads; `false` loses the charge landing between pads |
| `wValue` | double | 26.4 | Energy per electron-ion pair [eV] |
| `fanoFactor` | double | 0.16 | Fano factor of the ionization statistics |
| `collectionEfficiency` | double | 0.7 | Fraction of electrons surviving to the pads |
| `driftVelocity` | double | 3.011 | Drift velocity [cm/us] |
| `electronLifetime` | double | 3e6 | Electron lifetime [us] |
| `diffusionT` | double | 0.0160 | Transverse diffusion [cm/sqrt(cm)] |
| `diffusionL` | double | 0.0201 | Longitudinal diffusion [cm/sqrt(cm)] |
| `electronsPerCluster` | double | 20 | Electrons drifted as one macro-electron |
| `minClusters` | int | 10 | Minimum macro-electrons per energy deposit |
| `samplingRate` | double | 20.0 | ADC sampling rate [MHz] |
| `shapingTime` | double | 0.1 | CR-RC² shaping time [us] |
| `gain` | double | 0.5 | Peak ADC counts per collected electron |
| `noiseLevel` | double | 5.0 | Electronics noise [electrons] |
| `adcPedestal` | double | 100.0 | Baseline [ADC] |
| `adcRange` | double | 4096 | ADC dynamic range |
| `adcThreshold` | double | 12.0 | Zero suppression threshold above pedestal [ADC] |
| `roiPadding` | int | 5 | Ticks kept either side of a threshold crossing |
| `maxTicks` | int | 0 | Readout window length; 0 derives it from the drift length |
| `saveTruth` | bool | true | Carry MC truth on the waveforms |
| `writeWaveforms` | bool | false | Also write the waveforms to the output tree |
| `seed` | int | 0 | RNG seed; 0 seeds from the system |

### TPCHitFinderModule Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold` | double | 12.0 | Pulse finding threshold above pedestal [ADC] |
| `minPulseTicks` | int | 2 | Shortest pulse accepted |
| `minPulseADC` | double | 0.0 | Smallest pulse integral accepted [ADC] |
| `valleyFraction` | double | 0.7 | A valley this far below the lower of two peaks splits them |
| `chargeCalibration` | double | 1.0 | Multiplies the reconstructed charge |
| `lifetimeCorrection` | bool | true | Undo the attenuation over the reconstructed drift time |

### TPCClusterModule Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `epsilon` | double | 2.0 | DBSCAN neighbourhood radius [cm] |
| `minPoints` | int | 3 | Neighbours needed for a hit to be a core hit |
| `minHits` | int | 3 | Smallest cluster kept |
| `separatePlanes` | bool | true | Never merge hits from the two drift volumes |

## TPC response model

`TPCDigiModule` applies the following chain to the simulated energy deposits, one at a time:

1. **Ionization.** The deposit gives `E / W` electron-ion pairs on average, fluctuating with the Fano factor (Poisson below 20 pairs, where the Gaussian is a poor description).
2. **Drift.** The TPC is a cylinder about the `z` axis, so the readout planes are its end faces at `z = ±L/2` and the drift is along `z`. With `doubleSided` the cathode sits at `z = 0` and each deposit drifts to the nearer end. The charge is attenuated by `exp(-t/lifetime)` and by `collectionEfficiency`; diffusion spreads it by `sigma = D sqrt(d)` transversally and longitudinally.
3. **Macro-electrons.** Rather than drifting single electrons, the deposit's charge is split into `electronsPerCluster`-sized groups (at least `minClusters` of them) and each group is diffused as a unit — the sampling that makes the model fast enough to run over full events.
4. **Pad response.** Each group lands somewhere on the plane and its charge is shared over the pads within `padResponseRange` according to
   `PRF(dx, dy) = exp(-(|dx|^p + |dy|^p) / 2 sigma^p)`, `sigma = padPitch / padResponseWidth`.
   Normalized (the default) the response only decides how the charge is *shared*, so the collection efficiency is entirely in `collectionEfficiency`.
5. **Electronics.** Each pad's charge is convolved with a CR-RC² shaping function normalized to unit peak, scaled by `gain`, given a pedestal and Gaussian noise, then clipped to the ADC range and quantized.
6. **Zero suppression.** Only the regions where the waveform crosses `adcThreshold`, padded by `roiPadding` ticks, are kept. A channel hit at two well separated drift times therefore yields two `TPCWaveform` objects.

`TPCHitFinderModule` inverts the last steps: the pedestal-subtracted integral divided by the ADC-per-electron the digitizer reports gives the collected charge, and the charge-weighted pulse time, corrected for the delay the shaping introduces, gives the drift time and so `z`. Both modules agree on these constants because the digitizer publishes them in the `RecoStore` as `TPCConditions`.

## Output Structure

By default the output file starts life as a copy of the input, with a `Reco` tree added to it. The result therefore holds the simulation's `Events` and `Geometry` trees alongside the reconstruction's own products, and an analysis needs one file rather than two that have to be kept in step by hand. Put

```
/reco/global/copyInput false
```

in the macro to write the products on their own instead; that saves the disk the copy costs, at the price of having to carry the simulation file along with the result. Re-running the reconstruction on a file that already has a `Reco` tree replaces that tree — the earlier pass's products are gone, though both passes stay in the file's history. Two module sets therefore have to be listed in one macro to end up in one file; a second pass does not add to the first.

`Reco` has one entry per event, aligned entry by entry with `Events`.

Two more trees describe the file:

- `Schema` — one entry per branch of every tree in the file: its C++ type, the module instance that produced it, and the parameters that module was given
- `Provenance` — one entry per reconstruction pass: when it ran, on what input, and the configuration macro verbatim

They exist because the module set is not fixed: which products a file holds depends on how the reconstruction was configured, and the file is the only place that can say. `DumpSchema output.root` prints them, and the analysis framework reads them to report a missing product at start-up rather than mid-loop. They are plain `std::string` branches, so uproot reads them without a dictionary.

Products are attributed to modules automatically — the manager watches the output tree either side of each module's `Initialize()` — so a new module needs no extra code to appear in the record.

### Global settings

| Command | Default | Effect |
| --- | --- | --- |
| `/reco/global/inputTreeName` | `Events` | Tree to read from the input file |
| `/reco/global/outputTreeName` | `Reco` | Tree the products are written to |
| `/reco/global/copyInput` | `true` | Start the output as a copy of the input |

### Product branches

#### ECalDigiModule

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

#### TPCDigiModule

- **`TPCWaveforms`** (`std::vector<digi::TPCWaveform>`) — one zero-suppressed waveform per region of interest, written only when `writeWaveforms` is set:
  - `channel`, `plane`, `row`, `col` — readout channel address (`plane` 0 = `+z` end, 1 = `-z` end)
  - `padX`, `padY` — pad centre [cm]
  - `tickStart`, `pedestal` — tick of the first sample and the baseline it sits on
  - `adc` — ADC samples, pedestal included
  - `trueEnergy` — true deposited energy per sample [MeV]
  - `trackIDs`, `trackFractions` — contributing track IDs and their energy fractions

#### TPCHitFinderModule

- **`TPCHits`** (`std::vector<digi::TPCHit>`) — one hit per charge pulse:
  - `channel`, `plane`, `row`, `col` — readout channel address
  - `x`, `y` — pad centre [cm]; `z` — drift coordinate from the pulse time [cm]
  - `driftTime` — charge-weighted drift time [us]
  - `peakAmplitude`, `adcSum` — pulse height and pedestal-subtracted integral [ADC]
  - `charge` — collected charge [electrons]; `energy` — calibrated energy [MeV]
  - `sigmaT`, `sigmaZ` — pulse width in time [us] and along the drift [cm], with the shaping width taken out
  - `startTick`, `endTick` — tick range of the pulse
  - `trueEnergy`, `trackIDs`, `trackFractions` — MC truth

#### TPCClusterModule

- **`TPCClusters`** (`std::vector<digi::TPCCluster>`) — one entry per group of hits:
  - `clusterID`, `plane`, `nHits` — `plane` is -1 if the hits span both drift volumes
  - `x`, `y`, `z` — charge-weighted centroid [cm]
  - `charge`, `energy` — summed hit charge [electrons] and energy [MeV]
  - `rmsX`, `rmsY`, `rmsZ` — charge-weighted spread about the centroid [cm]
  - `dirX`, `dirY`, `dirZ` — principal axis; `length`, `width` — extent along it and transverse to it [cm]
  - `startX/Y/Z`, `endX/Y/Z` — the hits at either end of the principal axis [cm]
  - `hitIndices` — indices into the `TPCHits` collection
  - `trueEnergy`, `trackIDs`, `trackFractions` — MC truth, ordered by contribution

Data type definitions are in [common/include/DigiDataTypes.hh](../common/include/DigiDataTypes.hh).

## Analysis

The compiled analysis framework in [analysis/](../analysis/) reads reconstruction products directly, asking for them by name. That is the route for anything beyond a quick look.

For a quick look, these example macros plot the digitization output. They need `libDigiDataDict` and the `common/` headers, which is what sourcing the build environment provides:

```bash
ECalDigiAnalysis reco_output.root
TPCRecoAnalysis tpc_reco.root
```

To see what a file holds before writing anything:

```bash
DumpSchema tpc_reco.root            # products, types and producing modules
DumpSchema tpc_reco.root true       # and the macros that produced them
```

To read the products from python, flatten them first:

```bash
MakeNtuple tpc_reco.root ntuple.root
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
│   ├── ECalDigiModule.hh       # ECal digitization module
│   ├── TPCDigiModule.hh        # TPC readout simulation
│   ├── TPCHitFinderModule.hh   # TPC hit finding
│   ├── TPCClusterModule.hh     # TPC clustering
│   ├── TPCReadoutGeometry.hh   # TPC pad plane and pad response
│   └── TPCConditions.hh        # TPC drift and electronics constants
├── src/                        # Source files
│   ├── RecoManager.cc
│   ├── RecoModule.cc
│   ├── ModuleFactory.cc
│   ├── MacroParser.cc
│   ├── RecoStore.cc
│   ├── ECalDigiModule.cc
│   ├── TPCDigiModule.cc
│   ├── TPCHitFinderModule.cc
│   ├── TPCClusterModule.cc
│   └── TPCReadoutGeometry.cc
├── macros/
│   ├── ecal_digi.mac           # ECalDigiModule configuration
│   └── tpc_reco.mac            # TPC digitization, hit finding and clustering
└── utils/
    ├── ECalDigiAnalysis.C      # Example macro for ECal digi output
    ├── TPCRecoAnalysis.C       # Example macro for TPC reco output
    └── RecoExample.C           # Generic reconstruction output example
```

## Architecture

### Core Components

1. **RecoManager** — handles file I/O and module lifecycle (`Initialize → Execute → Finalize`)
2. **RecoModule** — base class providing a common interface, parameter access, and handles to input/output trees and the RecoStore
3. **ModuleFactory** — singleton registry; modules self-register at program startup via a static initializer
4. **RecoStore** — key-value store for passing objects between modules within an event
5. **MacroParser** — reads the configuration macro and instantiates modules
6. **ProductSchema** ([common/](../common/include/ProductSchema.hh)) — the record of what the output file holds and how it was made, written by the manager and read by the analysis

### Data Flow

```
Input ROOT File → copied to the output → RecoManager → Module::Initialize()
                  (Events, Geometry)                    ↓        (branches booked
                                        For each event: Module::Execute()   here are
                                                        → Fill the Reco tree  attributed
                                                       ↓                      to the module)
                                                       Module::Finalize()
                                                       ↓
                                        Schema + Provenance written
                                                       ↓
                     Output ROOT File: Events, Geometry, Reco, Schema, Provenance
```

## Adding a New Reconstruction Module

1. **Create header** in [include/](include/) — inherit from `RecoModule`, override `GetInputSpec`, `GetOutputSpec`, `Initialize`, `Execute`, `Finalize`
2. **Create implementation** in [src/](src/) — register with `ModuleFactory` via a static initializer (see `ECalDigiModule.cc` for the pattern)
3. CMakeLists.txt picks up new `.cc` files automatically via `file(GLOB sources src/*.cc)` — no changes needed
4. Add the module to your macro with `/reco/addModule`
