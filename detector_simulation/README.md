# GArSimulation

A simple Geant4 simulation featuring a high-pressure gaseous argon TPC with an electromagnetic calorimeter and muon tagger.

## Overview

This simulation models a detector with the following components:
- Cylindrical high-pressure gaseous argon Time Projection Chamber (TPC)
- Electromagnetic calorimeter (ECal) barrel with alternating layers of absorber and scintillator
- Circular ECal endcaps with the same layer configuration
- Muon tagger(MuID) barrel with layers of absorber and scintillator

## Features

- Fully configurable detector geometry through macro commands
- Customizable physics models (EM and hadronic)
- Particle gun with configurable properties
- Interface with GENIE and NuWro
- Energy deposition recording in all detector components
- Particle trajectory tracking

## Requirements

- Geant4 (10.7 or later recommended)
- CMake (3.16 or later)
- C++ compiler with C++17 support

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
./GArSimulation -m macros/run.mac
```

### Interactive Mode

Run the simulation with visualization:

```bash
./GArSimulation -v
```

## Configuring the Simulation

### Detector Geometry

You can modify the detector geometry using the following commands in a macro file:
```
/detector/TPCRadius 250 cm
/detector/TPCLength 500 cm
/detector/ECalLayers 42
/detector/MuIDLayers 3
```

### Particle Gun

Configure the particle source:
```
/gun/particleType mu-
/gun/energy 5 GeV
/gun/positionZ -50 cm
/gun/momentumDirection 0 0 1
```

### Physics Models

Select physics models:
```
/physics/emModel emstandard_opt4
/physics/hadronicModel FTFP_BERT
/physics/cutValue 1.0 mm
```
