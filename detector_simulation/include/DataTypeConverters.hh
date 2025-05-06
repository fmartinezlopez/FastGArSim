#ifndef DATATYPECONVERTERS_HH
#define DATATYPECONVERTERS_HH

#include "DataTypes.hh"       // Original Geant4 data types
#include "ROOTDataTypes.hh"   // ROOT-compatible data types

// Forward declarations of conversion functions
root::TPCHit ConvertTPCHit(const TPCHit& g4Hit);
root::ECalHit ConvertECalHit(const ECalHit& g4Hit);
root::MuIDHit ConvertMuIDHit(const MuIDHit& g4Hit);
root::TrajectoryPoint ConvertTrajectoryPoint(const TrajectoryPoint& g4Point);
root::MomentumPoint ConvertMomentumPoint(const MomentumPoint& g4MomPoint);
root::Trajectory ConvertTrajectory(const Trajectory& g4Trajectory);
root::Particle ConvertParticle(const Particle& g4Particle);
root::Event ConvertEvent(const Event& g4Event);

#endif