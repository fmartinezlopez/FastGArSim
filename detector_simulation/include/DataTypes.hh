#ifndef DATATYPES_HH
#define DATATYPES_HH

#include "G4ThreeVector.hh"
#include <vector>
#include <numeric>

struct TrajectoryPoint
{
  G4ThreeVector position;
  G4double time;

  // Copy constructor (automatically generated, but listed explicitly for clarity)
  TrajectoryPoint(const TrajectoryPoint& other) = default;
  
  // Assignment operator (automatically generated, but listed explicitly for clarity)
  TrajectoryPoint& operator=(const TrajectoryPoint& other) = default;
};

struct MomentumPoint
{
  G4ThreeVector momentum;

  // Copy constructor (automatically generated, but listed explicitly for clarity)
  MomentumPoint(const MomentumPoint& other) = default;
  
  // Assignment operator (automatically generated, but listed explicitly for clarity)
  MomentumPoint& operator=(const MomentumPoint& other) = default;
};

class Trajectory
{
public:
  Trajectory() = default;
  
  // Copy constructor
  Trajectory(const Trajectory& other)
    : points(other.points),
      mom_points(other.mom_points),
      start_vol_name(other.start_vol_name),
      stop_vol_name(other.stop_vol_name)
  {}
  
  // Assignment operator
  Trajectory& operator=(const Trajectory& other)
  {
    if (this != &other) {
      points = other.points;
      mom_points = other.mom_points;
      start_vol_name = other.start_vol_name;
      stop_vol_name = other.stop_vol_name;
    }
    return *this;
  }
  
  void AddPoint(const G4ThreeVector& pos, G4double t, G4String vol_name) {
    if (points.size() == 0) start_vol_name = vol_name;
    points.push_back({pos, t});
    stop_vol_name = vol_name;
  }

  void AddMomPoint(const G4ThreeVector& mom) {
    mom_points.push_back({mom});
  }

  const std::vector<TrajectoryPoint>& GetPoints() const { return points; }
  const std::vector<MomentumPoint>& GetMomPoints() const { return mom_points; }
  const G4String GetStartVolName() const { return start_vol_name; }
  const G4String GetStopVolName()  const { return stop_vol_name; }
  
private:
  std::vector<TrajectoryPoint> points;
  std::vector<MomentumPoint> mom_points;
  G4String start_vol_name;
  G4String stop_vol_name;
};

class TPCHit
{
public:
  TPCHit(const G4ThreeVector& pos, G4double edep, G4double step) 
    : position(pos), energyDeposit(edep), stepSize(step) {}
  
  // Copy constructor
  TPCHit(const TPCHit& other)
    : position(other.position),
      energyDeposit(other.energyDeposit),
      stepSize(other.stepSize)
  {}
  
  // Assignment operator
  TPCHit& operator=(const TPCHit& other)
  {
    if (this != &other) {
      position = other.position;
      energyDeposit = other.energyDeposit;
      stepSize = other.stepSize;
    }
    return *this;
  }
  
  const G4ThreeVector& GetPosition() const { return position; }
  G4double GetEnergyDeposit() const { return energyDeposit; }
  G4double GetStepSize() const { return stepSize; }
  
private:
  G4ThreeVector position;
  G4double energyDeposit;
  G4double stepSize;
};

class ECalHit
{
public:
  ECalHit(const G4ThreeVector& pos, G4double edep) 
    : position(pos), energyDeposit(edep) {}
  
  // Copy constructor
  ECalHit(const ECalHit& other)
    : position(other.position),
      energyDeposit(other.energyDeposit)
  {}
  
  // Assignment operator
  ECalHit& operator=(const ECalHit& other)
  {
    if (this != &other) {
      position = other.position;
      energyDeposit = other.energyDeposit;
    }
    return *this;
  }
  
  const G4ThreeVector& GetPosition() const { return position; }
  G4double GetEnergyDeposit() const { return energyDeposit; }
  
private:
  G4ThreeVector position;
  G4double energyDeposit;
};

class MuIDHit
{
public:
  MuIDHit(const G4ThreeVector& pos, G4double edep) 
    : position(pos), energyDeposit(edep) {}
  
  // Copy constructor
  MuIDHit(const MuIDHit& other)
    : position(other.position),
      energyDeposit(other.energyDeposit)
  {}
  
  // Assignment operator
  MuIDHit& operator=(const MuIDHit& other)
  {
    if (this != &other) {
      position = other.position;
      energyDeposit = other.energyDeposit;
    }
    return *this;
  }
  
  const G4ThreeVector& GetPosition() const { return position; }
  G4double GetEnergyDeposit() const { return energyDeposit; }
  
private:
  G4ThreeVector position;
  G4double energyDeposit;
};

class Particle
{
public:
  Particle(G4int id, G4int pdg, G4String process, G4int momid) : trackID(id), pdgCode(pdg), creatorProcess(process), motherID(momid) {}

  // Copy constructor
  Particle(const Particle& other)
    : trackID(other.trackID),
      pdgCode(other.pdgCode),
      creatorProcess(other.creatorProcess),
      endProcess(other.endProcess),
      motherID(other.motherID),
      trajectory(other.trajectory),
      tpcHits(other.tpcHits),
      ecalHits(other.ecalHits),
      muidHits(other.muidHits),
      sec_tpcHits(other.sec_tpcHits),
      sec_ecalHits(other.sec_ecalHits),
      sec_muidHits(other.sec_muidHits)
  {}
  
  // Assignment operator
  Particle& operator=(const Particle& other)
  {
    if (this != &other) {
      trackID = other.trackID;
      pdgCode = other.pdgCode;
      creatorProcess = other.creatorProcess;
      endProcess = other.endProcess;
      motherID = other.motherID;
      trajectory = other.trajectory;
      tpcHits = other.tpcHits;
      ecalHits = other.ecalHits;
      muidHits = other.muidHits;
      sec_tpcHits = other.sec_tpcHits;
      sec_ecalHits = other.sec_ecalHits;
      sec_muidHits = other.sec_muidHits;
    }
    return *this;
  }
  
  // Clone method for convenience
  Particle* Clone() const
  {
    return new Particle(*this);
  }

  void SetEndProcess(const G4String process) { endProcess = process; }

  void SetTrajectory(const Trajectory& traj) { trajectory = traj; }
  void AddTPCHit(const TPCHit& hit) { tpcHits.push_back(hit); }
  void AddECalHit(const ECalHit& hit) { ecalHits.push_back(hit); }
  void AddMuIDHit(const MuIDHit& hit) { muidHits.push_back(hit); }

  void AddSecondaryContribution(const Particle* part) {
    for (const auto& hit : part->GetTPCHits())  sec_tpcHits.push_back(hit);
    for (const auto& hit : part->GetECalHits()) sec_ecalHits.push_back(hit);
    for (const auto& hit : part->GetMuIDHits()) sec_muidHits.push_back(hit);
  }
  
  G4int GetTrackID() const { return trackID; }
  G4int GetPDGCode() const { return pdgCode; }
  G4String GetCreatorProcess() const { return creatorProcess; }
  G4String GetEndProcess() const { return endProcess; }
  G4int GetMotherID() const { return motherID; }

  const G4ThreeVector& GetStartPosition() const { return trajectory.GetPoints().front().position; }
  const G4ThreeVector& GetEndPosition() const { return trajectory.GetPoints().back().position; }
  const G4ThreeVector& GetMomentum() const { return trajectory.GetMomPoints().front().momentum; }

  const float GetTPCPathLength() const { return std::accumulate(tpcHits.begin(), tpcHits.end(), 0.0, [](float sum, const TPCHit& hit) { return sum + hit.GetStepSize(); }); }
  
  // Add const version of accessors
  Trajectory& GetTrajectory() { return trajectory; }
  const Trajectory& GetTrajectory() const { return trajectory; }
  const std::vector<TPCHit>& GetTPCHits() const { return tpcHits; }
  std::vector<TPCHit>& GetTPCHits() { return tpcHits; }
  const std::vector<ECalHit>& GetECalHits() const { return ecalHits; }
  std::vector<ECalHit>& GetECalHits() { return ecalHits; }
  const std::vector<MuIDHit>& GetMuIDHits() const { return muidHits; }
  std::vector<MuIDHit>& GetMuIDHits() { return muidHits; }

  //
  const std::vector<TPCHit>& GetSecondaryTPCHits() const { return sec_tpcHits; }
  std::vector<TPCHit>& GetSecondaryTPCHits() { return sec_tpcHits; }
  const std::vector<ECalHit>& GetSecondaryECalHits() const { return sec_ecalHits; }
  std::vector<ECalHit>& GetSecondaryECalHits() { return sec_ecalHits; }
  const std::vector<MuIDHit>& GetSecondaryMuIDHits() const { return sec_muidHits; }
  std::vector<MuIDHit>& GetSecondaryMuIDHits() { return sec_muidHits; }
  
private:
  G4int trackID;
  G4int pdgCode;
  G4String creatorProcess;
  G4String endProcess;
  G4int motherID;
  Trajectory trajectory;
  std::vector<TPCHit>  tpcHits;
  std::vector<ECalHit> ecalHits;
  std::vector<MuIDHit> muidHits;
  // Secondary collections (only for re-processed particles)
  std::vector<TPCHit>  sec_tpcHits;
  std::vector<ECalHit> sec_ecalHits;
  std::vector<MuIDHit> sec_muidHits;
};

class Event
{
public:
  Event(G4int id) : eventID(id) {}
  
  // Copy constructor
  Event(const Event& other)
    : eventID(other.eventID),
      particles(other.particles)
  {}
  
  // Assignment operator
  Event& operator=(const Event& other)
  {
    if (this != &other) {
      eventID = other.eventID;
      particles = other.particles;
    }
    return *this;
  }
  
  void AddParticle(const Particle& p) { particles.push_back(p); }
  G4int GetEventID() const { return eventID; }
  
  // Add const version
  const std::vector<Particle>& GetParticles() const { return particles; }
  std::vector<Particle>& GetParticles() { return particles; }
  
private:
  G4int eventID;
  std::vector<Particle> particles;
};

#endif