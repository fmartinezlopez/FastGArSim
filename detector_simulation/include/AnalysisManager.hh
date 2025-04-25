#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "DataTypes.hh"

#include "G4AnalysisManager.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"

class AnalysisManager
{
public:
  static AnalysisManager* GetInstance();
  virtual ~AnalysisManager();
  
  // Initialization methods for the beginning and end of run
  void Book();
  void Save();
  void Close();
  
  // Methods to record data
  void RecordEnergyDeposit(const G4Step* step);
  void RecordTrackInfo(const G4Track* track);

  // Methods for Event/Particle-level data
  void BeginEvent(G4int eventID);
  void EndEvent();
  void AddTPCHit(const G4Track* track, const G4ThreeVector& pos, G4double edep, G4double step);
  void AddECalHit(const G4Track* track, const G4ThreeVector& pos, G4double edep);
  void AddMuIDHit(const G4Track* track, const G4ThreeVector& pos, G4double edep);
  Particle* GetParticle(const G4Track* track);

  // Setter methods
  void SetOutputFileName(const G4String& name) { fOutputFileName = name; }
  void SetEnergyCut(G4double cut);

private:
  AnalysisManager();
  static G4ThreadLocal AnalysisManager* fInstance;
  
  // Helper methods
  void CreateNtuples();
  std::map<G4int, Particle*> SimplifyParticleCollection();
  void WriteEvent(const Event& event);
  void WriteParticle(const Particle& particle, G4int eventID);

  void DefineCommands();
  
  // Data members
  G4String fOutputFileName;
  G4double fEnergyCut;

  // Current event data
  Event* fCurrentEvent;
  std::map<G4int, Particle*> fTrackMap;  // Maps trackID to Particle

  //
  G4bool fWriteTrajectory;
  G4bool fWriteTPCHits;
  G4bool fWriteECalHits;
  G4bool fWriteMuIDHits;

  // Ntuple IDs
  G4int fParticleNtupleID;
  G4int fTrajectoryNtupleID;
  G4int fTPCHitNtupleID;
  G4int fTPCSecHitNtupleID;
  G4int fECalHitNtupleID;
  G4int fECalSecHitNtupleID;
  G4int fMuIDHitNtupleID;
  G4int fMuIDSecHitNtupleID;
  
  // Ntuple column indices
  struct {
    G4int eventID;
    G4int trackID;
    G4int pdgCode;
    G4int nTrajectoryPoints;
    G4int nTPCHits;
    G4int nECalHits;
    G4int processName;
  } fParticleNtuple;
  
  struct {
    G4int eventID;
    G4int trackID;
    G4int pointIndex;
    G4int posX;
    G4int posY;
    G4int posZ;
    G4int time;
  } fTrajectoryNtuple;
  
  struct {
    G4int eventID;
    G4int trackID;
    G4int hitIndex;
    G4int posX;
    G4int posY;
    G4int posZ;
    G4int energyDeposit;
    G4int stepSize;
  } fTPCHitNtuple;
  
  struct {
    G4int eventID;
    G4int trackID;
    G4int hitIndex;
    G4int posX;
    G4int posY;
    G4int posZ;
    G4int energyDeposit;
  } fECalHitNtuple;

  struct {
    G4int eventID;
    G4int trackID;
    G4int hitIndex;
    G4int posX;
    G4int posY;
    G4int posZ;
    G4int energyDeposit;
  } fMuIDHitNtuple;

  // Messenger for macro commands
  G4GenericMessenger* fMessenger;

};

#endif