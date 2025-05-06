#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "DataTypes.hh"
#include "ROOTDataTypes.hh"

#include "G4AnalysisManager.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

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

  // Getter methods
  G4bool GetWriteTrajectory() const { return fWriteTrajectory; }

private:
  AnalysisManager();
  static G4ThreadLocal AnalysisManager* fInstance;
  
  // Helper methods
  std::map<G4int, Particle*> SimplifyParticleCollection();
  void WriteEvent(const Event& g4Event);

  void DefineCommands();
  
  // Data members
  G4String fOutputFileName;
  G4double fEnergyCut;

  // Current event data (Geant4 types)
  Event* fCurrentEvent;
  std::map<G4int, Particle*> fTrackMap;  // Maps trackID to Particle

  // ROOT file and tree
  TFile* fRootFile;
  TTree* fEventTree;
  root::Event* fStoredEvent;  // ROOT Event object that will be stored in the tree
  
  // Configuration
  G4bool fWriteTrajectory;
  G4bool fWriteTPCHits;
  G4bool fWriteECalHits;
  G4bool fWriteMuIDHits;

  // Messenger for macro commands
  G4GenericMessenger* fMessenger;
};

#endif