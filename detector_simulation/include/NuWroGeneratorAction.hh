#ifndef NUWROGENERATORACTION_HH
#define NUWROGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include <vector>
#include <string>

// ROOT includes
#include "Rtypes.h"
#include "TString.h"

// Forward declarations for ROOT classes
class TFile;
class TTree;

class AnalysisManager;

class NuWroGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  NuWroGeneratorAction(const G4String& fileName);
  virtual ~NuWroGeneratorAction();
  
  virtual void GeneratePrimaries(G4Event*);
  
  // Method to get the current event number
  G4int GetCurrentEventNumber() const { return fCurrentEvent; }
  
  // Get the total number of events in the ROOT file
  G4int GetTotalEvents() const { return fTotalEvents; }

  void SetTargetVolumeByName(const G4String& volumeName);
  
private:

  AnalysisManager* fAnalysisManager;
  
  // Initialize ROOT file and tree
  void Initialize();
  
  // Load next event from ROOT file
  bool LoadNextEvent();

  // Logical volume to generate vertices in
  G4LogicalVolume* fTargetVolume;

  G4ThreeVector GetRandomPositionInVolume();

  // NuWro event tree
  TFile*   fNuWroFile;
  TTree*   fNuWroTree;
  
  // Event counters
  G4int fCurrentEvent;
  G4int fTotalEvents;
  
  // NuWro ROOT file name
  G4String fROOTFileName;
  
  // Particle gun for each primary particle
  G4ParticleGun* fParticleGun;

  static const int kMaxParticles = 100;
  
  // Variables to store NuWro event info
  Int_t fPDGCodes[kMaxParticles];  // PDG codes of final state particles
  Double_t fE[kMaxParticles];      // Energy of final state particles (MeV)
  Double_t fPx[kMaxParticles];     // Momenta of final state particle (MeV/c)
  Double_t fPy[kMaxParticles];
  Double_t fPz[kMaxParticles];

};

#endif