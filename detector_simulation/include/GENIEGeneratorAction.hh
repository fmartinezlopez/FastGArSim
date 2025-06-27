#ifndef GENIEGENERATORACTION_HH
#define GENIEGENERATORACTION_HH

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

class GENIEGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GENIEGeneratorAction(const G4String& fileName, G4int initialEvent);
  virtual ~GENIEGeneratorAction();
  
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

  G4int DetermineLeptonPDG(G4int nuPDG, G4bool isCC) const;
  G4ThreeVector GetRandomPositionInVolume();

  // GENIE event tree
  TFile*   fGenieFile;
  TTree*   fGenieTree;
  
  // Event counters
  G4int fInitialEvent;
  G4int fCurrentEvent;
  G4int fTotalEvents;
  
  // GENIE ROOT file name
  G4String fROOTFileName;
  
  // Particle gun for each primary particle
  G4ParticleGun* fParticleGun;

  static const int kMaxParticles = 100;
  
  // Variables to store GENIE event info
  int fNu;                               // Neutrino PDG code
  bool fCC;                              // Is it a CC event?
  double fELepton;                       // Energy of final state primary lepton (GeV)
  double fPxLepton;                      // Momentum of final state primary lepton (GeV/c)
  double fPyLepton;
  double fPzLepton;
  int fNFSI;                             // Number of final state particles in hadronic system
  Int_t fPDGCodes[kMaxParticles];        // PDG codes of final state particles in hadronic system
  Double_t fEHadron[kMaxParticles];      // Energy of final state particles in hadronic system (GeV)
  Double_t fPxHadron[kMaxParticles];     // Momenta of final state particle in hadronic system (GeV/c)
  Double_t fPyHadron[kMaxParticles];
  Double_t fPzHadron[kMaxParticles];
  double fVtxX;                          // Vertex position in detector coordinate system
  double fVtxY;
  double fVtxZ;
  double fVtxT;                          // Vertex time in detector coordinate system

};

#endif