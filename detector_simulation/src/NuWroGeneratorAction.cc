#include "NuWroGeneratorAction.hh"
#include "AnalysisManager.hh"

#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"

// ROOT includes
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

NuWroGeneratorAction::NuWroGeneratorAction(const G4String& fileName, G4int initialEvent)
 : G4VUserPrimaryGeneratorAction(),
   fAnalysisManager(nullptr),
   fTargetVolume(nullptr),
   fNuWroFile(nullptr),
   fNuWroTree(nullptr),
   fCurrentEvent(0),
   fTotalEvents(0),
   fROOTFileName(fileName),
   fInitialEvent(initialEvent),
   fParticleGun(nullptr)
{

  // Ensure ROOT classes are properly loaded
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libRIO");

  SetTargetVolumeByName("TPC_log");

  // Get analysis manager
  fAnalysisManager = AnalysisManager::GetInstance();

  for (int i = 0; i < kMaxParticles; i++) {
    fPDGCodes[i] = 0;
    fE[i] = 0.0;
    fPx[i] = 0.0;
    fPy[i] = 0.0;
    fPz[i] = 0.0;
  }

  // Initialize particle gun
  fParticleGun = new G4ParticleGun(1);

  // Initialize ROOT file and tree
  Initialize();
  
}

NuWroGeneratorAction::~NuWroGeneratorAction()
{
  // Clean up
  delete fParticleGun;
  
  if (fNuWroFile) {
    fNuWroFile->Close();
    delete fNuWroFile;
  }
}

void NuWroGeneratorAction::Initialize()
{

  // Open NuWro output ROOT file
  fNuWroFile = new TFile(fROOTFileName.c_str(), "READ");
  if (!fNuWroFile || fNuWroFile->IsZombie()) {
    G4cerr << "Error: Could not open NuWro ROOT file " << fROOTFileName << G4endl;
    delete fNuWroFile;
    fNuWroFile = nullptr;
    fTotalEvents = 0;
    return;
  }
  
  // Get NuWro treeout event tree
  // Only need the tree with post-FSI
  fNuWroTree = (TTree*)fNuWroFile->Get("treeout");
  if (!fNuWroTree) {
    G4cerr << "Error: Could not find NuWro event tree in file " << fROOTFileName << G4endl;
    fNuWroFile->Close();
    delete fNuWroFile;
    fNuWroFile = nullptr;
    fTotalEvents = 0;
    return;
  }

  // Get number of events in the file
  fTotalEvents = fNuWroTree->GetEntries();
  G4cout << "NuWro input file contains " << fTotalEvents << " events." << G4endl;
  
  // Set up branch addresses
  fNuWroTree->SetBranchAddress("e.post.pdg",  fPDGCodes);
  fNuWroTree->SetBranchAddress("e.post.t",    fE);
  fNuWroTree->SetBranchAddress("e.post.x",    fPx);
  fNuWroTree->SetBranchAddress("e.post.y",    fPy);
  fNuWroTree->SetBranchAddress("e.post.z",    fPz);
  
  // Reset event counter
  fCurrentEvent = fInitialEvent;
}

bool NuWroGeneratorAction::LoadNextEvent()
{
  
  if (!fNuWroFile || !fNuWroTree) {
    G4cerr << "Error: NuWro ROOT file or tree not available." << G4endl;
    return false;
  }
  
  // Check if we have more events
  if (fCurrentEvent >= fTotalEvents) {
    G4cout << "Warning: Reached the end of NuWro events. Restarting from the beginning." << G4endl;
    fCurrentEvent = fInitialEvent;
  }
  
  // Get next event
  fNuWroTree->GetEntry(fCurrentEvent);
  fCurrentEvent++;
  
  return true;
}

void NuWroGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Load the next event from the NuWro file
  if (!LoadNextEvent()) {
    G4cerr << "Failed to load event from NuWro file. Generating a default particle." << G4endl;
    // Generate a default particle (e.g., 1 GeV muon) as fallback
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(1.0*GeV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    return;
  }
  
  // Create a primary vertex at the neutrino interaction point
  //G4ThreeVector vtxPos(0, 0, 0);
  G4ThreeVector vtxPos = GetRandomPositionInVolume();
  G4double vtxTime = 0.0;
  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(vtxPos, vtxTime);

  //G4cout << "Creating interaction vertex at " << vtxPos/cm << " cm" << G4endl;

  // Get particle definition
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // Create final state particles from the NuWro event
  for (int i = 0; i < kMaxParticles; ++i) {

    if (fPDGCodes[i] == 0) break;
    
    G4ParticleDefinition* particleDef = particleTable->FindParticle(fPDGCodes[i]);
    
    if (!particleDef) {
      G4cerr << "Warning: Unknown PDG code " << fPDGCodes[i] << " skipped." << G4endl;
      continue;
    }
    
    // Give units to energy and momentum
    G4double e  = fE[i]  * MeV;
    G4double px = fPz[i] * MeV;
    G4double py = fPy[i] * MeV;
    G4double pz = fPx[i] * MeV;

    //G4cout << "Particle  " << i << " PDG: " << fPDGCodes[i] << "\n"
    //       << "     E:   " << e/GeV  << " GeV \n"
    //       << "     px:  " << px/GeV << " GeV/c \n"
    //       << "     py:  " << py/GeV << " GeV/c \n"
    //       << "     pz:  " << pz/GeV << " GeV/c \n"
    //       << G4endl;
    
    // Create primary particle
    G4PrimaryParticle* particle = new G4PrimaryParticle(particleDef);
    particle->SetMomentum(px, py, pz);
    particle->SetTotalEnergy(e);
    
    // Add particle to vertex
    vertex->SetPrimary(particle);

  }
  
  // Add vertex to event
  anEvent->AddPrimaryVertex(vertex);
  
  // Print event summary
  //G4cout << "Generated Geant4 event from NuWro event #" << fCurrentEvent-1 << "\n"
  //       << G4endl;
}

void NuWroGeneratorAction::SetTargetVolumeByName(const G4String& volumeName)
{
  // Get the logical volume store
  G4LogicalVolumeStore* logicalVolumeStore = G4LogicalVolumeStore::GetInstance();
  
  // Iterate through all logical volumes to find the one with the matching name
  for (G4LogicalVolume* logicalVolume : *logicalVolumeStore) {
    if (logicalVolume->GetName() == volumeName) {
      fTargetVolume = logicalVolume;
      G4cout << "Target volume for vertex generation set to logical volume: " 
             << volumeName << G4endl;
      return;
    }
  }
  
  // If we get here, the volume wasn't found
  G4cerr << "Error: Could not find logical volume named '" << volumeName 
         << "' for vertex generation." << G4endl;
  fTargetVolume = nullptr;
}

G4ThreeVector NuWroGeneratorAction::GetRandomPositionInVolume()
{
  if (!fTargetVolume) {
    G4cerr << "Error: Target logical volume not set for random vertex generation!" << G4endl;
    return G4ThreeVector(0, 0, 0);
  }
  
  // Get the solid directly from the logical volume
  G4VSolid* solid = fTargetVolume->GetSolid();
  
  // Get the bounding box
  G4ThreeVector min, max;
  solid->BoundingLimits(min, max);
  
  // Generate random position within the bounding box
  G4ThreeVector position;
  G4bool isInside = false;
  
  // Keep generating until we get a point inside the actual solid
  while (!isInside) {
    // Generate random position within the bounding box
    position.setX(min.x() + G4UniformRand() * (max.x() - min.x()));
    position.setY(min.y() + G4UniformRand() * (max.y() - min.y()));
    position.setZ(min.z() + G4UniformRand() * (max.z() - min.z()));
    
    // Check if the point is inside the solid
    isInside = solid->Inside(position);
  }
  
  return position;
}
