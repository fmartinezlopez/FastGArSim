#include "GENIEGeneratorAction.hh"
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

GENIEGeneratorAction::GENIEGeneratorAction(const G4String& fileName)
 : G4VUserPrimaryGeneratorAction(),
   fAnalysisManager(nullptr),
   fTargetVolume(nullptr),
   fGenieFile(nullptr),
   fGenieTree(nullptr),
   fCurrentEvent(0),
   fTotalEvents(0),
   fROOTFileName(fileName),
   fParticleGun(nullptr),
   fNu(0), fCC(false),
   fELepton(0.0), fPxLepton(0.0), fPyLepton(0.0), fPzLepton(0.0),
   fNFSI(0),
   //fPDGCodes(nullptr),
   //fEHadron(nullptr), fPxHadron(nullptr), fPyHadron(nullptr), fPzHadron(nullptr),
   fVtxX(0.0), fVtxY(0.0), fVtxZ(0.0), fVtxT(0.0)
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
    fEHadron[i] = 0.0;
    fPxHadron[i] = 0.0;
    fPyHadron[i] = 0.0;
    fPzHadron[i] = 0.0;
  }

  // Initialize particle gun
  fParticleGun = new G4ParticleGun(1);

  // Initialize ROOT file and tree
  Initialize();
  
}

GENIEGeneratorAction::~GENIEGeneratorAction()
{
  // Clean up
  delete fParticleGun;
  
  if (fGenieFile) {
    fGenieFile->Close();
    delete fGenieFile;
  }
}

void GENIEGeneratorAction::Initialize()
{

  // Open GENIE output ROOT file
  fGenieFile = new TFile(fROOTFileName.c_str(), "READ");
  if (!fGenieFile || fGenieFile->IsZombie()) {
    G4cerr << "Error: Could not open GENIE ROOT file " << fROOTFileName << G4endl;
    delete fGenieFile;
    fGenieFile = nullptr;
    fTotalEvents = 0;
    return;
  }
  
  // Get GENIE gst event tree
  fGenieTree = (TTree*)fGenieFile->Get("gst");
  if (!fGenieTree) {
    G4cerr << "Error: Could not find GENIE event tree in file " << fROOTFileName << G4endl;
    fGenieFile->Close();
    delete fGenieFile;
    fGenieFile = nullptr;
    fTotalEvents = 0;
    return;
  }

  // Get number of events in the file
  fTotalEvents = fGenieTree->GetEntries();
  G4cout << "GENIE input file contains " << fTotalEvents << " events." << G4endl;
  
  // Set up branch addresses
  fGenieTree->SetBranchAddress("neu",   &fNu);
  fGenieTree->SetBranchAddress("cc",    &fCC);
  fGenieTree->SetBranchAddress("El",    &fELepton);
  fGenieTree->SetBranchAddress("pxl",   &fPxLepton);
  fGenieTree->SetBranchAddress("pyl",   &fPyLepton);
  fGenieTree->SetBranchAddress("pzl",   &fPzLepton);
  fGenieTree->SetBranchAddress("nf",    &fNFSI);
  fGenieTree->SetBranchAddress("pdgf",  fPDGCodes);
  fGenieTree->SetBranchAddress("Ef",    fEHadron);
  fGenieTree->SetBranchAddress("pxf",   fPxHadron);
  fGenieTree->SetBranchAddress("pyf",   fPyHadron);
  fGenieTree->SetBranchAddress("pzf",   fPzHadron);
  fGenieTree->SetBranchAddress("vtxx",  &fVtxX);
  fGenieTree->SetBranchAddress("vtxy",  &fVtxY);
  fGenieTree->SetBranchAddress("vtxz",  &fVtxZ);
  fGenieTree->SetBranchAddress("vtxt",  &fVtxT);
  
  // Reset event counter
  fCurrentEvent = 0;
}

bool GENIEGeneratorAction::LoadNextEvent()
{
  
  if (!fGenieFile || !fGenieTree) {
    G4cerr << "Error: GENIE ROOT file or tree not available." << G4endl;
    return false;
  }
  
  // Check if we have more events
  if (fCurrentEvent >= fTotalEvents) {
    G4cout << "Warning: Reached the end of GENIE events. Restarting from the beginning." << G4endl;
    fCurrentEvent = 0;
  }
  
  // Get next event
  fGenieTree->GetEntry(fCurrentEvent);
  fCurrentEvent++;
  
  return true;
}

void GENIEGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Load the next event from the GENIE file
  if (!LoadNextEvent()) {
    G4cerr << "Failed to load event from GENIE file. Generating a default particle." << G4endl;
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
  G4double vtxTime = 0;

  //G4cout << "Creating interaction vertex at " << vtxPos/m << G4endl;
  
  G4PrimaryVertex* vertex = new G4PrimaryVertex(vtxPos, vtxTime);

  // Get particle definition
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // Create lepton from the GENIE event
  G4int leptonPDG = DetermineLeptonPDG(fNu, fCC);
  G4ParticleDefinition* leptonDef = particleTable->FindParticle(leptonPDG);

  // Give units to energy and momentum
  G4double el  = fELepton  * GeV;
  G4double pxl = fPzLepton * GeV;
  G4double pyl = fPyLepton * GeV;
  G4double pzl = fPxLepton * GeV;

  //G4cout << "Lepton PDG: " << leptonPDG << "\n"
  //       << "       E:   " << el/GeV    << " GeV \n"
  //       << "       px:  " << pxl/GeV   << " GeV/c \n"
  //       << "       py:  " << pyl/GeV   << " GeV/c \n"
  //       << "       pz:  " << pzl/GeV   << " GeV/c \n"
  //       << G4endl;

  G4PrimaryParticle* lepton = new G4PrimaryParticle(leptonDef);
  lepton->SetMomentum(pxl, pyl, pzl);
  lepton->SetTotalEnergy(el);
  
  // Add lepton to vertex
  vertex->SetPrimary(lepton);

  //G4cout << "Number of FSI hadrons:  " << fNFSI << G4endl;

  // Create final state hadrons from the GENIE event
  for (int i = 0; i < fNFSI; ++i) {
    
    G4ParticleDefinition* hadronDef = particleTable->FindParticle(fPDGCodes[i]);
    
    if (!hadronDef) {
      G4cerr << "Warning: Unknown PDG code " << fPDGCodes[i] << " skipped." << G4endl;
      continue;
    }
    
    // Give units to energy and momentum
    G4double ef  = fEHadron[i]  * GeV;
    G4double pxf = fPzHadron[i] * GeV;
    G4double pyf = fPyHadron[i] * GeV;
    G4double pzf = fPxHadron[i] * GeV;

    //G4cout << "Hadron " << i << " PDG: " << fPDGCodes[i] << "\n"
    //       << "     E:   " << ef/GeV  << " GeV \n"
    //       << "     px:  " << pxf/GeV << " GeV/c \n"
    //       << "     py:  " << pyf/GeV << " GeV/c \n"
    //       << "     pz:  " << pzf/GeV << " GeV/c \n"
    //       << G4endl;
    
    // Create primary particle
    G4PrimaryParticle* hadron = new G4PrimaryParticle(hadronDef);
    hadron->SetMomentum(pxf, pyf, pzf);
    hadron->SetTotalEnergy(ef);
    
    // Add particle to vertex
    vertex->SetPrimary(hadron);

  }
  
  // Add vertex to event
  anEvent->AddPrimaryVertex(vertex);
  
  // Print event summary
  //G4cout << "Generated Geant4 event from GENIE event #" << fCurrentEvent-1 << "\n"
  //       << "Nu flavour: " << fNu << "\n"
  //       << (fCC ? "CC" : "NC") << "\n"
  //       << G4endl;
}

void GENIEGeneratorAction::SetTargetVolumeByName(const G4String& volumeName)
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

G4ThreeVector GENIEGeneratorAction::GetRandomPositionInVolume()
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

G4int GENIEGeneratorAction::DetermineLeptonPDG(G4int nuPDG, G4bool isCC) const
{
  // Check if this is a NC or CC interaction
  if (!isCC) {
    // For neutral current (NC) interactions, the outgoing lepton is a neutrino of the same type
    return nuPDG;
  }
  
  // For charged current (CC) interactions:
  // ve->e-, vmu->mu-, vtau->tau-
  // Anti-neutrinos produce positive leptons
  switch (std::abs(nuPDG)) {
    case 12: // electron neutrino
      return (nuPDG > 0) ? 11 : -11;  // e- or e+
    case 14: // muon neutrino
      return (nuPDG > 0) ? 13 : -13;  // mu- or mu+
    case 16: // tau neutrino
      return (nuPDG > 0) ? 15 : -15;  // tau- or tau+
    default:
      G4cerr << "Warning: Unknown neutrino PDG code " << nuPDG << G4endl;
      return 0;
  }
}