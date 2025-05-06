#include "AnalysisManager.hh"
#include "DataTypeConverters.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Event.hh"

// Initialize static instance - thread local
G4ThreadLocal AnalysisManager* AnalysisManager::fInstance = nullptr;

// Get singleton instance
AnalysisManager* AnalysisManager::GetInstance()
{
  if (!fInstance) {
    // Double-check locking pattern with mutex for thread safety during initialization
    static G4Mutex analysisMutex = G4MUTEX_INITIALIZER;
    G4AutoLock lock(&analysisMutex);
    if (!fInstance) {
      fInstance = new AnalysisManager();
    }
  }
  return fInstance;
}

// Constructor
AnalysisManager::AnalysisManager()
: fOutputFileName("output"),
  fEnergyCut(0.001*MeV),
  fCurrentEvent(nullptr),
  fRootFile(nullptr),
  fEventTree(nullptr),
  fStoredEvent(nullptr),
  fWriteTrajectory(false),
  fWriteTPCHits(true),
  fWriteECalHits(true),
  fWriteMuIDHits(true)
{
  G4cout << "AnalysisManager created!" << G4endl;
  DefineCommands();
}

// Destructor
AnalysisManager::~AnalysisManager()
{
  delete fCurrentEvent;
  delete fStoredEvent;
  for (auto& pair : fTrackMap) {
    delete pair.second;
  }
  delete fMessenger;
  
  // ROOT file is closed in Close() method
}

// Initialize analysis - call at the beginning of the run
void AnalysisManager::Book()
{
  // Create ROOT file and tree
  G4String filename = fOutputFileName + ".root";
  fRootFile = new TFile(filename.c_str(), "RECREATE");
  if (!fRootFile || fRootFile->IsZombie()) {
    G4cerr << "Error: Cannot open ROOT file " << filename << G4endl;
    return;
  }
  
  fEventTree = new TTree("Events", "Event Data");
  fStoredEvent = new root::Event();
  fEventTree->Branch("Event", &fStoredEvent);
  
  G4cout << "AnalysisManager: Created ROOT file " << filename << G4endl;
}


// Save data - call at the end of the run
void AnalysisManager::Save()
{
  if (fRootFile && fEventTree) {
    fRootFile->cd();
    fEventTree->Write();
    G4cout << "AnalysisManager: Saved " << fEventTree->GetEntries() 
           << " events to ROOT file" << G4endl;
  }
}

// Close file - call at the end of the run after Save()
void AnalysisManager::Close()
{
  if (fRootFile) {
    fRootFile->Close();
    delete fRootFile;
    fRootFile = nullptr;
    G4cout << "AnalysisManager: Closed ROOT file" << G4endl;
  }
}

// Begin a new event
void AnalysisManager::BeginEvent(G4int eventID)
{
  if (eventID < 10 || (eventID < 100 && eventID%10 == 0) || eventID%100 == 0) {
    G4cout << "Processing event: " << eventID << G4endl;
  }
  
  // Only clean up if we have previous data
  if (fCurrentEvent) {
    delete fCurrentEvent;
    fCurrentEvent = nullptr;
    
    for (auto& pair : fTrackMap) {
      delete pair.second;
    }
    fTrackMap.clear();
  }
  
  // Create new event
  fCurrentEvent = new ::Event(eventID);
}

// End current event and write to ROOT file
void AnalysisManager::EndEvent()
{
  G4int eventID = fCurrentEvent->GetEventID();
  
  std::map<G4int, ::Particle*> ParticleMap = SimplifyParticleCollection();

  if (fCurrentEvent) {
    // Add particles to the event
    for (const auto& pair : ParticleMap) {
      fCurrentEvent->AddParticle(*pair.second);
    }
    
    WriteEvent(*fCurrentEvent);
  }
}

// Get Particle object for a given track ID
::Particle* AnalysisManager::GetParticle(const G4Track* track)
{
  G4int trackID = track->GetTrackID();

  auto it = fTrackMap.find(trackID);
  if (it != fTrackMap.end()) {
    return it->second;
  }

  G4int pdgCode = track->GetDefinition()->GetPDGEncoding();
  G4String creatorProcess = "primary"; // if particle doesn't have a creator process then it's a primary
  G4int momID = track->GetParentID();

  if (track->GetCreatorProcess()) {
    creatorProcess = track->GetCreatorProcess()->GetProcessName();
  }
  
  ::Particle* particle = new ::Particle(trackID, pdgCode, creatorProcess, momID);
  fTrackMap[trackID] = particle;
  return particle;
}

// Add a TPC hit to a particle
void AnalysisManager::AddTPCHit(const G4Track* track, const G4ThreeVector& pos, G4double edep, G4double step)
{
  ::Particle* particle = GetParticle(track);
  particle->AddTPCHit(::TPCHit(pos, edep, step));
}

// Add an ECal hit to a particle
void AnalysisManager::AddECalHit(const G4Track* track, const G4ThreeVector& pos, G4double edep)
{
  ::Particle* particle = GetParticle(track);
  particle->AddECalHit(::ECalHit(pos, edep));
}

// Add a MuID hit to a particle
void AnalysisManager::AddMuIDHit(const G4Track* track, const G4ThreeVector& pos, G4double edep)
{
  ::Particle* particle = GetParticle(track);
  particle->AddMuIDHit(::MuIDHit(pos, edep));
}

// Record energy deposit data from a step
void AnalysisManager::RecordEnergyDeposit(const G4Step* step)
{
  // Skip if no energy was deposited
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= fEnergyCut) return;
  
  // Get volume information
  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (!volume) return;
  
  // Get track and position information
  G4Track* track = step->GetTrack();
  G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
  G4double stepSize = step->GetStepLength();
  
  // Check volume type and record appropriate hit
  G4String volumeName = volume->GetName();
  
  if (volumeName.find("TPC") != std::string::npos) {
    // This is a TPC hit
    AddTPCHit(track, position, edep, stepSize);
  }
  else if (volumeName.find("ECal") != std::string::npos && volumeName.find("Scintillator") != std::string::npos) {
    // This is an ECal hit
    AddECalHit(track, position, edep);
  }
  else if (volumeName.find("MuID") != std::string::npos && volumeName.find("Scintillator") != std::string::npos) {
    // This is a MuID hit
    AddMuIDHit(track, position, edep);
  }
}

// Record track information (example of another data collection method)
void AnalysisManager::RecordTrackInfo(const G4Track* track)
{
  G4ThreeVector position = track->GetPosition();
  G4ThreeVector momentum = track->GetMomentum();
  G4double time = track->GetGlobalTime();
  G4String volumeName = track->GetVolume()->GetLogicalVolume()->GetName();
  
  // Update the PDG code for the particle
  ::Particle* particle = GetParticle(track);
  
  // Add trajectory point
  particle->GetTrajectory().AddPoint(position, time, volumeName);

  // Add momentum point
  particle->GetTrajectory().AddMomPoint(momentum);
}

std::map<G4int, ::Particle*> AnalysisManager::SimplifyParticleCollection()
{
  // Vector with TrackIDs of particles to write
  std::vector<G4int> particles_to_write;

  for (const auto& pair : fTrackMap) {
    ::Particle* particle = pair.second;
    G4bool write = false;

    // Check for primary particles
    G4String process = particle->GetCreatorProcess();
    if (process == "primary") {
      write = true;
    }

    // Check track length in TPC
    G4double tpc_length = particle->GetTPCPathLength()/cm;
    if (tpc_length > 1.0) {
      write = true;
    }

    // Check if particle starts in TPC but stops in ECal
    G4String vol_start = particle->GetTrajectory().GetStartVolName();
    G4String vol_stop  = particle->GetTrajectory().GetStopVolName();
    if (vol_start.find("TPC") != std::string::npos && vol_stop.find("ECal") != std::string::npos) {
      write = true;
    }

    // If particle passes any of the checks add to collection
    if (write) particles_to_write.push_back(pair.first);
  }

  // Create output collection
  std::map<G4int, ::Particle*> particle_map;
  // Create TrackID map to climb particle tree
  std::map<G4int, G4int> daughter_mother_map;

  for (const auto& pair : fTrackMap) {
    ::Particle* particle = pair.second;

    if (std::find(particles_to_write.begin(), particles_to_write.end(), pair.first) != particles_to_write.end()) {
      particle_map[pair.first] = particle->Clone();
    }
    else {
      // Add entry to daughter -> mother map
      G4int trackID = particle->GetTrackID();
      G4int momID = particle->GetMotherID();
      daughter_mother_map[trackID] = momID;

      // Initialize values for search
      G4bool found = false;
      G4int currentID = trackID;
      G4int nIterations = 0;

      while (!found) {
        // Shouldn't happen, but just in case...
        if (nIterations > 1000) {
          G4cerr << "Warning: Particle couldn't be backtracked!" << G4endl;
          break;
        }

        G4int currentMom = daughter_mother_map[currentID];
        if (std::find(particles_to_write.begin(), particles_to_write.end(), currentMom) != particles_to_write.end()) {
          found = true;
          particle_map[currentMom]->AddSecondaryContribution(particle);
        }
        currentID = currentMom;
        nIterations++;
      }
    }
  }

  return particle_map;
}

void AnalysisManager::WriteEvent(const Event& g4Event)
{
  if (fEventTree && fStoredEvent) {
    // Convert Geant4 Event to ROOT Event
    *fStoredEvent = ConvertEvent(g4Event);
    
    // Fill the tree
    fEventTree->Fill();
    
    G4int eventID = g4Event.GetEventID();
    if (eventID < 10 || (eventID < 100 && eventID%10 == 0) || eventID%100 == 0) {
      G4cout << "Event " << eventID << " written to ROOT file with " 
             << g4Event.GetParticles().size() << " particles" << G4endl;
    }
  }
}

void AnalysisManager::SetEnergyCut(G4double cut)
{
    fEnergyCut = cut;
    G4cout << "Energy cut set to " << fEnergyCut/MeV << " MeV" << G4endl;
}

void AnalysisManager::DefineCommands()
{
    // Initialize messenger for macro commands
    fMessenger = new G4GenericMessenger(this, "/analysis/", "Analysis manager commands");

    fMessenger->DeclarePropertyWithUnit("EnergyCut", "MeV", fEnergyCut, "Set energy cut");
}