#include "AnalysisManager.hh"
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
  fWriteTrajectory(false),
  fWriteTPCHits(true),
  fWriteECalHits(true),
  fWriteMuIDHits(true),
  fParticleNtupleID(-1),
  fTrajectoryNtupleID(-1),
  fTPCHitNtupleID(-1),
  fECalHitNtupleID(-1),
  fMuIDHitNtupleID(-1)
{

  G4cout << "AnalysisManager created!" << G4endl;
  DefineCommands();

}

// Destructor
AnalysisManager::~AnalysisManager()
{
  delete fCurrentEvent;
  for (auto& pair : fTrackMap) {
    delete pair.second;
  }
  delete fMessenger;
}

// Initialize analysis - call at the beginning of the run
void AnalysisManager::Book()
{
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  // Default settings
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName(fOutputFileName);
  
  // Create ntuple structure - this only needs to be done once
  CreateNtuples();
  
  // Open output file - G4AnalysisManager handles the threading details
  analysisManager->OpenFile();
}

// Create all ntuples
void AnalysisManager::CreateNtuples()
{
  auto analysisManager = G4AnalysisManager::Instance();
  
  // Particle ntuple
  fParticleNtupleID = analysisManager->CreateNtuple("Particles", "Particle Information");
  fParticleNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
  fParticleNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
  fParticleNtuple.pdgCode = analysisManager->CreateNtupleIColumn("PDGCode");
  fParticleNtuple.nTrajectoryPoints = analysisManager->CreateNtupleIColumn("NTrajectoryPoints");
  fParticleNtuple.nTPCHits = analysisManager->CreateNtupleIColumn("NTPCHits");
  fParticleNtuple.nECalHits = analysisManager->CreateNtupleIColumn("NECalHits");
  fParticleNtuple.processName = analysisManager->CreateNtupleSColumn("ProcessName");
  analysisManager->FinishNtuple();
  
  // Trajectory ntuple
  if (fWriteTrajectory) {
    fTrajectoryNtupleID = analysisManager->CreateNtuple("Trajectories", "Trajectory Points");
    fTrajectoryNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fTrajectoryNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fTrajectoryNtuple.pointIndex = analysisManager->CreateNtupleIColumn("PointIndex");
    fTrajectoryNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fTrajectoryNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fTrajectoryNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fTrajectoryNtuple.time = analysisManager->CreateNtupleDColumn("Time");
    analysisManager->FinishNtuple();
  }
  
  // TPC Hit ntuple
  if (fWriteTPCHits) {
    fTPCHitNtupleID = analysisManager->CreateNtuple("TPCHits", "TPC Hit Information");
    fTPCHitNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fTPCHitNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fTPCHitNtuple.hitIndex = analysisManager->CreateNtupleIColumn("HitIndex");
    fTPCHitNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fTPCHitNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fTPCHitNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fTPCHitNtuple.energyDeposit = analysisManager->CreateNtupleDColumn("EnergyDeposit");
    fTPCHitNtuple.stepSize = analysisManager->CreateNtupleDColumn("StepSize");
    analysisManager->FinishNtuple();

    fTPCSecHitNtupleID = analysisManager->CreateNtuple("TPCSecHits", "TPC Secondary Hit Information");
    fTPCHitNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fTPCHitNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fTPCHitNtuple.hitIndex = analysisManager->CreateNtupleIColumn("HitIndex");
    fTPCHitNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fTPCHitNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fTPCHitNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fTPCHitNtuple.energyDeposit = analysisManager->CreateNtupleDColumn("EnergyDeposit");
    fTPCHitNtuple.stepSize = analysisManager->CreateNtupleDColumn("StepSize");
    analysisManager->FinishNtuple();
  }
  
  // ECal Hit ntuple
  if (fWriteECalHits) {
    fECalHitNtupleID = analysisManager->CreateNtuple("ECalHits", "ECal Hit Information");
    fECalHitNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fECalHitNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fECalHitNtuple.hitIndex = analysisManager->CreateNtupleIColumn("HitIndex");
    fECalHitNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fECalHitNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fECalHitNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fECalHitNtuple.energyDeposit = analysisManager->CreateNtupleDColumn("EnergyDeposit");
    analysisManager->FinishNtuple();

    fECalSecHitNtupleID = analysisManager->CreateNtuple("ECalSecHits", "ECal Secondary Hit Information");
    fECalHitNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fECalHitNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fECalHitNtuple.hitIndex = analysisManager->CreateNtupleIColumn("HitIndex");
    fECalHitNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fECalHitNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fECalHitNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fECalHitNtuple.energyDeposit = analysisManager->CreateNtupleDColumn("EnergyDeposit");
    analysisManager->FinishNtuple();
  }

  // MuID Hit ntuple
  if (fWriteMuIDHits) {
    fMuIDHitNtupleID = analysisManager->CreateNtuple("MuIDHits", "MuID Hit Information");
    fMuIDHitNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fMuIDHitNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fMuIDHitNtuple.hitIndex = analysisManager->CreateNtupleIColumn("HitIndex");
    fMuIDHitNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fMuIDHitNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fMuIDHitNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fMuIDHitNtuple.energyDeposit = analysisManager->CreateNtupleDColumn("EnergyDeposit");
    analysisManager->FinishNtuple();

    fMuIDSecHitNtupleID = analysisManager->CreateNtuple("MuIDSecHits", "MuID Secondary Hit Information");
    fMuIDHitNtuple.eventID = analysisManager->CreateNtupleIColumn("EventID");
    fMuIDHitNtuple.trackID = analysisManager->CreateNtupleIColumn("TrackID");
    fMuIDHitNtuple.hitIndex = analysisManager->CreateNtupleIColumn("HitIndex");
    fMuIDHitNtuple.posX = analysisManager->CreateNtupleDColumn("PositionX");
    fMuIDHitNtuple.posY = analysisManager->CreateNtupleDColumn("PositionY");
    fMuIDHitNtuple.posZ = analysisManager->CreateNtupleDColumn("PositionZ");
    fMuIDHitNtuple.energyDeposit = analysisManager->CreateNtupleDColumn("EnergyDeposit");
    analysisManager->FinishNtuple();
  }
  
  // You can add more ntuples for other data types here

}

// Save data - call at the end of the run
void AnalysisManager::Save()
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
}

// Close file - call at the end of the run after Save()
void AnalysisManager::Close()
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->CloseFile();
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
  fCurrentEvent = new Event(eventID);

}

// End current event and write to ntuple
void AnalysisManager::EndEvent()
{

  G4int eventID = fCurrentEvent->GetEventID();
  //G4cout << "Finishing event " << eventID << G4endl;

  std::map<G4int, Particle*> ParticleMap = SimplifyParticleCollection();

  if (fCurrentEvent) {
    // Write Event data to ntuples
    for (const auto& pair : ParticleMap) {
      fCurrentEvent->AddParticle(*pair.second);
    }
    
    WriteEvent(*fCurrentEvent);
  
  }
}

// Get Particle object for a given track ID
Particle* AnalysisManager::GetParticle(const G4Track* track)
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
  
  Particle* particle = new Particle(trackID, pdgCode, creatorProcess, momID);
  fTrackMap[trackID] = particle;
  return particle;
}

// Add a TPC hit to a particle
void AnalysisManager::AddTPCHit(const G4Track* track, const G4ThreeVector& pos, G4double edep, G4double step)
{
  Particle* particle = GetParticle(track);
  particle->AddTPCHit(TPCHit(pos, edep, step));
}

// Add an ECal hit to a particle
void AnalysisManager::AddECalHit(const G4Track* track, const G4ThreeVector& pos, G4double edep)
{
  Particle* particle = GetParticle(track);
  particle->AddECalHit(ECalHit(pos, edep));
}

// Add a MuID hit to a particle
void AnalysisManager::AddMuIDHit(const G4Track* track, const G4ThreeVector& pos, G4double edep)
{
  Particle* particle = GetParticle(track);
  particle->AddMuIDHit(MuIDHit(pos, edep));
}

// Record energy deposit data from a step
void AnalysisManager::RecordEnergyDeposit(const G4Step* step)
{
  // Skip if no energy was deposited
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= fEnergyCut) return;
  //G4cout << "Edep: " << edep/MeV << " MeV" << G4endl;
  
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
  G4double time = track->GetGlobalTime();
  G4String volumeName = track->GetVolume()->GetLogicalVolume()->GetName();
  
  // Update the PDG code for the particle
  Particle* particle = GetParticle(track);
  
  // Add trajectory point
  particle->GetTrajectory().AddPoint(position, time, volumeName);

}

std::map<G4int, Particle*> AnalysisManager::SimplifyParticleCollection()
{

  // Vector with TrackIDs of particles to write
  std::vector<G4int> particles_to_write;

  for (const auto& pair : fTrackMap) {
    Particle* particle = pair.second;
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
  std::map<G4int, Particle*> particle_map;
  // Create TrackID map to climb particle tree
  std::map<G4int, G4int> daughter_mother_map;

  for (const auto& pair : fTrackMap) {
    Particle* particle = pair.second;

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

void AnalysisManager::WriteEvent(const Event& event)
{
  G4int eventID = event.GetEventID();
  
  for (const auto& particle : event.GetParticles()) {
    WriteParticle(particle, eventID);
  }
}

// Write a particle and all its data to ntuples
void AnalysisManager::WriteParticle(const Particle& particle, G4int eventID)
{
  auto analysisManager = G4AnalysisManager::Instance();
  G4int trackID = particle.GetTrackID();
  
  // Write particle information
  analysisManager->FillNtupleIColumn(fParticleNtupleID, fParticleNtuple.eventID, eventID);
  analysisManager->FillNtupleIColumn(fParticleNtupleID, fParticleNtuple.trackID, trackID);
  analysisManager->FillNtupleIColumn(fParticleNtupleID, fParticleNtuple.pdgCode, particle.GetPDGCode());
  analysisManager->FillNtupleIColumn(fParticleNtupleID, fParticleNtuple.nTrajectoryPoints, 
                                    particle.GetTrajectory().GetPoints().size());
  analysisManager->FillNtupleIColumn(fParticleNtupleID, fParticleNtuple.nTPCHits, 
                                    particle.GetTPCHits().size());
  analysisManager->FillNtupleIColumn(fParticleNtupleID, fParticleNtuple.nECalHits, 
                                    particle.GetECalHits().size());
  analysisManager->FillNtupleSColumn(fParticleNtupleID, fParticleNtuple.processName, particle.GetCreatorProcess());
  analysisManager->AddNtupleRow(fParticleNtupleID);
  
  // Write trajectory points
  if (fWriteTrajectory) {
    const auto& trajectoryPoints = particle.GetTrajectory().GetPoints();
    for (size_t i = 0; i < trajectoryPoints.size(); ++i) {
      const auto& point = trajectoryPoints[i];
      analysisManager->FillNtupleIColumn(fTrajectoryNtupleID, fTrajectoryNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fTrajectoryNtupleID, fTrajectoryNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fTrajectoryNtupleID, fTrajectoryNtuple.pointIndex, i);
      analysisManager->FillNtupleDColumn(fTrajectoryNtupleID, fTrajectoryNtuple.posX, point.position.x()/cm);
      analysisManager->FillNtupleDColumn(fTrajectoryNtupleID, fTrajectoryNtuple.posY, point.position.y()/cm);
      analysisManager->FillNtupleDColumn(fTrajectoryNtupleID, fTrajectoryNtuple.posZ, point.position.z()/cm);
      analysisManager->FillNtupleDColumn(fTrajectoryNtupleID, fTrajectoryNtuple.time, point.time/ns);
      analysisManager->AddNtupleRow(fTrajectoryNtupleID);
    }
  }
  
  // Write TPC hits
  if (fWriteTPCHits) {
    const auto& tpcHits = particle.GetTPCHits();
    for (size_t i = 0; i < tpcHits.size(); ++i) {
      const auto& hit = tpcHits[i];
      analysisManager->FillNtupleIColumn(fTPCHitNtupleID, fTPCHitNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fTPCHitNtupleID, fTPCHitNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fTPCHitNtupleID, fTPCHitNtuple.hitIndex, i);
      analysisManager->FillNtupleDColumn(fTPCHitNtupleID, fTPCHitNtuple.posX, hit.GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(fTPCHitNtupleID, fTPCHitNtuple.posY, hit.GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(fTPCHitNtupleID, fTPCHitNtuple.posZ, hit.GetPosition().z()/cm);
      analysisManager->FillNtupleDColumn(fTPCHitNtupleID, fTPCHitNtuple.energyDeposit, hit.GetEnergyDeposit()/MeV);
      analysisManager->FillNtupleDColumn(fTPCHitNtupleID, fTPCHitNtuple.stepSize, hit.GetStepSize()/cm);
      analysisManager->AddNtupleRow(fTPCHitNtupleID);
    }

    const auto& tpcSecHits = particle.GetSecondaryTPCHits();
    for (size_t i = 0; i < tpcSecHits.size(); ++i) {
      const auto& hit = tpcSecHits[i];
      analysisManager->FillNtupleIColumn(fTPCSecHitNtupleID, fTPCHitNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fTPCSecHitNtupleID, fTPCHitNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fTPCSecHitNtupleID, fTPCHitNtuple.hitIndex, i);
      analysisManager->FillNtupleDColumn(fTPCSecHitNtupleID, fTPCHitNtuple.posX, hit.GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(fTPCSecHitNtupleID, fTPCHitNtuple.posY, hit.GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(fTPCSecHitNtupleID, fTPCHitNtuple.posZ, hit.GetPosition().z()/cm);
      analysisManager->FillNtupleDColumn(fTPCSecHitNtupleID, fTPCHitNtuple.energyDeposit, hit.GetEnergyDeposit()/MeV);
      analysisManager->FillNtupleDColumn(fTPCSecHitNtupleID, fTPCHitNtuple.stepSize, hit.GetStepSize()/cm);
      analysisManager->AddNtupleRow(fTPCSecHitNtupleID);
    }
  }
  
  // Write ECal hits
  if (fWriteECalHits) {
    const auto& ecalHits = particle.GetECalHits();
    for (size_t i = 0; i < ecalHits.size(); ++i) {
      const auto& hit = ecalHits[i];
      analysisManager->FillNtupleIColumn(fECalHitNtupleID, fECalHitNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fECalHitNtupleID, fECalHitNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fECalHitNtupleID, fECalHitNtuple.hitIndex, i);
      analysisManager->FillNtupleDColumn(fECalHitNtupleID, fECalHitNtuple.posX, hit.GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(fECalHitNtupleID, fECalHitNtuple.posY, hit.GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(fECalHitNtupleID, fECalHitNtuple.posZ, hit.GetPosition().z()/cm);
      analysisManager->FillNtupleDColumn(fECalHitNtupleID, fECalHitNtuple.energyDeposit, hit.GetEnergyDeposit()/MeV);
      analysisManager->AddNtupleRow(fECalHitNtupleID);
    }

    const auto& ecalSecHits = particle.GetSecondaryECalHits();
    for (size_t i = 0; i < ecalSecHits.size(); ++i) {
      const auto& hit = ecalSecHits[i];
      analysisManager->FillNtupleIColumn(fECalSecHitNtupleID, fECalHitNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fECalSecHitNtupleID, fECalHitNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fECalSecHitNtupleID, fECalHitNtuple.hitIndex, i);
      analysisManager->FillNtupleDColumn(fECalSecHitNtupleID, fECalHitNtuple.posX, hit.GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(fECalSecHitNtupleID, fECalHitNtuple.posY, hit.GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(fECalSecHitNtupleID, fECalHitNtuple.posZ, hit.GetPosition().z()/cm);
      analysisManager->FillNtupleDColumn(fECalSecHitNtupleID, fECalHitNtuple.energyDeposit, hit.GetEnergyDeposit()/MeV);
      analysisManager->AddNtupleRow(fECalSecHitNtupleID);
    }
  }

  // Write MuID hits
  if (fWriteMuIDHits) {
    const auto& muidHits = particle.GetMuIDHits();
    for (size_t i = 0; i < muidHits.size(); ++i) {
      const auto& hit = muidHits[i];
      analysisManager->FillNtupleIColumn(fMuIDHitNtupleID, fMuIDHitNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fMuIDHitNtupleID, fMuIDHitNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fMuIDHitNtupleID, fMuIDHitNtuple.hitIndex, i);
      analysisManager->FillNtupleDColumn(fMuIDHitNtupleID, fMuIDHitNtuple.posX, hit.GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(fMuIDHitNtupleID, fMuIDHitNtuple.posY, hit.GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(fMuIDHitNtupleID, fMuIDHitNtuple.posZ, hit.GetPosition().z()/cm);
      analysisManager->FillNtupleDColumn(fMuIDHitNtupleID, fMuIDHitNtuple.energyDeposit, hit.GetEnergyDeposit()/MeV);
      analysisManager->AddNtupleRow(fMuIDHitNtupleID);
    }

    const auto& muidSecHits = particle.GetSecondaryMuIDHits();
    for (size_t i = 0; i < muidSecHits.size(); ++i) {
      const auto& hit = muidSecHits[i];
      analysisManager->FillNtupleIColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.eventID, eventID);
      analysisManager->FillNtupleIColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.trackID, trackID);
      analysisManager->FillNtupleIColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.hitIndex, i);
      analysisManager->FillNtupleDColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.posX, hit.GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.posY, hit.GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.posZ, hit.GetPosition().z()/cm);
      analysisManager->FillNtupleDColumn(fMuIDSecHitNtupleID, fMuIDHitNtuple.energyDeposit, hit.GetEnergyDeposit()/MeV);
      analysisManager->AddNtupleRow(fMuIDSecHitNtupleID);
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
