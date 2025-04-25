#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr),
  fMessenger(nullptr),
  fParticleName("mu-"),
  fEnergy(1.0*GeV),
  fPositionZ(0.0),
  fMomentumDirection(G4ThreeVector(1.0, 0.0, 0.0))
{
    // Create particle gun
    fParticleGun = new G4ParticleGun(1);
    
    // Set default particle
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(fParticleName);
    fParticleGun->SetParticleDefinition(particle);
    
    // Set default kinematic properties
    fParticleGun->SetParticleEnergy(fEnergy);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, fPositionZ));
    fParticleGun->SetParticleMomentumDirection(fMomentumDirection);
    
    // Define command interface
    DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{

    // Set default particle
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(fParticleName);
    fParticleGun->SetParticleDefinition(particle);
    
    // Set default kinematic properties
    fParticleGun->SetParticleEnergy(fEnergy);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, fPositionZ));
    fParticleGun->SetParticleMomentumDirection(fMomentumDirection);

    // Generate primary particles
    fParticleGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::SetParticleType(G4String particleName)
{
    // Find the particle by name
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
    
    if (particle) {
        fParticleName = particleName;
        fParticleGun->SetParticleDefinition(particle);
        G4cout << "Particle type set to " << particleName << G4endl;
    } else {
        G4cout << "WARNING: Particle " << particleName << " not found!" << G4endl;
    }
}

void PrimaryGeneratorAction::SetParticleEnergy(G4double energy)
{
    fEnergy = energy;
    fParticleGun->SetParticleEnergy(fEnergy);
    G4cout << "Particle energy set to " << G4BestUnit(fEnergy, "Energy") << G4endl;
}

void PrimaryGeneratorAction::SetPositionZ(G4double z)
{
    fPositionZ = z;
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, fPositionZ));
    G4cout << "Particle position Z set to " << G4BestUnit(fPositionZ, "Length") << G4endl;
}

void PrimaryGeneratorAction::SetMomentumDirection(G4ThreeVector direction)
{
    fMomentumDirection = direction.unit();
    fParticleGun->SetParticleMomentumDirection(fMomentumDirection);
    G4cout << "Particle momentum direction set to (" 
           << fMomentumDirection.x() << ", " 
           << fMomentumDirection.y() << ", " 
           << fMomentumDirection.z() << ")" << G4endl;
}

void PrimaryGeneratorAction::DefineCommands()
{
    // Create command directory
    fMessenger = new G4GenericMessenger(this, "/gun/", "Primary Generator Control");
    
    // Define particle type command
    fMessenger->DeclareProperty("particleType", fParticleName,
                                "Set particle type (e-, e+, mu-, mu+, pi+, proton, etc.)");
    fMessenger->DeclarePropertyWithUnit("energy", "GeV", fEnergy,
                                        "Set particle energy");
    fMessenger->DeclarePropertyWithUnit("positionZ", "cm", fPositionZ,
                                        "Set particle starting position Z coordinate");
    
    // Define momentum direction command with vector parameters
    auto& dirCmd = fMessenger->DeclareProperty("momentumDirection", fMomentumDirection,
                                         "Set particle momentum direction (x,y,z)");
    dirCmd.SetParameterName("x", "y", "z", true);
    dirCmd.SetRange("x != 0 || y != 0 || z != 0");
}