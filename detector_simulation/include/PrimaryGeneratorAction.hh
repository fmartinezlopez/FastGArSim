#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
    
    // Method from base class
    virtual void GeneratePrimaries(G4Event* event);
    
    // Setter methods for particle gun properties
    void SetParticleType(G4String particleName);
    void SetParticleEnergy(G4double energy);
    void SetPositionZ(G4double z);
    void SetMomentumDirection(G4ThreeVector direction);
    
    // Getter methods
    G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    
private:
    G4ParticleGun* fParticleGun;
    G4GenericMessenger* fMessenger;
    
    // Configurable parameters with default values
    G4String fParticleName;
    G4double fEnergy;
    G4double fPositionZ;
    G4ThreeVector fMomentumDirection;
    
    // Helper methods
    void DefineCommands();
};

#endif