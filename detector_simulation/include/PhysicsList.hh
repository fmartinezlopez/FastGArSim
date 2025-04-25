#ifndef PHYSICS_LIST_HH
#define PHYSICS_LIST_HH

#include "G4VModularPhysicsList.hh"
#include "G4GenericMessenger.hh"

class G4VPhysicsConstructor;

class PhysicsList : public G4VModularPhysicsList
{
public:
    PhysicsList();
    virtual ~PhysicsList();
    
    // Mandatory methods
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    virtual void SetCuts();
    
    // Method to add specific physics components
    void AddPhysicsList(const G4String& name);
    
private:
    // Physics lists to include
    G4VPhysicsConstructor* fEmPhysicsList;
    G4VPhysicsConstructor* fHadronPhysicsList;
    G4VPhysicsConstructor* fDecayPhysicsList;
    
    // Messenger for macro commands
    G4GenericMessenger* fMessenger;
    
    // Physics list options
    G4String fEmName;
    G4String fHadronName;
    
    // Configurable parameters
    G4double fCutValue;
};

#endif