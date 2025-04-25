#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization(DetectorConstruction* detConstruction);
    virtual ~ActionInitialization();
    
    virtual void BuildForMaster() const;
    virtual void Build() const;
    void UpdatePrimaryGeneratorAction() const;

    // Setter methods
    void SetGeneratorType(const G4String& type) { fGeneratorType = type; }
    void SetGenieFileName(const G4String& fileName) { fGenieFileName = fileName; }
    void SetNuWroFileName(const G4String& fileName) { fNuWroFileName = fileName; }
    
private:
    DetectorConstruction* fDetConstruction;
    G4String fGeneratorType;   // Type of generator: "particle", "genie", or "nuwro"
    G4String fGenieFileName;   // GENIE ROOT file name (if using genie)
    G4String fNuWroFileName;   // NuWro ROOT file name (if using nuwro)
};

#endif