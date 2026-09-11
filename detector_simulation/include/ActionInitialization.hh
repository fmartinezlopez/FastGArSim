#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"
#include "G4VStateDependent.hh"
#include "globals.hh"

class DetectorConstruction;
class G4VUserPrimaryGeneratorAction;

// The primary generator is configured from the run macro, which is read long
// after the actions have been registered. How the final configuration reaches
// the generator depends on the run manager:
//
//   - multi-threaded: Build() runs on each worker at the start of a run, so it
//     already sees the finished configuration and nothing else is needed. The
//     master thread must never call SetUserAction() for a generator -- doing so
//     is a fatal G4MTRunManager error.
//   - sequential: Build() runs once, when the actions are registered and before
//     the macro has been read. The generator is therefore rebuilt from Notify()
//     when the state changes at the start of a run, which is the first moment
//     the configuration is known to be complete.
//
class ActionInitialization : public G4VUserActionInitialization,
                             public G4VStateDependent
{
public:
    ActionInitialization(DetectorConstruction* detector);
    virtual ~ActionInitialization();
    
    virtual void BuildForMaster() const;
    virtual void Build() const;

    // Called by G4StateManager on every state change
    virtual G4bool Notify(G4ApplicationState requestedState);

    // Setter methods
    void SetGeneratorType(const G4String& type) { fGeneratorType = type; fGeneratorChanged = true; }
    void SetGenieFileName(const G4String& fileName) { fGenieFileName = fileName; fGeneratorChanged = true; }
    void SetNuWroFileName(const G4String& fileName) { fNuWroFileName = fileName; fGeneratorChanged = true; }
    void SetInitialEvent(const G4int& index) { fInitialEvent = index; fGeneratorChanged = true; }
    
private:
    // Builds the generator described by the current settings, reporting a
    // missing input file and falling back to the particle gun
    G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction() const;

    DetectorConstruction* fDetectorConstruction;

    G4String fGeneratorType;   // Type of generator: "particle", "genie", or "nuwro"
    G4String fGenieFileName;   // GENIE ROOT file name (if using genie)
    G4String fNuWroFileName;   // NuWro ROOT file name (if using nuwro)
    G4int    fInitialEvent;

    G4bool fGeneratorChanged;  // Settings changed since the generator was built
    G4bool fRunStarted;        // A run has already been started
    // The generator this class installed itself, in sequential mode, and which
    // it therefore has to delete when it replaces it. The run manager owns the
    // one currently in use.
    G4VUserPrimaryGeneratorAction* fInstalledGenerator;
};

#endif
