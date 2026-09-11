#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "GunGeneratorAction.hh"
#include "GENIEGeneratorAction.hh"
#include "NuWroGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"

#include "G4RunManager.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector)
: G4VUserActionInitialization(),
  G4VStateDependent(),
  fDetectorConstruction(detector),
  fGeneratorType("particle"),
  fGenieFileName(""),
  fNuWroFileName(""),
  fInitialEvent(0),
  fGeneratorChanged(false),
  fRunStarted(false),
  fInstalledGenerator(nullptr)
{
    // Nothing else to initialize
}

ActionInitialization::~ActionInitialization()
{
    // fInstalledGenerator is either already deleted or owned by the run manager
}

G4VUserPrimaryGeneratorAction* ActionInitialization::CreatePrimaryGeneratorAction() const
{
    if (fGeneratorType == "genie") {
        // Check if GENIE file name is set
        if (fGenieFileName.empty()) {
            G4cerr << "Error: GENIE generator selected but no input file specified." << G4endl;
            G4cerr << "Falling back to default particle gun generator." << G4endl;
            return new GunGeneratorAction();
        }
        G4cout << "Using GENIE primary generator with file: " << fGenieFileName << G4endl;
        return new GENIEGeneratorAction(fGenieFileName, fInitialEvent);
    }

    if (fGeneratorType == "nuwro") {
        // Check if NuWro file name is set
        if (fNuWroFileName.empty()) {
            G4cerr << "Error: NuWro generator selected but no input file specified." << G4endl;
            G4cerr << "Falling back to default particle gun generator." << G4endl;
            return new GunGeneratorAction();
        }
        G4cout << "Using NuWro primary generator with file: " << fNuWroFileName << G4endl;
        return new NuWroGeneratorAction(fNuWroFileName, fInitialEvent);
    }

    // Default or explicitly specified "particle" type
    G4cout << "Using standard particle gun generator" << G4endl;
    return new GunGeneratorAction();
}

void ActionInitialization::BuildForMaster() const
{
    // Create and register run action for master thread
    RunAction* runAction = new RunAction(fDetectorConstruction);
    SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
    // Register the primary generator action. On a worker thread this runs at
    // the start of a run and the settings are final; in sequential mode it
    // runs before the macro has been read and Notify() rebuilds it later.
    SetUserAction(CreatePrimaryGeneratorAction());
    
    // Create and register run action
    RunAction* runAction = new RunAction(fDetectorConstruction);
    SetUserAction(runAction);
    
    // Create and register event action
    EventAction* eventAction = new EventAction(runAction);
    SetUserAction(eventAction);

    // Create and register tracking action
    TrackingAction* trackingAction = new TrackingAction(eventAction, runAction);
    SetUserAction(trackingAction);
    
    // Create and register stepping action
    SteppingAction* steppingAction = new SteppingAction(eventAction, runAction);
    SetUserAction(steppingAction);
}

G4bool ActionInitialization::Notify(G4ApplicationState requestedState)
{
    // The geometry is closed at the start of every run, so this is the last
    // point at which the macro can still have changed the configuration
    if (requestedState != G4State_GeomClosed) return true;

    const G4bool firstRun = !fRunStarted;
    fRunStarted = true;

    if (!fGeneratorChanged) return true;

    // Only the sequential run manager keeps the generator on this thread. In
    // multi-threaded mode every worker builds its own in Build(), and calling
    // SetUserAction() here would abort the job.
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (runManager == nullptr ||
        runManager->GetRunManagerType() != G4RunManager::sequentialRM) {
        // The workers build their generator when they start, so the settings
        // are picked up for the first run but never again afterwards
        if (!firstRun) {
            G4cerr << "Warning: the primary generator cannot be changed after the first run"
                   << " in multi-threaded mode. The change will be ignored." << G4endl;
        }
        fGeneratorChanged = false;
        return true;
    }

    G4VUserPrimaryGeneratorAction* generator = CreatePrimaryGeneratorAction();
    runManager->SetUserAction(generator);

    // The run manager has dropped the generator installed by the previous run
    delete fInstalledGenerator;
    fInstalledGenerator = generator;
    fGeneratorChanged = false;

    return true;
}
