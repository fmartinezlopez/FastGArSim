#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "GENIEGeneratorAction.hh"
#include "NuWroGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
: G4VUserActionInitialization(),
  fDetConstruction(detConstruction),
  fGeneratorType("particle"),
  fGenieFileName(""),
  fNuWroFileName(""),
  fInitialEvent(0)
{
    // Nothing else to initialize
}

ActionInitialization::~ActionInitialization()
{
    // No pointers to delete
}

void ActionInitialization::BuildForMaster() const
{
    // Create and register run action for master thread
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);
}

void ActionInitialization::Build() const
{

    // Create primary generator action based on the selected type
    G4VUserPrimaryGeneratorAction* primaryGeneratorAction = nullptr;

    if (fGeneratorType == "genie") {
        // Check if GENIE file name is set
        if (fGenieFileName.empty()) {
            G4cerr << "Error: GENIE generator selected but no input file specified." << G4endl;
            G4cerr << "Falling back to default particle gun generator." << G4endl;
            primaryGeneratorAction = new PrimaryGeneratorAction();
        } else {
            G4cout << "Using GENIE primary generator with file: " << fGenieFileName << G4endl;
            primaryGeneratorAction = new GENIEGeneratorAction(fGenieFileName, fInitialEvent);
        }
    } else if (fGeneratorType == "nuwro") {
        // Check if NuWro file name is set
        if (fNuWroFileName.empty()) {
            G4cerr << "Error: NuWro generator selected but no input file specified." << G4endl;
            G4cerr << "Falling back to default particle gun generator." << G4endl;
            primaryGeneratorAction = new PrimaryGeneratorAction();
        } else {
            G4cout << "Using NuWro primary generator with file: " << fNuWroFileName << G4endl;
            primaryGeneratorAction = new NuWroGeneratorAction(fNuWroFileName, fInitialEvent);
        }
    } else {
        // Default or explicitly specified "particle" type
        G4cout << "Using standard particle gun generator" << G4endl;
        primaryGeneratorAction = new PrimaryGeneratorAction();
    }

    // Register the primary generator action
    SetUserAction(primaryGeneratorAction);
    
    // Create and register run action
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);
    
    // Create and register event action
    EventAction* eventAction = new EventAction(runAction);
    SetUserAction(eventAction);

    // Create and register tracking action
    TrackingAction* trackingAction = new TrackingAction(eventAction, runAction);
    SetUserAction(trackingAction);
    
    // Create and register stepping action
    SteppingAction* steppingAction = new SteppingAction(eventAction, runAction, fDetConstruction);
    SetUserAction(steppingAction);
}

void ActionInitialization::UpdatePrimaryGeneratorAction() const
{
    // Create new generator based on current settings
    G4VUserPrimaryGeneratorAction* primaryGeneratorAction = nullptr;

    if (fGeneratorType == "genie") {
        // Check if GENIE file name is set
        if (fGenieFileName.empty()) {
            G4cerr << "Error: GENIE generator selected but no input file specified." << G4endl;
            G4cerr << "Falling back to default particle gun generator." << G4endl;
            primaryGeneratorAction = new PrimaryGeneratorAction();
        } else {
            G4cout << "Using GENIE primary generator with file: " << fGenieFileName << G4endl;
            primaryGeneratorAction = new GENIEGeneratorAction(fGenieFileName, fInitialEvent);
        }
    } else if (fGeneratorType == "nuwro") {
        // Check if NuWro file name is set
        if (fNuWroFileName.empty()) {
            G4cerr << "Error: NuWro generator selected but no input file specified." << G4endl;
            G4cerr << "Falling back to default particle gun generator." << G4endl;
            primaryGeneratorAction = new PrimaryGeneratorAction();
        } else {
            G4cout << "Using NuWro primary generator with file: " << fNuWroFileName << G4endl;
            primaryGeneratorAction = new NuWroGeneratorAction(fNuWroFileName, fInitialEvent);
        }
    } else {
        // Default or explicitly specified "particle" type
        G4cout << "Using standard particle gun generator" << G4endl;
        primaryGeneratorAction = new PrimaryGeneratorAction();
    }
    
    // Update the primary generator action
    SetUserAction(primaryGeneratorAction);
}