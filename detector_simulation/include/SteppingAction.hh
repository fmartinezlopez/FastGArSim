#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class EventAction;
class RunAction;
class DetectorConstruction;
class AnalysisManager;

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(EventAction* eventAction, RunAction* runAction, DetectorConstruction* detConstruction);
    virtual ~SteppingAction();
    
    // Method from base class
    virtual void UserSteppingAction(const G4Step* step);
    
private:
    EventAction* fEventAction;
    RunAction* fRunAction;
    DetectorConstruction* fDetConstruction;
    AnalysisManager* fAnalysisManager;
    
    // Logical volumes for identification
    G4LogicalVolume* fTPCLogical;
    G4LogicalVolume* fECalLogical;
    G4LogicalVolume* fMuIDLogical;
};

#endif