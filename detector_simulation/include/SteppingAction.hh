#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class EventAction;
class RunAction;
class AnalysisManager;

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(EventAction* eventAction, RunAction* runAction);
    virtual ~SteppingAction();
    
    // Method from base class
    virtual void UserSteppingAction(const G4Step* step);
    
private:
    EventAction* fEventAction;
    RunAction* fRunAction;
    AnalysisManager* fAnalysisManager;

};

#endif