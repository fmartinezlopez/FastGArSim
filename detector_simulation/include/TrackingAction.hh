#ifndef TRACKINGACTION_HH
#define TRACKINGACTION_HH

#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class EventAction;
class RunAction;
class AnalysisManager;
class G4Track;

class TrackingAction : public G4UserTrackingAction
{
  public:
    TrackingAction(EventAction* eventAction, RunAction* runAction);
    virtual ~TrackingAction();
  
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
    
  private:
    RunAction*    fRunAction;
    EventAction*  fEventAction;
    AnalysisManager* fAnalysisManager;
};

#endif