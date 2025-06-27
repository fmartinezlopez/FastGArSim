#include "TrackingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

TrackingAction::TrackingAction(EventAction* eventAction,
                               RunAction* runAction)
: G4UserTrackingAction(),
  fEventAction(eventAction),
  fRunAction(runAction),
  fAnalysisManager(nullptr)
{
    // Get analysis manager
    fAnalysisManager = AnalysisManager::GetInstance();
}

TrackingAction::~TrackingAction()
{
    // No pointers to delete
}

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
    if (fRunAction->GetSaveOutput()) {
        fAnalysisManager->RecordTrackInfo(track);
    }
}

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
    if (fRunAction->GetSaveOutput()) {
        if (!fAnalysisManager->GetWriteTrajectory()) {
            fAnalysisManager->RecordTrackInfo(track, true);
        }
    }
}