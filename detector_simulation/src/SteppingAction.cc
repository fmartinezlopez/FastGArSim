#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4TrackingManager.hh"
#include "G4TrajectoryPoint.hh"

SteppingAction::SteppingAction(EventAction* eventAction,
                               RunAction* runAction,
                               DetectorConstruction* detConstruction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(runAction),
  fDetConstruction(detConstruction),
  fAnalysisManager(nullptr),
  fTPCLogical(nullptr),
  fECalLogical(nullptr),
  fMuIDLogical(nullptr)
{

    // Get analysis manager
    fAnalysisManager = AnalysisManager::GetInstance();

    // Get "sensitive" logical volumes for identification
    fTPCLogical  = fDetConstruction->GetTPCSensLogical();
    fECalLogical = fDetConstruction->GetECalSensLogical();
    fMuIDLogical = fDetConstruction->GetMuIDSensLogical();

}

SteppingAction::~SteppingAction()
{
    // No pointers to delete
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (fRunAction->GetSaveOutput()) {
        // Record energy deposit and check for TPC/ECal hits
        fAnalysisManager->RecordEnergyDeposit(step);
        fAnalysisManager->RecordTrackInfo(step->GetTrack());
    }
}