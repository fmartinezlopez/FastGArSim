#include "EventAction.hh"
#include "RunAction.hh"
#include "AnalysisManager.hh"

#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fAnalysisManager(nullptr)
{
    // Get analysis manager
    fAnalysisManager = AnalysisManager::GetInstance();
}

EventAction::~EventAction()
{
    // No pointers to delete
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    // Start a new event
    fAnalysisManager->BeginEvent(event->GetEventID());
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Finish the event and write data
    fAnalysisManager->EndEvent();
}