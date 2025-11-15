#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction(DetectorConstruction* detector)
: G4UserRunAction(),
  fSaveOutput(true), fOutputFileName("gar_simulation"),
  fDetectorConstruction(detector),
  fAnalysisManager(nullptr)
{

    // Define command interface
    DefineCommands();

    // Get analysis manager
    fAnalysisManager = AnalysisManager::GetInstance();

    // Enable file merging
    G4AnalysisManager::Instance()->SetNtupleMerging(true);
    
    // For ROOT files, set default file type explicitly
    G4AnalysisManager::Instance()->SetDefaultFileType("root");
    
    // Set custom parameters if needed
    fAnalysisManager->SetOutputFileName(fOutputFileName);
}

RunAction::~RunAction()
{
    delete fMessenger;
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    
    // Inform the user
    G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

    // Update the filename in the analysis manager
    fAnalysisManager->SetOutputFileName(fOutputFileName);
    
    // Initialize analysis
    if (fSaveOutput) {
        fAnalysisManager->Book();

        // Record geometry AFTER Book() is called
        if (fDetectorConstruction) {
            fDetectorConstruction->RecordGeometry();
        }
    }

}

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;
    
    // Print run summary
    G4cout << "\n--------------------End of Run " << run->GetRunID() << " summary--------------------\n"
           << " Number of events processed: " << nofEvents << "\n"
           << "------------------------------------------------------------" << G4endl;
    
    // Finalize analysis
    if (fSaveOutput) {
        fAnalysisManager->Save();
        fAnalysisManager->Close();
    }
}

void RunAction::SetSaveOutput(G4bool save)
{
    fSaveOutput = save;
    G4cout << "Save output set to " << fSaveOutput << G4endl;
}

void RunAction::SetOutputFileName(G4String& name)
{
    fOutputFileName = name;
    G4cout << "Out name set to " << fOutputFileName << G4endl;
}

void RunAction::DefineCommands()
{
    // Initialize messenger for macro commands
    fMessenger = new G4GenericMessenger(this, "/run/", "Run control commands");

    fMessenger->DeclareProperty("OutputFileName", fOutputFileName,
                                "Set output file name");
    fMessenger->DeclareProperty("SaveOutput", fSaveOutput,
                                "Save output file");

}