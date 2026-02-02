#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"
#include "G4GDMLParser.hh"

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"

// Print usage instructions
void PrintUsage()
{
    G4cerr << "Usage: GArSimulation [options] [macro_file] [arguments]\n"
           << "Options:\n"
           << " -m MACRO     Execute macro_file with UI session\n"
           << " -u UISESSION Launch specific UI session (e.g., Qt, Xm, tcsh)\n"
           << " -v           Enable visualization\n"
           << " -g           Save current geometry as GDML file\n"
           << " -h           Print this help message\n"
           << G4endl;
}

int main(int argc, char** argv)
{
    // Parse command line arguments
    G4String macro;
    G4String session;
    G4bool visualization = false;
    G4bool geometry = false;
    
    for (G4int i = 1; i < argc; i++)
    {
        G4String arg = argv[i];
        if (arg == "-m" && i+1 < argc)
        {
            macro = argv[++i];
        }
        else if (arg == "-u" && i+1 < argc)
        {
            session = argv[++i];
        }
        else if (arg == "-v")
        {
            visualization = true;
        }
        else if (arg == "-g")
        {
            geometry = true;
        }
        else if (arg == "-h")
        {
            PrintUsage();
            return 0;
        }
        else if (!macro.size())
        {
            // First unknown argument is assumed to be macro file
            macro = arg;
        }
        else
        {
            G4cerr << "Unknown option: " << arg << G4endl;
            PrintUsage();
            return 1;
        }
    }
    
    // Initialize random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4Random::setTheSeed(time(NULL));

    // Construct the run manager
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    
    // Set mandatory initialization classes
    // Detector construction
    auto detConstruction = new DetectorConstruction();
    runManager->SetUserInitialization(detConstruction);
    
    // Physics list
    auto physicsList = new PhysicsList();
    runManager->SetUserInitialization(physicsList);
    
    // Action initialization
    auto actionInitialization = new ActionInitialization(detConstruction);
    runManager->SetUserInitialization(actionInitialization);

    // Create the messenger to allow macro control of generator type
    PrimaryGeneratorMessenger* messenger = new PrimaryGeneratorMessenger(actionInitialization);
    
    // Initialize visualization
    G4VisManager* visManager = nullptr;
    if (visualization)
    {
        // Visualization manager
        visManager = new G4VisExecutive();
        visManager->Initialize();
    }
    
    // Get the UI manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if (macro.size())
    {
        // Batch mode - execute the specified macro
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macro);
    }
    else
    {
        // Interactive mode
        G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
        if (visualization)
        {
            // Execute default visualization macro
            UImanager->ApplyCommand("/control/execute macros/vis.mac");
        }
        ui->SessionStart();
        delete ui;
    }
    
    // Export to GDML after initialization
    if (geometry)
    {
        G4GDMLParser parser;
        const G4String filename = "GArGeometry.gdml";
        std::remove(filename.c_str());  // delete if exists
        parser.Write(filename, detConstruction->GetWorldVolume());
    }

    // Cleanup
    if (visManager) delete visManager;
    delete runManager;
    delete messenger;
    
    return 0;
}
