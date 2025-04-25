#include "PrimaryGeneratorMessenger.hh"
#include "ActionInitialization.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(ActionInitialization* actionInit)
: G4UImessenger(),
  fActionInitialization(actionInit)
{
  // Create directories for commands
  fDirectory = new G4UIdirectory("/generator/");
  fDirectory->SetGuidance("Primary generator control commands.");
  
  // Command to select generator type
  fGeneratorTypeCmd = new G4UIcmdWithAString("/generator/select", this);
  fGeneratorTypeCmd->SetGuidance("Select the primary generator type.");
  fGeneratorTypeCmd->SetGuidance("  Choice : particle, genie, nuwro");
  fGeneratorTypeCmd->SetParameterName("Type", false);
  fGeneratorTypeCmd->SetCandidates("particle genie nuwro");
  fGeneratorTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  // Command to set GENIE input file name
  fGenieFileNameCmd = new G4UIcmdWithAString("/generator/genieFile", this);
  fGenieFileNameCmd->SetGuidance("Set the GENIE input ROOT file name.");
  fGenieFileNameCmd->SetParameterName("FileName", false);
  fGenieFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set NuWro input file name
  fNuWroFileNameCmd = new G4UIcmdWithAString("/generator/nuwroFile", this);
  fNuWroFileNameCmd->SetGuidance("Set the NuWro input ROOT file name.");
  fNuWroFileNameCmd->SetParameterName("FileName", false);
  fNuWroFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fGeneratorTypeCmd;
  delete fGenieFileNameCmd;
  delete fNuWroFileNameCmd;
  delete fDirectory;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fGeneratorTypeCmd) {
    fActionInitialization->SetGeneratorType(newValue);
    fActionInitialization->UpdatePrimaryGeneratorAction();
  }
  else if (command == fGenieFileNameCmd) {
    fActionInitialization->SetGenieFileName(newValue);
    fActionInitialization->UpdatePrimaryGeneratorAction();
  }
  else if (command == fNuWroFileNameCmd) {
    fActionInitialization->SetNuWroFileName(newValue);
    fActionInitialization->UpdatePrimaryGeneratorAction();
  }
}