#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
 : G4UImessenger(),
   fDetector(det)
{
  fDetectorDir = new G4UIdirectory("/detector/");
  fDetectorDir->SetGuidance("Detector construction control");
  
  // TPC commands
  fTPCRadiusCmd = new G4UIcmdWithADoubleAndUnit("/detector/TPCRadius", this);
  fTPCRadiusCmd->SetGuidance("Set TPC radius");
  fTPCRadiusCmd->SetParameterName("Radius", false);
  fTPCRadiusCmd->SetUnitCategory("Length");
  fTPCRadiusCmd->SetRange("Radius>0.");
  fTPCRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTPCLengthCmd = new G4UIcmdWithADoubleAndUnit("/detector/TPCLength", this);
  fTPCLengthCmd->SetGuidance("Set TPC length");
  fTPCLengthCmd->SetParameterName("Length", false);
  fTPCLengthCmd->SetUnitCategory("Length");
  fTPCLengthCmd->SetRange("Length>0.");
  fTPCLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // ECal commands
  fECalAbsorberThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalAbsorberThickness", this);
  fECalAbsorberThicknessCmd->SetGuidance("Set ECal absorber thickness");
  fECalAbsorberThicknessCmd->SetParameterName("Thickness", false);
  fECalAbsorberThicknessCmd->SetUnitCategory("Length");
  fECalAbsorberThicknessCmd->SetRange("Thickness>0.");
  fECalAbsorberThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalScintillatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalScintillatorThickness", this);
  fECalScintillatorThicknessCmd->SetGuidance("Set ECal scintillator thickness");
  fECalScintillatorThicknessCmd->SetParameterName("Thickness", false);
  fECalScintillatorThicknessCmd->SetUnitCategory("Length");
  fECalScintillatorThicknessCmd->SetRange("Thickness>0.");
  fECalScintillatorThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalLayersCmd = new G4UIcmdWithAnInteger("/detector/ECalLayers", this);
  fECalLayersCmd->SetGuidance("Set number of layers in ECal");
  fECalLayersCmd->SetParameterName("Layers", false);
  fECalLayersCmd->SetRange("Layers>0");
  fECalLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // MuID commands
  fMuIDAbsorberThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/MuIDAbsorberThickness", this);
  fMuIDAbsorberThicknessCmd->SetGuidance("Set MuID absorber thickness");
  fMuIDAbsorberThicknessCmd->SetParameterName("Thickness", false);
  fMuIDAbsorberThicknessCmd->SetUnitCategory("Length");
  fMuIDAbsorberThicknessCmd->SetRange("Thickness>0.");
  fMuIDAbsorberThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMuIDScintillatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/MuIDScintillatorThickness", this);
  fMuIDScintillatorThicknessCmd->SetGuidance("Set MuID scintillator thickness");
  fMuIDScintillatorThicknessCmd->SetParameterName("Thickness", false);
  fMuIDScintillatorThicknessCmd->SetUnitCategory("Length");
  fMuIDScintillatorThicknessCmd->SetRange("Thickness>0.");
  fMuIDScintillatorThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMuIDLayersCmd = new G4UIcmdWithAnInteger("/detector/MuIDLayers", this);
  fMuIDLayersCmd->SetGuidance("Set number of layers in MuID");
  fMuIDLayersCmd->SetParameterName("Layers", false);
  fMuIDLayersCmd->SetRange("Layers>0");
  fMuIDLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fTPCRadiusCmd;
  delete fTPCLengthCmd;
  delete fECalAbsorberThicknessCmd;
  delete fECalScintillatorThicknessCmd;
  delete fECalLayersCmd;
  delete fMuIDAbsorberThicknessCmd;
  delete fMuIDScintillatorThicknessCmd;
  delete fMuIDLayersCmd;
  delete fDetectorDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fTPCRadiusCmd) {
    fDetector->SetTPCRadius(fTPCRadiusCmd->GetNewDoubleValue(newValue));
  } 
  else if (command == fTPCLengthCmd) {
    fDetector->SetTPCLength(fTPCLengthCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalAbsorberThicknessCmd) {
    fDetector->SetECalAbsorberThickness(fECalAbsorberThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalScintillatorThicknessCmd) {
    fDetector->SetECalScintillatorThickness(fECalScintillatorThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalLayersCmd) {
    fDetector->SetECalLayers(fECalLayersCmd->GetNewIntValue(newValue));
  }
  else if (command == fMuIDAbsorberThicknessCmd) {
    fDetector->SetMuIDAbsorberThickness(fMuIDAbsorberThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fMuIDScintillatorThicknessCmd) {
    fDetector->SetMuIDScintillatorThickness(fMuIDScintillatorThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fMuIDLayersCmd) {
    fDetector->SetMuIDLayers(fMuIDLayersCmd->GetNewIntValue(newValue));
  }
}