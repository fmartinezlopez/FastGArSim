#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
 : G4UImessenger(),
   fDetector(det)
{
  fDetectorDir = new G4UIdirectory("/detector/");
  fDetectorDir->SetGuidance("Detector construction control");

  // Geometry selection command
  fGeometryCmd = new G4UIcmdWithAString("/detector/setGeometry", this);
  fGeometryCmd->SetGuidance("Select detector geometry type");
  fGeometryCmd->SetGuidance("  gar : GAr TPC + ECal + MuID");
  fGeometryCmd->SetGuidance("  lar : LAr TPC");
  fGeometryCmd->SetParameterName("GeometryType", false);
  fGeometryCmd->SetCandidates("gar lar");
  fGeometryCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
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

  // LArTPC commands
  fLArNModulesXCmd = new G4UIcmdWithAnInteger("/detector/LArNModulesX", this);
  fLArNModulesXCmd->SetGuidance("Set number of LAr modules in X direction");
  fLArNModulesXCmd->SetParameterName("LArNModuleX", false);
  fLArNModulesXCmd->SetRange("LArNModuleX>0");
  fLArNModulesXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArNModulesYCmd = new G4UIcmdWithAnInteger("/detector/LArNModulesY", this);
  fLArNModulesYCmd->SetGuidance("Set number of LAr modules in Y direction");
  fLArNModulesYCmd->SetParameterName("LArNModuleY", false);
  fLArNModulesYCmd->SetRange("LArNModuleY>0");
  fLArNModulesYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArNModulesZCmd = new G4UIcmdWithAnInteger("/detector/LArNModulesZ", this);
  fLArNModulesZCmd->SetGuidance("Set number of LAr modules in Z direction");
  fLArNModulesZCmd->SetParameterName("LArNModuleZ", false);
  fLArNModulesZCmd->SetRange("LArNModuleZ>0");
  fLArNModulesZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArModuleLengthCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArModuleLength", this);
  fLArModuleLengthCmd->SetGuidance("Set LAr Module length");
  fLArModuleLengthCmd->SetParameterName("LArModuleLength", false);
  fLArModuleLengthCmd->SetUnitCategory("Length");
  fLArModuleLengthCmd->SetRange("LArModuleLength>0.");
  fLArModuleLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArModuleWidthCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArModuleWidth", this);
  fLArModuleWidthCmd->SetGuidance("Set LAr Module width");
  fLArModuleWidthCmd->SetParameterName("LArModuleWidth", false);
  fLArModuleWidthCmd->SetUnitCategory("Length");
  fLArModuleWidthCmd->SetRange("LArModuleWidth>0.");
  fLArModuleWidthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArModuleDepthCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArModuleDepth", this);
  fLArModuleDepthCmd->SetGuidance("Set LAr Module depth");
  fLArModuleDepthCmd->SetParameterName("LArModuleDepth", false);
  fLArModuleDepthCmd->SetUnitCategory("Length");
  fLArModuleDepthCmd->SetRange("LArModuleDepth>0.");
  fLArModuleDepthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArModuleGapCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArModuleGap", this);
  fLArModuleGapCmd->SetGuidance("Set LAr Module Gap");
  fLArModuleGapCmd->SetParameterName("LArModuleGap", false);
  fLArModuleGapCmd->SetUnitCategory("Length");
  fLArModuleGapCmd->SetRange("LArModuleGap>0.");
  fLArModuleGapCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArInsulationThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArInsulationThickness", this);
  fLArInsulationThicknessCmd->SetGuidance("Set LAr Module Insulation Thickness");
  fLArInsulationThicknessCmd->SetParameterName("LArInsulationThickness", false);
  fLArInsulationThicknessCmd->SetUnitCategory("Length");
  fLArInsulationThicknessCmd->SetRange("LArInsulationThickness>0.");
  fLArInsulationThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArCryostatThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArCryostatThickness", this);
  fLArCryostatThicknessCmd->SetGuidance("Set LAr Cryostat Thickness");
  fLArCryostatThicknessCmd->SetParameterName("LArCryostatThickness", false);
  fLArCryostatThicknessCmd->SetUnitCategory("Length");
  fLArCryostatThicknessCmd->SetRange("LArCryostatThickness>0.");
  fLArCryostatThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArEnableMuonWindowCmd = new G4UIcmdWithABool("/detector/LArEnableMuonWindow", this);
  fLArEnableMuonWindowCmd->SetGuidance("Enable LAr Muon Window");
  fLArEnableMuonWindowCmd->SetParameterName("LArEnableMuonWindow", false);
  fLArEnableMuonWindowCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLArMuonWindowThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/LArMuonWindowThickness", this);
  fLArMuonWindowThicknessCmd->SetGuidance("Set LAr Muon Window Thickness");
  fLArMuonWindowThicknessCmd->SetParameterName("LArMuonWindowThickness", false);
  fLArMuonWindowThicknessCmd->SetUnitCategory("Length");
  fLArMuonWindowThicknessCmd->SetRange("LArMuonWindowThickness>0.");
  fLArMuonWindowThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fGeometryCmd;
  delete fTPCRadiusCmd;
  delete fTPCLengthCmd;
  delete fECalAbsorberThicknessCmd;
  delete fECalScintillatorThicknessCmd;
  delete fECalLayersCmd;
  delete fMuIDAbsorberThicknessCmd;
  delete fMuIDScintillatorThicknessCmd;
  delete fMuIDLayersCmd;
  delete fLArNModulesXCmd;
  delete fLArNModulesYCmd;
  delete fLArNModulesZCmd;
  delete fLArModuleLengthCmd;
  delete fLArModuleWidthCmd;
  delete fLArModuleDepthCmd;
  delete fLArModuleGapCmd;
  delete fLArInsulationThicknessCmd;
  delete fLArCryostatThicknessCmd;
  delete fLArEnableMuonWindowCmd;
  delete fLArMuonWindowThicknessCmd;
  delete fDetectorDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fGeometryCmd) {
    fDetector->SetGeometryType(newValue);
  }
  else if (command == fTPCRadiusCmd) {
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
  else if (command == fLArNModulesXCmd) {
    fDetector->SetLArNModulesX(fLArNModulesXCmd->GetNewIntValue(newValue));
  }
  else if (command == fLArNModulesYCmd) {
    fDetector->SetLArNModulesY(fLArNModulesYCmd->GetNewIntValue(newValue));
  }
  else if (command == fLArNModulesZCmd) {
    fDetector->SetLArNModulesZ(fLArNModulesZCmd->GetNewIntValue(newValue));
  }
  else if (command == fLArModuleLengthCmd) {
    fDetector->SetLArModuleLength(fLArModuleLengthCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLArModuleWidthCmd) {
    fDetector->SetLArModuleWidth(fLArModuleWidthCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLArModuleDepthCmd) {
    fDetector->SetLArModuleDepth(fLArModuleDepthCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLArModuleGapCmd) {
    fDetector->SetLArModuleGap(fLArModuleGapCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLArInsulationThicknessCmd) {
    fDetector->SetLArInsulationThickness(fLArInsulationThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLArCryostatThicknessCmd) {
    fDetector->SetLArCryostatThickness(fLArCryostatThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLArEnableMuonWindowCmd) {
    fDetector->SetLArEnableMuonWindow(fLArEnableMuonWindowCmd->GetNewBoolValue(newValue));
  }
  else if (command == fLArMuonWindowThicknessCmd) {
    fDetector->SetLArMuonWindowThickness(fLArMuonWindowThicknessCmd->GetNewDoubleValue(newValue));
  }
}