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
  fECalHGAbsorberThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalHGAbsorberThickness", this);
  fECalHGAbsorberThicknessCmd->SetGuidance("Set ECal HG absorber thickness");
  fECalHGAbsorberThicknessCmd->SetParameterName("Thickness", false);
  fECalHGAbsorberThicknessCmd->SetUnitCategory("Length");
  fECalHGAbsorberThicknessCmd->SetRange("Thickness>0.");
  fECalHGAbsorberThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalHGScintillatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalHGScintillatorThickness", this);
  fECalHGScintillatorThicknessCmd->SetGuidance("Set ECal HG scintillator thickness");
  fECalHGScintillatorThicknessCmd->SetParameterName("Thickness", false);
  fECalHGScintillatorThicknessCmd->SetUnitCategory("Length");
  fECalHGScintillatorThicknessCmd->SetRange("Thickness>0.");
  fECalHGScintillatorThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalHGBoardThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalHGBoardThickness", this);
  fECalHGBoardThicknessCmd->SetGuidance("Set ECal HG PCB thickness");
  fECalHGBoardThicknessCmd->SetParameterName("Thickness", false);
  fECalHGBoardThicknessCmd->SetUnitCategory("Length");
  fECalHGBoardThicknessCmd->SetRange("Thickness>0.");
  fECalHGBoardThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalLGAbsorberThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalLGAbsorberThickness", this);
  fECalLGAbsorberThicknessCmd->SetGuidance("Set ECal LG absorber thickness");
  fECalLGAbsorberThicknessCmd->SetParameterName("Thickness", false);
  fECalLGAbsorberThicknessCmd->SetUnitCategory("Length");
  fECalLGAbsorberThicknessCmd->SetRange("Thickness>0.");
  fECalLGAbsorberThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalLGScintillatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/ECalLGScintillatorThickness", this);
  fECalLGScintillatorThicknessCmd->SetGuidance("Set ECal LG scintillator thickness");
  fECalLGScintillatorThicknessCmd->SetParameterName("Thickness", false);
  fECalLGScintillatorThicknessCmd->SetUnitCategory("Length");
  fECalLGScintillatorThicknessCmd->SetRange("Thickness>0.");
  fECalLGScintillatorThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalBarrelHGLayersCmd = new G4UIcmdWithAnInteger("/detector/ECalBarrelHGLayers", this);
  fECalBarrelHGLayersCmd->SetGuidance("Set number of HG layers in ECal barrel");
  fECalBarrelHGLayersCmd->SetParameterName("Layers", false);
  fECalBarrelHGLayersCmd->SetRange("Layers>0");
  fECalBarrelHGLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalBarrelLGLayersCmd = new G4UIcmdWithAnInteger("/detector/ECalBarrelLGLayers", this);
  fECalBarrelLGLayersCmd->SetGuidance("Set number of LG layers in ECal barrel");
  fECalBarrelLGLayersCmd->SetParameterName("Layers", false);
  fECalBarrelLGLayersCmd->SetRange("Layers>0");
  fECalBarrelLGLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalEndcapHGLayersCmd = new G4UIcmdWithAnInteger("/detector/ECalEndcapHGLayers", this);
  fECalEndcapHGLayersCmd->SetGuidance("Set number of HG layers in ECal end cap");
  fECalEndcapHGLayersCmd->SetParameterName("Layers", false);
  fECalEndcapHGLayersCmd->SetRange("Layers>0");
  fECalEndcapHGLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fECalEndcapLGLayersCmd = new G4UIcmdWithAnInteger("/detector/ECalEndcapLGLayers", this);
  fECalEndcapLGLayersCmd->SetGuidance("Set number of LG layers in ECal end cap");
  fECalEndcapLGLayersCmd->SetParameterName("Layers", false);
  fECalEndcapLGLayersCmd->SetRange("Layers>0");
  fECalEndcapLGLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

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

  // Material commands
  fPressureCmd = new G4UIcmdWithADoubleAndUnit("/detector/GasPressure", this);
  fPressureCmd->SetGuidance("Set gas pressure");
  fPressureCmd->SetParameterName("GasPressure", false);
  fPressureCmd->SetUnitCategory("Pressure");
  fPressureCmd->SetRange("GasPressure>0.");
  fPressureCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fGeometryCmd;
  delete fTPCRadiusCmd;
  delete fTPCLengthCmd;
  delete fECalHGAbsorberThicknessCmd;
  delete fECalHGScintillatorThicknessCmd;
  delete fECalHGBoardThicknessCmd;
  delete fECalLGAbsorberThicknessCmd;
  delete fECalLGScintillatorThicknessCmd;
  delete fECalBarrelHGLayersCmd;
  delete fECalBarrelLGLayersCmd;
  delete fECalEndcapHGLayersCmd;
  delete fECalEndcapLGLayersCmd;
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
  delete fPressureCmd;
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
  else if (command == fECalHGAbsorberThicknessCmd) {
    fDetector->SetECalHGAbsorberThickness(fECalHGAbsorberThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalHGScintillatorThicknessCmd) {
    fDetector->SetECalHGScintillatorThickness(fECalHGScintillatorThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalHGBoardThicknessCmd) {
    fDetector->SetECalHGBoardThickness(fECalHGBoardThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalLGAbsorberThicknessCmd) {
    fDetector->SetECalLGAbsorberThickness(fECalLGAbsorberThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalLGScintillatorThicknessCmd) {
    fDetector->SetECalLGScintillatorThickness(fECalLGScintillatorThicknessCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fECalBarrelHGLayersCmd) {
    fDetector->SetECalBarrelHGLayers(fECalBarrelHGLayersCmd->GetNewIntValue(newValue));
  }
  else if (command == fECalBarrelLGLayersCmd) {
    fDetector->SetECalBarrelLGLayers(fECalBarrelLGLayersCmd->GetNewIntValue(newValue));
  }
  else if (command == fECalEndcapHGLayersCmd) {
    fDetector->SetECalEndcapHGLayers(fECalEndcapHGLayersCmd->GetNewIntValue(newValue));
  }
  else if (command == fECalEndcapLGLayersCmd) {
    fDetector->SetECalEndcapLGLayers(fECalEndcapLGLayersCmd->GetNewIntValue(newValue));
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
  else if (command == fPressureCmd) {
    fDetector->SetPressure(fPressureCmd->GetNewDoubleValue(newValue));
  }
}
