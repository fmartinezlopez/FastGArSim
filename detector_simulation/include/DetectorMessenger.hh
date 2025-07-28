#ifndef DETECTORMESSENGER_HH
#define DETECTORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction*);
    virtual ~DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DetectorConstruction* fDetector;
    
    G4UIdirectory*             fDetectorDir;
    G4UIcmdWithAString*        fGeometryCmd;
    G4UIcmdWithADoubleAndUnit* fTPCRadiusCmd;
    G4UIcmdWithADoubleAndUnit* fTPCLengthCmd;
    G4UIcmdWithADoubleAndUnit* fECalAbsorberThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fECalScintillatorThicknessCmd;
    G4UIcmdWithAnInteger*      fECalLayersCmd;
    G4UIcmdWithADoubleAndUnit* fMuIDAbsorberThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fMuIDScintillatorThicknessCmd;
    G4UIcmdWithAnInteger*      fMuIDLayersCmd;
    G4UIcmdWithAnInteger*      fLArNModulesXCmd;
    G4UIcmdWithAnInteger*      fLArNModulesYCmd;
    G4UIcmdWithAnInteger*      fLArNModulesZCmd;
    G4UIcmdWithADoubleAndUnit* fLArModuleLengthCmd;
    G4UIcmdWithADoubleAndUnit* fLArModuleWidthCmd;
    G4UIcmdWithADoubleAndUnit* fLArModuleDepthCmd;
    G4UIcmdWithADoubleAndUnit* fLArModuleGapCmd;
    G4UIcmdWithADoubleAndUnit* fLArInsulationThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fLArCryostatThicknessCmd;
    G4UIcmdWithABool*          fLArEnableMuonWindowCmd;
    G4UIcmdWithADoubleAndUnit* fLArMuonWindowThicknessCmd;
};

#endif