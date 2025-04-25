#ifndef DETECTORMESSENGER_HH
#define DETECTORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction*);
    virtual ~DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DetectorConstruction* fDetector;
    
    G4UIdirectory*             fDetectorDir;
    G4UIcmdWithADoubleAndUnit* fTPCRadiusCmd;
    G4UIcmdWithADoubleAndUnit* fTPCLengthCmd;
    G4UIcmdWithADoubleAndUnit* fECalAbsorberThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fECalScintillatorThicknessCmd;
    G4UIcmdWithAnInteger*      fECalLayersCmd;
    G4UIcmdWithADoubleAndUnit* fMuIDAbsorberThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fMuIDScintillatorThicknessCmd;
    G4UIcmdWithAnInteger*      fMuIDLayersCmd;
};

#endif