#ifndef GUNGENERATORMESSENGER_HH
#define GUNGENERATORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class GunGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;

class GunGeneratorMessenger : public G4UImessenger
{
public:
  GunGeneratorMessenger(GunGeneratorAction* gunGen);
  virtual ~GunGeneratorMessenger();
  
  virtual void SetNewValue(G4UIcommand*, G4String);
  
private:
  GunGeneratorAction* fGunGeneratorAction;
  
  G4UIdirectory*             fDirectory;
  G4UIcmdWithAString*        fParticleTypeCmd;
  G4UIcmdWithADoubleAndUnit* fMomentumCmd;
  G4UIcmdWithADoubleAndUnit* fMomentumSpreadCmd;
  G4UIcmdWithAString*        fMomentumDistCmd;
  G4UIcmdWith3VectorAndUnit* fPositionCmd;
  G4UIcmdWith3VectorAndUnit* fPositionSpreadCmd;
  G4UIcmdWithAString*        fPositionDistCmd;
  G4UIcmdWithADoubleAndUnit* fPositionRMaxCmd;
  G4UIcmdWithADoubleAndUnit* fXZAngleCmd;
  G4UIcmdWithADoubleAndUnit* fXZAngleSpreadCmd;
  G4UIcmdWithAString*        fXZAngleDistCmd;
  G4UIcmdWithADoubleAndUnit* fXYAngleCmd;
  G4UIcmdWithADoubleAndUnit* fXYAngleSpreadCmd;
  G4UIcmdWithAString*        fXYAngleDistCmd;

};

#endif
