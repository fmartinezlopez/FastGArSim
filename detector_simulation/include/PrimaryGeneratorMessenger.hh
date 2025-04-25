#ifndef PRIMARYGENERATORMESSENGER_HH
#define PRIMARYGENERATORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class ActionInitialization;
class G4UIdirectory;
class G4UIcmdWithAString;

class PrimaryGeneratorMessenger : public G4UImessenger
{
public:
  PrimaryGeneratorMessenger(ActionInitialization* actionInit);
  virtual ~PrimaryGeneratorMessenger();
  
  virtual void SetNewValue(G4UIcommand*, G4String);
  
private:
  ActionInitialization* fActionInitialization;
  
  G4UIdirectory*        fDirectory;
  G4UIcmdWithAString*   fGeneratorTypeCmd;
  G4UIcmdWithAString*   fGenieFileNameCmd;
  G4UIcmdWithAString*   fNuWroFileNameCmd;
};

#endif