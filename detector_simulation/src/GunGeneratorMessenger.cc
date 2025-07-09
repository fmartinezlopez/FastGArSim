#include "GunGeneratorMessenger.hh"
#include "GunGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

GunGeneratorMessenger::GunGeneratorMessenger(GunGeneratorAction* gunGen)
: G4UImessenger(),
  fGunGeneratorAction(gunGen)
{
  // Create directories for commands
  fDirectory = new G4UIdirectory("/particle/");
  fDirectory->SetGuidance("Particle gun control commands.");
  
  // Command to select particle type
  fParticleTypeCmd = new G4UIcmdWithAString("/particle/particleType", this);
  fParticleTypeCmd->SetGuidance("Select the primary particle type.");
  fParticleTypeCmd->SetParameterName("Type", false);
  fParticleTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  // Command to set momentum
  fMomentumCmd = new G4UIcmdWithADoubleAndUnit("/particle/momentum", this);
  fMomentumCmd->SetGuidance("Set the primary particle momentum.");
  fMomentumCmd->SetParameterName("Momentum", false);
  fMomentumCmd->SetUnitCategory("Energy");
  fMomentumCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set momentum spread
  fMomentumSpreadCmd = new G4UIcmdWithADoubleAndUnit("/particle/momentumSpread", this);
  fMomentumSpreadCmd->SetGuidance("Set the primary particle momentum spread.");
  fMomentumSpreadCmd->SetParameterName("MomentumSpread", false);
  fMomentumSpreadCmd->SetUnitCategory("Energy");
  fMomentumSpreadCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set momentum distribution type
  fMomentumDistCmd = new G4UIcmdWithAString("/particle/momentumDistribution", this);
  fMomentumDistCmd->SetGuidance("Select the momentum distribution.");
  fMomentumDistCmd->SetGuidance("  Choice : uniform, gaussian");
  fMomentumDistCmd->SetParameterName("MomentumDist", false);
  fMomentumDistCmd->SetCandidates("uniform gaussian");
  fMomentumDistCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set initial position
  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/particle/position", this);
  fPositionCmd->SetGuidance("Set the starting position of the particle.");
  fPositionCmd->SetParameterName("X", "Y", "Z", false);
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set initial position spread
  fPositionSpreadCmd = new G4UIcmdWith3VectorAndUnit("/particle/positionSpread", this);
  fPositionSpreadCmd->SetGuidance("Set the starting position spread of the particle.");
  fPositionSpreadCmd->SetParameterName("DX", "DY", "DZ", false);
  fPositionSpreadCmd->SetUnitCategory("Length");
  fPositionSpreadCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set initial position distribution
  fPositionDistCmd = new G4UIcmdWithAString("/particle/positionDistribution", this);
  fPositionDistCmd->SetGuidance("Select the initial position distribution.");
  fPositionDistCmd->SetGuidance("  Choice : uniform, gaussian");
  fPositionDistCmd->SetParameterName("PositionDist", false);
  fPositionDistCmd->SetCandidates("uniform gaussian");
  fPositionDistCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set XZ angle
  fXZAngleCmd = new G4UIcmdWithADoubleAndUnit("/particle/angleXZ", this);
  fXZAngleCmd->SetGuidance("Set the primary particle initial XZ angle.");
  fXZAngleCmd->SetParameterName("angleXZ", false);
  fXZAngleCmd->SetUnitCategory("Angle");
  fXZAngleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set XZ angle spread
  fXZAngleSpreadCmd = new G4UIcmdWithADoubleAndUnit("/particle/angleXZSpread", this);
  fXZAngleSpreadCmd->SetGuidance("Set the primary particle initial XZ angle spread.");
  fXZAngleSpreadCmd->SetParameterName("angleXZSpread", false);
  fXZAngleSpreadCmd->SetUnitCategory("Angle");
  fXZAngleSpreadCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set XZ angle distribution type
  fXZAngleDistCmd = new G4UIcmdWithAString("/particle/angleXZDistribution", this);
  fXZAngleDistCmd->SetGuidance("Select the XZ angle distribution.");
  fXZAngleDistCmd->SetGuidance("  Choice : uniform, gaussian");
  fXZAngleDistCmd->SetParameterName("angleXZDist", false);
  fXZAngleDistCmd->SetCandidates("uniform gaussian");
  fXZAngleDistCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set XY angle
  fXYAngleCmd = new G4UIcmdWithADoubleAndUnit("/particle/angleXY", this);
  fXYAngleCmd->SetGuidance("Set the primary particle initial XY angle.");
  fXYAngleCmd->SetParameterName("angleXY", false);
  fXYAngleCmd->SetUnitCategory("Angle");
  fXYAngleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set XY angle spread
  fXYAngleSpreadCmd = new G4UIcmdWithADoubleAndUnit("/particle/angleXYSpread", this);
  fXYAngleSpreadCmd->SetGuidance("Set the primary particle initial XY angle spread.");
  fXYAngleSpreadCmd->SetParameterName("angleXYSpread", false);
  fXYAngleSpreadCmd->SetUnitCategory("Angle");
  fXYAngleSpreadCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Command to set XY angle distribution type
  fXYAngleDistCmd = new G4UIcmdWithAString("/particle/angleXYDistribution", this);
  fXYAngleDistCmd->SetGuidance("Select the XY angle distribution.");
  fXYAngleDistCmd->SetGuidance("  Choice : uniform, gaussian");
  fXYAngleDistCmd->SetParameterName("angleXYDist", false);
  fXYAngleDistCmd->SetCandidates("uniform gaussian");
  fXYAngleDistCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

GunGeneratorMessenger::~GunGeneratorMessenger()
{
  delete fParticleTypeCmd;
  delete fMomentumCmd;
  delete fMomentumSpreadCmd;
  delete fMomentumDistCmd;
  delete fPositionCmd;
  delete fPositionSpreadCmd;
  delete fPositionDistCmd;
  delete fXZAngleCmd;
  delete fXZAngleSpreadCmd;
  delete fXZAngleDistCmd;
  delete fXYAngleCmd;
  delete fXYAngleSpreadCmd;
  delete fXYAngleDistCmd;
  delete fDirectory;
}

void GunGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fParticleTypeCmd) {
    fGunGeneratorAction->SetParticleType(newValue);
  }
  else if (command == fMomentumCmd) {
    G4cout << "fMomentumCmd" << G4endl;
    //G4double value = G4UIcommand::ConvertToDimensionedDouble(newValue);
    G4double value = fMomentumCmd->GetNewDoubleValue(newValue);
    G4cout << "Input value: " << G4BestUnit(value, "Energy") << G4endl;
    fGunGeneratorAction->SetMomentum(value);
  }
  else if (command == fMomentumSpreadCmd) {
    G4double value = fMomentumSpreadCmd->GetNewDoubleValue(newValue);
    fGunGeneratorAction->SetMomentumSpread(value);
  }
  else if (command == fMomentumDistCmd) {
    fGunGeneratorAction->SetMomentumDist(newValue);
  }
  else if (command == fPositionCmd) {
    G4ThreeVector value = fPositionCmd->GetNew3VectorValue(newValue);
    fGunGeneratorAction->SetPosition(value);
  }
  else if (command == fPositionSpreadCmd) {
    G4ThreeVector value = fPositionSpreadCmd->GetNew3VectorValue(newValue);
    fGunGeneratorAction->SetPositionSpread(value);
  }
  else if (command == fPositionDistCmd) {
    fGunGeneratorAction->SetPositionDist(newValue);
  }
  else if (command == fXZAngleCmd) {
    G4double value = fXZAngleCmd->GetNewDoubleValue(newValue);
    fGunGeneratorAction->SetAngleXZ(value);
  }
  else if (command == fXZAngleSpreadCmd) {
    G4double value = fXZAngleSpreadCmd->GetNewDoubleValue(newValue);
    fGunGeneratorAction->SetAngleXZSpread(value);
  }
  else if (command == fXZAngleDistCmd) {
    fGunGeneratorAction->SetAngleXZDist(newValue);
  }
  else if (command == fXYAngleCmd) {
    G4double value = fXYAngleCmd->GetNewDoubleValue(newValue);
    fGunGeneratorAction->SetAngleXY(value);
  }
  else if (command == fXYAngleSpreadCmd) {
    G4double value = fXYAngleSpreadCmd->GetNewDoubleValue(newValue);
    fGunGeneratorAction->SetAngleXYSpread(value);
  }
  else if (command == fXYAngleDistCmd) {
    fGunGeneratorAction->SetAngleXYDist(newValue);
  }
}