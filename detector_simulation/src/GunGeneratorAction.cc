#include "GunGeneratorAction.hh"
#include "GunGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"

GunGeneratorAction::GunGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr),
  fMessenger(nullptr),
  fParticleName("mu+"),
  fMomentum(0.25*GeV),
  fMomentumSpread(0.00*GeV),
  fMomentumDist("uniform"),
  fPosition(G4ThreeVector(0.0, 0.0, 0.0)*m),
  fPositionSpread(G4ThreeVector(0.0, 0.0, 0.0)*m),
  fPositionDist("uniform"),
  fXZAngle(0.00*deg),
  fXZAngleSpread(0.00*deg),
  fXZAngleDist("uniform"),
  fXYAngle(0.00*deg),
  fXYAngleSpread(0.00*deg),
  fXYAngleDist("uniform")
{
    // Create particle gun
    fParticleGun = new G4ParticleGun(1);
    
    // Define command interface
    DefineCommands();

    // Initialize the particle gun with default values
    Update();
}

GunGeneratorAction::~GunGeneratorAction()
{
    delete fParticleGun;
    delete fMessenger;
}

void GunGeneratorAction::GeneratePrimaries(G4Event* event)
{   
    Update();

    G4cout << "Generating " << fParticleName 
           << " with momentum " << G4BestUnit(fParticleGun->GetParticleMomentum(), "Energy")
           << " at position (" << G4BestUnit(fParticleGun->GetParticlePosition().x(), "Length") 
           << ", " << G4BestUnit(fParticleGun->GetParticlePosition().y(), "Length")
           << ", " << G4BestUnit(fParticleGun->GetParticlePosition().z(), "Length") << ")"
           << " with direction (" << fParticleGun->GetParticleMomentumDirection().x() 
           << ", " << fParticleGun->GetParticleMomentumDirection().y()
           << ", " << fParticleGun->GetParticleMomentumDirection().z() << ")" << G4endl;
    
    // Generate primary vertex
    fParticleGun->GeneratePrimaryVertex(event);
}

void GunGeneratorAction::Update()
{
    // Set particle type
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(fParticleName);
    if (particle) {
        fParticleGun->SetParticleDefinition(particle);
    }
    
    // Set momentum
    G4double momentum = RandomScalar(fMomentum, fMomentumSpread, fMomentumDist);
    fParticleGun->SetParticleMomentum(momentum);

    // Set initial position
    G4ThreeVector position = RandomVector(fPosition, fPositionSpread, fPositionDist);
    fParticleGun->SetParticlePosition(position);

    // Set direction
    G4double angleXZ = RandomScalar(fXZAngle, fXZAngleSpread, fXZAngleDist);
    G4double angleXY = RandomScalar(fXYAngle, fXYAngleSpread, fXYAngleDist);
    G4ThreeVector dir(sin(angleXZ)*sin(angleXY), sin(angleXZ)*cos(angleXY), cos(angleXZ));
    fParticleGun->SetParticleMomentumDirection(dir);
}

void GunGeneratorAction::SetParticleType(G4String particleName)
{
    // Find the particle by name
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
    
    if (particle) {
        fParticleName = particleName;
        G4cout << "Particle type set to " << fParticleName << G4endl;
    } else {
        G4cout << "WARNING: Particle " << particleName << " not found!" << G4endl;
    }
}

void GunGeneratorAction::SetMomentum(G4double momentum)
{
    fMomentum = momentum;
    G4cout << "Particle momentum set to " << G4BestUnit(fMomentum, "Energy") << G4endl;
}

void GunGeneratorAction::SetMomentumSpread(G4double momentumSpread)
{
    fMomentumSpread = momentumSpread;
    G4cout << "Particle momentum spread set to " << G4BestUnit(fMomentumSpread, "Energy") << G4endl;
}

void GunGeneratorAction::SetMomentumDist(G4String momentumDist)
{
    fMomentumDist = momentumDist;
    G4cout << "Particle momentum distribution set to " << fMomentumDist << G4endl;
}

void GunGeneratorAction::SetPosition(G4ThreeVector position)
{
    fPosition = position;
    G4cout << "Particle position set to (" 
           << G4BestUnit(fPosition.x(), "Length") << ", " 
           << G4BestUnit(fPosition.y(), "Length") << ", " 
           << G4BestUnit(fPosition.z(), "Length") << ")" << G4endl;
}

void GunGeneratorAction::SetPositionSpread(G4ThreeVector positionSpread)
{
    fPositionSpread = positionSpread;
    G4cout << "Particle position spread set to (" 
           << G4BestUnit(fPositionSpread.x(), "Length") << ", " 
           << G4BestUnit(fPositionSpread.y(), "Length") << ", " 
           << G4BestUnit(fPositionSpread.z(), "Length") << ")" << G4endl;
}

void GunGeneratorAction::SetPositionDist(G4String positionDist)
{
    fPositionDist = positionDist;
    G4cout << "Particle position distribution set to " << fPositionDist << G4endl;
}

void GunGeneratorAction::SetAngleXZ(G4double angleXZ)
{
    fXZAngle = angleXZ;
    G4cout << "Particle XZ angle set to " << G4BestUnit(fXZAngle, "Angle") << G4endl;
}

void GunGeneratorAction::SetAngleXZSpread(G4double angleXZSpread)
{
    fXZAngleSpread = angleXZSpread;
    G4cout << "Particle XZ angle spread set to " << G4BestUnit(fXZAngleSpread, "Angle") << G4endl;
}

void GunGeneratorAction::SetAngleXZDist(G4String angleXZDist)
{
    fXZAngleDist = angleXZDist;
    G4cout << "Particle XZ angle distribution set to " << fXZAngleDist << G4endl;
}

void GunGeneratorAction::SetAngleXY(G4double angleXY)
{
    fXYAngle = angleXY;
    G4cout << "Particle XY angle set to " << G4BestUnit(fXYAngle, "Angle") << G4endl;
}

void GunGeneratorAction::SetAngleXYSpread(G4double angleXYSpread)
{
    fXYAngleSpread = angleXYSpread;
    G4cout << "Particle XY angle spread set to " << G4BestUnit(fXYAngleSpread, "Angle") << G4endl;
}

void GunGeneratorAction::SetAngleXYDist(G4String angleXYDist)
{
    fXYAngleDist = angleXYDist;
    G4cout << "Particle XY angle distribution set to " << fXYAngleDist << G4endl;
}

G4double GunGeneratorAction::RandomScalar(G4double central_value, G4double spread, G4String dist)
{
    G4double ret;
    if (dist == "uniform") {
        ret = central_value + (2.0*G4UniformRand() - 1.0) * spread;
    } else if (dist == "gaussian") {
        ret = G4RandGauss::shoot(central_value, spread);
    }
    return ret;
}

G4ThreeVector GunGeneratorAction::RandomVector(G4ThreeVector central_value, G4ThreeVector spread, G4String dist)
{
    G4ThreeVector ret;
    ret.setX(RandomScalar(central_value.x(), spread.x(), dist));
    ret.setY(RandomScalar(central_value.y(), spread.y(), dist));
    ret.setZ(RandomScalar(central_value.z(), spread.z(), dist));
    return ret;
}

void GunGeneratorAction::DefineCommands()
{
    // Create a new messenger
    fMessenger = new GunGeneratorMessenger(this);
}