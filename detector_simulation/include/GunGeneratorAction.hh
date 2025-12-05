#ifndef GUNGENERATORACTION_HH
#define GUNGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4ParticleGun;
class GunGeneratorMessenger;

class GunGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    GunGeneratorAction();
    virtual ~GunGeneratorAction();
    
    // Method from base class
    virtual void GeneratePrimaries(G4Event* event);
    void Update();
    
    // Setter methods for particle gun properties
    void SetParticleType(G4String particleName);
    void SetMomentum(G4double momentum);
    void SetMomentumSpread(G4double momentumSpread);
    void SetMomentumDist(G4String momentumDist);
    void SetPosition(G4ThreeVector position);
    void SetPositionSpread(G4ThreeVector positionSpread);
    void SetPositionDist(G4String positionDist);
    void SetPositionRMax(G4double rmax);
    void SetAngleXZ(G4double angleXZ);
    void SetAngleXZSpread(G4double angleXZSpread);
    void SetAngleXZDist(G4String angleXZDist);
    void SetAngleXY(G4double angleXY);
    void SetAngleXYSpread(G4double angleXYSpread);
    void SetAngleXYDist(G4String angleXYDist);
    
    // Getter methods
    G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    
private:
    G4ParticleGun* fParticleGun;
    GunGeneratorMessenger* fMessenger;

    // Random numbers
    G4double RandomScalar(G4double central_value, G4double spread, G4String dist);
    G4ThreeVector RandomVector(G4ThreeVector central_value, G4ThreeVector spread, G4String dist, G4double rmax=-999.0);
    G4ThreeVector RandomVectorInCylinder(G4ThreeVector central_value, G4ThreeVector spread, G4double rmax);
    G4double GetDistOverlap(G4ThreeVector center, G4ThreeVector spread, G4String distribution, G4double rmax);
    
    // Configurable parameters with default values
    G4String fParticleName;
    G4double fMomentum;
    G4double fMomentumSpread;
    G4String fMomentumDist;
    G4ThreeVector fPosition;
    G4ThreeVector fPositionSpread;
    G4String fPositionDist;
    G4double fPositionRMax;
    G4double fXZAngle;
    G4double fXZAngleSpread;
    G4String fXZAngleDist;
    G4double fXYAngle;
    G4double fXYAngleSpread;
    G4String fXYAngleDist;

    bool firstEvent = true;

    // Helper methods
    void DefineCommands();
};

#endif
