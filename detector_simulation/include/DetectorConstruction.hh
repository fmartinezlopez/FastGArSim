#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Cache.hh"
#include "G4GenericMessenger.hh"

class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    G4bool UpdateGeometry();
    
    // Getter methods
    G4LogicalVolume* GetTPCSensLogical()  const { return fTPCLogical; }
    G4LogicalVolume* GetECalSensLogical() const { return fECalScintillatorLogical; }
    G4LogicalVolume* GetMuIDSensLogical() const { return fMuIDScintillatorLogical; }
    
    // Messenger methods for configurable parameters
    void SetTPCRadius(G4double radius);
    void SetTPCLength(G4double length);
    void SetECalAbsorberThickness(G4double thickness);
    void SetECalScintillatorThickness(G4double thickness);
    void SetECalLayers(G4int layers);
    void SetMuIDAbsorberThickness(G4double thickness);
    void SetMuIDScintillatorThickness(G4double thickness);
    void SetMuIDLayers(G4int layers);
    
private:

    G4bool fGeometryInitialized;

    // Methods to create detector components
    void DefineMaterials();
    G4VPhysicalVolume* ConstructWorld();
    void ConstructFieldEnclosure();
    void ConstructTPC();
    void ConstructECal();
    void ConstructMuID();
    
    // Helper construction methods
    void ConstructSamplingBarrel(G4String baseName,                       // base name for volumes
                                 G4double barrelLength,                   // total length of barrel
                                 G4double barrelInnerDistance,            // apothem of inner polygon
                                 G4int numSides,                          // number of sides
                                 G4double totalThickness,                 // total thickness of barrel
                                 G4double layerThickness,                 // layer thickness
                                 G4double layerAbsorberThickness,         // absorber thickness
                                 G4int numLayers,                         // number of layers
                                 G4Material* absorberMaterial,            // absorber material
                                 G4Material* scintillatorMaterial,        // scintillator material
                                 G4LogicalVolume* parentVolume,           // parent logical volume
                                 G4LogicalVolume** outVolume,
                                 G4LogicalVolume** outScintillatorVolume,  // scintillator logical volume
                                 G4Colour visColor
                                 );

    void ConstructSamplingEndcap(G4String baseName,
                                 G4double endcapStart,
                                 G4double endcapRadius,
                                 G4double totalThickness,
                                 G4double layerThickness,
                                 G4double layerAbsorberThickness,
                                 G4int numLayers,
                                 G4Material* absorberMaterial,
                                 G4Material* scintillatorMaterial,
                                 G4LogicalVolume* parentVolume,
                                 G4LogicalVolume** outVolume,
                                 G4LogicalVolume** outScintillatorVolume,
                                 G4Colour visColor
                                 );

    // Helper methods
    void DefineCommands();
    
    // Materials
    G4Material* fWorldMaterial;
    G4Material* fTPCMaterial;
    G4Material* fECalAbsorberMaterial;
    G4Material* fECalScintillatorMaterial;
    G4Material* fMuIDScintillatorMaterial;
    G4Material* fMuIDAbsorberMaterial;
    
    // Logical volumes
    G4LogicalVolume* fWorldLogical;
    G4LogicalVolume* fFieldLogical;
    G4LogicalVolume* fTPCLogical;
    G4LogicalVolume* fECalBarrelLogical;
    G4LogicalVolume* fECalEndcapsLogical;
    G4LogicalVolume* fECalScintillatorLogical;
    G4LogicalVolume* fMuIDLogical;
    G4LogicalVolume* fMuIDScintillatorLogical;
    
    // Physical volumes
    G4VPhysicalVolume* fWorldPhysical;
    G4VPhysicalVolume* fFieldPhysical;
    G4VPhysicalVolume* fTPCPhysical;

    // Magnetic field
    G4UniformMagField* fMagneticField;
    G4double fMagneticFieldStrength;
    
    // Configuration parameters with default values
    G4double fTPCRadius;
    G4double fTPCLength;
    G4double fECalBarrelGap;
    G4double fECalEndcapGap;
    G4double fECalAbsorberThickness;
    G4double fECalScintillatorThickness;
    G4int    fECalNumSides;
    G4int    fECalLayers;
    G4double fMuIDBarrelGap;
    G4double fMuIDAbsorberThickness;
    G4double fMuIDScintillatorThickness;
    G4int    fMuIDNumSides;
    G4int    fMuIDLayers;

    // Derived configuration parameters
    G4double fECalLayerThickness;
    G4double fMuIDLayerThickness;
    G4double fECalTotalThickness;
    G4double fMuIDTotalThickness;
    
    // Messenger for macro commands
    DetectorMessenger* fMessenger;
};

#endif