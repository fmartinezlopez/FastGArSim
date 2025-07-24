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

    // Enum for geometry types
    enum GeometryType {
        kGArLike, // magnetised GAr TPC + ECal + MuID
        kLArLike  // modular LAr TPC
    };

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
    void SetGeometryType(G4String type);
    void SetTPCRadius(G4double radius);
    void SetTPCLength(G4double length);
    void SetECalAbsorberThickness(G4double thickness);
    void SetECalScintillatorThickness(G4double thickness);
    void SetECalLayers(G4int layers);
    void SetMuIDAbsorberThickness(G4double thickness);
    void SetMuIDScintillatorThickness(G4double thickness);
    void SetMuIDLayers(G4int layers);
    void SetLArNModulesX(G4int nx);
    void SetLArNModulesY(G4int ny);
    void SetLArNModulesZ(G4int nz);
    void SetLArModuleLength(G4double length);
    void SetLArModuleWidth(G4double width);
    void SetLArModuleDepth(G4double depth);
    void SetLArModuleGap(G4double gap);
    void SetLArInsulationThickness(G4double thickness);
    void SetLArCryostatThickness(G4double thickness);
    void SetLArEnableMuonWindow(G4bool enable);
    void SetLArMuonWindowThickness(G4double thickness);

private:

    G4bool fGeometryInitialized;
    GeometryType fGeometryType;

    void ComputeDerivedQuantities();

    // Methods to create detector components
    void DefineMaterials();
    G4VPhysicalVolume* ConstructWorld();

    // GAr detector construction methods
    void ConstructGArDetector();
    void ConstructFieldEnclosure();
    void ConstructTPC();
    void ConstructECal();
    void ConstructMuID();

    // LAr detector construction methods
    void ConstructLArDetector();
    
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
    G4Material* fGArTPCMaterial;
    G4Material* fECalAbsorberMaterial;
    G4Material* fECalScintillatorMaterial;
    G4Material* fMuIDScintillatorMaterial;
    G4Material* fMuIDAbsorberMaterial;
    G4Material* fLArTPCMaterial;
    G4Material* fLArCryostatMaterial;
    G4Material* fLArInsulationMaterial;
    G4Material* fLArMuonWindowMaterial;
    
    // Logical volumes
    G4LogicalVolume* fWorldLogical;
    G4LogicalVolume* fFieldLogical;
    G4LogicalVolume* fTPCLogical;
    G4LogicalVolume* fECalBarrelLogical;
    G4LogicalVolume* fECalEndcapsLogical;
    G4LogicalVolume* fECalScintillatorLogical;
    G4LogicalVolume* fMuIDLogical;
    G4LogicalVolume* fMuIDScintillatorLogical;
    G4LogicalVolume* fLArTPCLogical;
    
    // Physical volumes
    G4VPhysicalVolume* fWorldPhysical;
    G4VPhysicalVolume* fFieldPhysical;
    G4VPhysicalVolume* fTPCPhysical;
    G4VPhysicalVolume* fLArTPCPhysical;

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
    G4int    fLArNModulesX;
    G4int    fLArNModulesY;
    G4int    fLArNModulesZ;
    G4double fLArModuleLength;
    G4double fLArModuleWidth;
    G4double fLArModuleDepth;
    G4double fLArModuleGap;
    G4double fLArInsulationThickness;
    G4double fLArCryostatThickness;
    G4bool   fLArEnableMuonWindow;
    G4double fLArMuonWindowThickness;

    // Derived configuration parameters
    G4double fECalLayerThickness;
    G4double fMuIDLayerThickness;
    G4double fECalTotalThickness;
    G4double fMuIDTotalThickness;
    G4double fLArTotalLength;
    G4double fLArTotalWidth;
    G4double fLArTotalDepth;
    
    // Messenger for macro commands
    DetectorMessenger* fMessenger;
};

#endif