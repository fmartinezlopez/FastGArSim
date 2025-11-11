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
    G4VPhysicalVolume* GetWorldVolume() const { return fWorldPhysical; }
    G4LogicalVolume* GetTPCSensLogical()  const { return fTPCLogical; }
    G4LogicalVolume* GetECalSensLogical() const { return fECalScintillatorLogical; }
    G4LogicalVolume* GetMuIDSensLogical() const { return fMuIDScintillatorLogical; }
    
    // Messenger methods for configurable parameters
    void SetGeometryType(G4String type);
    void SetMagneticFieldStrength(G4double bfield);
    void SetTPCRadius(G4double radius);
    void SetTPCLength(G4double length);
    void SetECalHGAbsorberThickness(G4double thickness);
    void SetECalHGScintillatorThickness(G4double thickness);
    void SetECalHGBoardThickness(G4double thickness);
    void SetECalLGAbsorberThickness(G4double thickness);
    void SetECalLGScintillatorThickness(G4double thickness);
    void SetECalBarrelHGLayers(G4int layers);
    void SetECalBarrelLGLayers(G4int layers);
    void SetECalEndcapHGLayers(G4int layers);
    void SetECalEndcapLGLayers(G4int layers);
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
    void SetPressure(G4double pressure);

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
    void ConstructSamplingBarrel(G4String baseName,                        // base name for volumes
                                 G4double barrelLength,                    // total length of barrel
                                 G4double barrelInnerDistance,             // apothem of inner polygon
                                 G4int numSides,                           // number of sides
                                 G4double totalThickness,                  // total thickness of barrel
                                 G4double layerHGThickness,                // HG layer thickness
                                 G4double layerHGAbsorberThickness,        // HG absorber thickness
                                 G4double layerHGScintillatorThickness,    // HG scintillator thickness
                                 G4double layerLGThickness,                // LG layer thickness
                                 G4double layerLGAbsorberThickness,        // LG absorber thickness
                                 G4int numHGLayers,                        // number of HG layers
                                 G4int numLGLayers,                        // number of LG layers
                                 G4Material* absorberMaterial,             // absorber material
                                 G4Material* scintillatorMaterial,         // scintillator material
                                 G4Material* boardMaterial,                // PCB material
                                 G4LogicalVolume* parentVolume,            // parent logical volume
                                 G4LogicalVolume** outVolume,              // logical volume
                                 G4LogicalVolume** outScintillatorVolume,  // scintillator logical volume
                                 G4Colour visColor                         // colour for visualisation
                                 );

    void ConstructSamplingEndcap(G4String baseName,                        // base name for volumes
                                 G4double endcapStart,                     // detector start in Z
                                 G4double endcapRadius,                    // radius for circular endcap
                                 G4double totalThickness,                  // total thickness of barrel
                                 G4double layerHGThickness,                // HG layer thickness
                                 G4double layerHGAbsorberThickness,        // HG absorber thickness
                                 G4double layerHGScintillatorThickness,    // HG scintillator thickness
                                 G4double layerLGThickness,                // LG layer thickness
                                 G4double layerLGAbsorberThickness,        // LG absorber thickness
                                 G4int numHGLayers,                        // number of HG layers
                                 G4int numLGLayers,                        // number of LG layers
                                 G4Material* absorberMaterial,             // absorber material
                                 G4Material* scintillatorMaterial,         // scintillator material
                                 G4Material* boardMaterial,                // PCB material
                                 G4LogicalVolume* parentVolume,            // parent logical volume
                                 G4LogicalVolume** outVolume,              // logical volume
                                 G4LogicalVolume** outScintillatorVolume,  // scintillator logical volume
                                 G4Colour visColor                         // colour for visualisation
                                 );

    // Helper methods
    void DefineCommands();
    
    // Materials
    G4Material* fWorldMaterial;
    G4Material* fGArTPCMaterial;
    G4Material* fECalAbsorberMaterial;
    G4Material* fECalScintillatorMaterial;
    G4Material* fECalPCBMaterial;
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
    
    // Detector configuration parameters (GAr)
    G4double fTPCRadius;                       // radius of TPC
    G4double fTPCLength;                       // length of TPC
    G4double fECalBarrelGap;                   // distance between TPC and inner apothem of ECal barrel
    G4double fECalEndcapGap;                   // distance between TPC and start of ECal endcap (drift direction)
    G4double fECalHGAbsorberThickness;         // absorber thickness in high-granularity ECal layers
    G4double fECalHGScintillatorThickness;     // scintillator thickness in high-granularity ECal layers
    G4double fECalHGBoardThickness;            // PCB thickness in high-granularity ECal layers (for tiles only)
    G4double fECalLGAbsorberThickness;         // absorber thickness in low-granularity ECal layers
    G4double fECalLGScintillatorThickness;     // scintillator thickness in low-granularity ECal layers
    G4int    fECalNumSides;                    // number of sides for ECal barrel polyhedron
    G4int    fECalBarrelHGLayers;              // number of ECal barrel high-granularity layers
    G4int    fECalBarrelLGLayers;              // number of ECal barrel low-granularity layers
    G4int    fECalEndcapHGLayers;              // number of ECal end cap high-granularity layers
    G4int    fECalEndcapLGLayers;              // number of ECal end cap low-granularity layers
    G4double fMuIDBarrelGap;                   // distance between TPC and inner apothem of MuID barrel
    G4double fMuIDAbsorberThickness;           // absorber thickness in MuID layers
    G4double fMuIDScintillatorThickness;       // scintillator thickness in MuID layers
    G4int    fMuIDNumSides;                    // number of sides for MuID barrel polyhedron
    G4int    fMuIDLayers;                      // number of MuID layers

    // Detector configuration parameters (LAr)
    G4int    fLArNModulesX;                    // number of TPC modules in X direction (drift direction)
    G4int    fLArNModulesY;                    // number of TPC modules in Y direction (vertical direction)
    G4int    fLArNModulesZ;                    // number of TPC modules in Z direction (beam direction)
    G4double fLArModuleLength;                 // length of TPC modules (drift direction)
    G4double fLArModuleWidth;                  // width of TPC modules (vertical direction)
    G4double fLArModuleDepth;                  // depth of TPC modules (beam direction)
    G4double fLArModuleGap;                    // distance between TPC modules
    G4double fLArInsulationThickness;          // thickness of optical isolation layer
    G4double fLArCryostatThickness;            // thickness of steel cryostat
    G4bool   fLArEnableMuonWindow;             // add low-density downstream muon window?
    G4double fLArMuonWindowThickness;          // low-density muon window thickness

    // Derived configuration parameters (GAr)
    G4double fECalHGLayerThickness;            // thickness of high-granularity ECal layers
    G4double fECalLGLayerThickness;            // thickness of low-granularity ECal layers
    G4double fMuIDLayerThickness;              // thickness of MuID layers
    G4double fECalBarrelTotalThickness;        // ECal barrel total thickness
    G4double fECalEndcapTotalThickness;        // ECal end cap total thickness
    G4double fMuIDTotalThickness;              // MuID total thickness

    // Derived configuration parameters (LAr)
    G4double fLArTotalLength;                  // LAr TPC array total length (drift direction)
    G4double fLArTotalWidth;                   // LAr TPC array total width (vertical direction)
    G4double fLArTotalDepth;                   // LAr TPC array total depth (beam direction)

    // Material properties
    G4double fPressure;                        // pressure in HPgTPC
	G4double fRefPressure;                     // reference pressure for gas density
	G4double fTemperature;                     // pressure in HPgTPC
    G4double fGasDensity;                      // gas density at reference conditions
    G4double fFoamDensity;                     // muon window density (LAr)

    // Messenger for macro commands
    DetectorMessenger* fMessenger;
};

#endif
