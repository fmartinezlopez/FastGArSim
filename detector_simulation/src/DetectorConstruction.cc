#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PhysicalConstants.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"

#include <cmath>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fWorldLogical(nullptr), fFieldLogical(nullptr), fTPCLogical(nullptr),
  fECalBarrelLogical(nullptr), fECalEndcapsLogical(nullptr), fECalScintillatorLogical(nullptr),
  fMuIDLogical(nullptr), fMuIDScintillatorLogical(nullptr),
  fWorldPhysical(nullptr), fFieldPhysical(nullptr), fTPCPhysical(nullptr),
  fMagneticField(nullptr), fMagneticFieldStrength(0.5*tesla),
  fTPCRadius(250.0*cm), fTPCLength(500.0*cm),
  fECalBarrelGap(20.0*cm), fECalEndcapGap(25.0*cm),
  fECalAbsorberThickness(5.0*mm), fECalScintillatorThickness(10.0*mm), fECalNumSides(12), fECalLayers(42),
  fMuIDBarrelGap(50.0*cm),
  fMuIDAbsorberThickness(10.0*cm), fMuIDScintillatorThickness(2.0*cm), fMuIDNumSides(21), fMuIDLayers(3),
  fMessenger(nullptr),
  fGeometryInitialized(false)
{

    // Define command interface
    DefineCommands();
    
    // Define materials
    DefineMaterials();

    G4cout << "TPC radius: " << G4BestUnit(fTPCRadius, "Length") << "\n"
           << "TPC length: " << G4BestUnit(fTPCLength, "Length") << "\n"
           << G4endl;

    // Define some numerical quantities that will be used extensively
    fECalLayerThickness = fECalAbsorberThickness + fECalScintillatorThickness;
    fMuIDLayerThickness = fMuIDAbsorberThickness + fMuIDScintillatorThickness;

    fECalTotalThickness = fECalLayerThickness*fECalLayers;
    fMuIDTotalThickness = fMuIDLayerThickness*fMuIDLayers;
}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

G4bool DetectorConstruction::UpdateGeometry()
{
    // Recalculate derived quantities
    fECalLayerThickness = fECalAbsorberThickness + fECalScintillatorThickness;
    fMuIDLayerThickness = fMuIDAbsorberThickness + fMuIDScintillatorThickness;
    fECalTotalThickness = fECalLayerThickness*fECalLayers;
    fMuIDTotalThickness = fMuIDLayerThickness*fMuIDLayers;
    
    G4cout << "Updating geometry..." << G4endl;
    
    // Tell G4RunManager to rebuild geometry
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    return true;
}

void DetectorConstruction::DefineMaterials()
{

    // World material definition
	G4NistManager *nistManager = G4NistManager::Instance();
	
    // Define materials for different detector components
    fWorldMaterial = nistManager->FindOrBuildMaterial("G4_AIR");

	// Reference temperature and pressure for gas
	G4double fPressure    =  10.13; // bar
	G4double fAtmPressure =   1.01; // bar
	G4double fTemperature = 294.26; // K

	// Gas mixture density
	G4double p10Density = 0.001677*g/cm3; // at atm pressure (used in GArSoft -- need to check)
    //G4double p10Density = 0.00159*g/cm3; // at atm pressure

	// Define needed elements
	G4double z, fractionmass;
	G4int nel, natoms;

	G4Element* H  = new G4Element("Hydrogen", "H",  z = 1.,  1.008*g/mole);
	G4Element* C  = new G4Element("Carbon",   "C",  z = 6.,  12.011*g/mole);
	G4Element* Ar = new G4Element("Argon",    "Ar", z = 18., 39.948*g/mole);

	// Define gas material
	fTPCMaterial = new G4Material("P10", p10Density*fPressure/fAtmPressure, nel = 3, kStateGas, fTemperature*kelvin, fPressure*bar);

	// Add elements to material
	fTPCMaterial->AddElement(H,  fractionmass = 0.011);
	fTPCMaterial->AddElement(C,  fractionmass = 0.032);
	fTPCMaterial->AddElement(Ar, fractionmass = 0.957);

    // Materials for ECal
    fECalAbsorberMaterial = nistManager->FindOrBuildMaterial("G4_Pb");
    fECalScintillatorMaterial = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    // Materials for MuID
    fMuIDAbsorberMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    fMuIDScintillatorMaterial = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Construct world volume
    fWorldPhysical = ConstructWorld();
    ConstructFieldEnclosure();

    // Construct detector components
    ConstructTPC();
    ConstructECal();
    ConstructMuID();

    G4UserLimits* limitsTPCSteps = new G4UserLimits(5.0*mm);
    G4UserLimits* limitsCaloSteps = new G4UserLimits(5.0*mm);

    fTPCLogical->SetUserLimits(limitsTPCSteps);
    fECalBarrelLogical->SetUserLimits(limitsCaloSteps);
    fECalEndcapsLogical->SetUserLimits(limitsCaloSteps);
    fMuIDLogical->SetUserLimits(limitsCaloSteps);
    
    fGeometryInitialized = true;
    return fWorldPhysical;
}

G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{
    // Calculate world size based on detector dimensions
    G4double worldSizeXY = 2.5 * (fTPCRadius + fECalTotalThickness + fMuIDTotalThickness);
    G4double worldSizeZ  = 2.5 * (fTPCLength/2 + fECalTotalThickness + fMuIDTotalThickness);
    
    // Create world box
    G4Box* worldSolid = new G4Box("World", worldSizeXY, worldSizeXY, worldSizeZ);
    fWorldLogical = new G4LogicalVolume(worldSolid, fWorldMaterial, "World_log");
    
    // Set world visualization attributes
    G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    worldVisAtt->SetVisibility(false);
    fWorldLogical->SetVisAttributes(worldVisAtt);
    
    // Place world volume
    fWorldPhysical = new G4PVPlacement(0, G4ThreeVector(), fWorldLogical, "World_phys", 0, false, 0);
    
    return fWorldPhysical;
}

void DetectorConstruction::ConstructFieldEnclosure()
{
    // Calculate field enclosure size based on detector dimensions
    G4double fieldEnclosureLength = fTPCLength + 2*fECalEndcapGap + 2*fECalTotalThickness;
    G4double fieldEnclosureRadius = fTPCRadius + fECalBarrelGap + fECalTotalThickness + fMuIDBarrelGap;

    // Create field enclosure cylinder
    G4Tubs* fieldSolid = new G4Tubs("FieldEnclosure", 0, fieldEnclosureRadius, fieldEnclosureLength/2, 0, twopi);
    fFieldLogical = new G4LogicalVolume(fieldSolid, fWorldMaterial, "FieldEnclosure_log");

    // Set field enclosure visualization attributes
    G4VisAttributes* fieldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    fieldVisAtt->SetVisibility(false);
    fFieldLogical->SetVisAttributes(fieldVisAtt);

    // Place field enclosure volume
    fFieldPhysical = new G4PVPlacement(0, G4ThreeVector(), fFieldLogical, "FieldEnclosure_phys", fWorldLogical, false, 0);

}

void DetectorConstruction::ConstructTPC()
{

    // Create TPC cylindrical volume
    G4Tubs* tpcSolid = new G4Tubs("TPC", 0, fTPCRadius, fTPCLength/2, 0, twopi);
    fTPCLogical = new G4LogicalVolume(tpcSolid, fTPCMaterial, "TPC_log");
    
    // Set TPC visualization attributes
    G4VisAttributes* tpcVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    tpcVisAtt->SetVisibility(true);
    fTPCLogical->SetVisAttributes(tpcVisAtt);
    
    // Place TPC in world
    fTPCPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fTPCLogical, "TPC_phys", fFieldLogical, false, 0);
}

void DetectorConstruction::ConstructECal()
{

    G4double ecalBarrelTotalLength = fTPCLength + 2*fECalEndcapGap + 2*fECalTotalThickness;
    G4double ecalBarrelInnerDistance = fTPCRadius + fECalBarrelGap;

    ConstructSamplingBarrel("ECal", ecalBarrelTotalLength, ecalBarrelInnerDistance, fECalNumSides, fECalTotalThickness, fECalLayerThickness, fECalAbsorberThickness, fECalLayers, fECalAbsorberMaterial, fECalScintillatorMaterial, fFieldLogical, &fECalBarrelLogical, &fECalScintillatorLogical, G4Colour(0.0, 1.0, 0.0, 0.3));

    G4double ecalEndcapRadius = ecalBarrelInnerDistance / std::cos(pi / fECalNumSides);
    G4double ecalEndcapStart = fTPCLength/2 + fECalEndcapGap;

    ConstructSamplingEndcap("ECal", ecalEndcapStart, ecalEndcapRadius, fECalTotalThickness, fECalLayerThickness, fECalAbsorberThickness, fECalLayers, fECalAbsorberMaterial, fECalScintillatorMaterial, fFieldLogical, &fECalEndcapsLogical, &fECalScintillatorLogical, G4Colour(0.0, 1.0, 0.0, 0.3));
}

void DetectorConstruction::ConstructMuID()
{

    G4double muidBarrelTotalLength = fTPCLength + 2*fECalEndcapGap + 2*fECalTotalThickness;
    G4double muidBarrelInnerDistance = fTPCRadius + fECalBarrelGap + fECalTotalThickness + fMuIDBarrelGap;

    ConstructSamplingBarrel("MuID", muidBarrelTotalLength, muidBarrelInnerDistance, fMuIDNumSides, fMuIDTotalThickness, fMuIDLayerThickness, fMuIDAbsorberThickness, fMuIDLayers, fMuIDAbsorberMaterial, fMuIDScintillatorMaterial, fWorldLogical, &fMuIDLogical, &fMuIDScintillatorLogical, G4Colour(1.0, 0.0, 0.0, 0.3));

}

void DetectorConstruction::ConstructSamplingBarrel(G4String baseName,
                                                   G4double barrelLength,
                                                   G4double barrelInnerDistance,
                                                   G4int numSides,
                                                   G4double totalThickness,
                                                   G4double layerThickness,
                                                   G4double layerAbsorberThickness,
                                                   G4int numLayers,
                                                   G4Material* absorberMaterial,
                                                   G4Material* scintillatorMaterial,
                                                   G4LogicalVolume* parentVolume,
                                                   G4LogicalVolume** outVolume,
                                                   G4LogicalVolume** outScintillatorVolume,
                                                   G4Colour visColor)
{

    // Define barrel dimensions
    G4double cosineAngle = std::cos(pi / numSides);
    G4double innerApothem = barrelInnerDistance;
    G4double outerApothem = innerApothem + totalThickness;

    G4double innerRadius = innerApothem / cosineAngle;
    G4double outerRadius = outerApothem / cosineAngle;
    G4double halfLength = barrelLength / 2;

    // Create barrel
    G4double phiStart = 0;
    G4double phiTotal = twopi;
    G4double zPlanes[] = {-halfLength, halfLength};
    G4double rInner[]  = {innerRadius, innerRadius};
    G4double rOuter[]  = {outerRadius, outerRadius};
    
    G4Polyhedra* barrelSolid = new G4Polyhedra(baseName+"_barrel", phiStart, phiTotal, numSides, 2, zPlanes, rInner, rOuter);
    G4LogicalVolume* barrelLogical = new G4LogicalVolume(barrelSolid, fWorldMaterial, baseName+"_barrel_log");
    
    // Visualization attributes
    G4VisAttributes* barrelVisAtt = new G4VisAttributes(visColor);
    barrelVisAtt->SetVisibility(true);
    barrelLogical->SetVisAttributes(barrelVisAtt);
    
    // Place barrel
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), barrelLogical, baseName+"_barrel_phys", parentVolume, false, 0);

    *outVolume = barrelLogical;

    // Calculate segment angle
    G4double segmentAngle = phiTotal / numSides;

    // 
    G4VisAttributes* segmentVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    segmentVisAtt->SetVisibility(false);

    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    scintillatorVisAtt->SetVisibility(false);

    G4VisAttributes* absorberVisAtt = new G4VisAttributes(visColor);
    absorberVisAtt->SetVisibility(false);
    
    // Create segments (one for each face)
    for (G4int segmentIndex = 0; segmentIndex < numSides; segmentIndex++) {
        
        // Segment phi range
        G4double segmentPhiStart = phiStart + segmentIndex * segmentAngle;
        
        // Create segment polyhedron
        G4Polyhedra* segmentSolid = new G4Polyhedra(baseName+"_barrel_Segment", segmentPhiStart, segmentAngle, 1, 2, zPlanes, rInner, rOuter);
        G4LogicalVolume* segmentLogical = new G4LogicalVolume(segmentSolid, fWorldMaterial, baseName+"_barrel_Segment_log");

        segmentLogical->SetVisAttributes(segmentVisAtt);
        
        // Place segment
        new G4PVPlacement(nullptr, G4ThreeVector(), segmentLogical, baseName+"_barrel_Segment_phys", barrelLogical, false, segmentIndex, true);
        
        // Create layers within each segment
        for (G4int layerIndex = 0; layerIndex < numLayers; layerIndex++) {

            // Calculate inner and outer radii for this layer
            G4double layerInnerRadius = innerRadius + layerIndex * layerThickness;
            G4double absorberOuterRadius = layerInnerRadius + layerAbsorberThickness;
            G4double layerOuterRadius = layerInnerRadius + layerThickness;
            
            // Initialize layer radii arrays
            G4double layerRInner[]    = {layerInnerRadius, layerInnerRadius};
            G4double absorberROuter[] = {absorberOuterRadius, absorberOuterRadius};
            G4double layerROuter[]    = {layerOuterRadius, layerOuterRadius};
            
            // Create absorber layer
            G4Polyhedra* absorberSolid = new G4Polyhedra(baseName+"_barrel_Absorber", segmentPhiStart, segmentAngle, 1, 2, zPlanes, layerRInner, absorberROuter);
            G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid, absorberMaterial, baseName+"_barrel_Absorber_log");

            absorberLogical->SetVisAttributes(absorberVisAtt);
            
            // Place absorber layer
            new G4PVPlacement(nullptr, G4ThreeVector(), absorberLogical, baseName+"_barrel_Absorber_phys", segmentLogical, false, layerIndex, true);
            
            // Create scintillator layer
            G4Polyhedra* scintillatorSolid = new G4Polyhedra(baseName+"_barrel_Scintillator", segmentPhiStart, segmentAngle, 1, 2, zPlanes, absorberROuter, layerROuter);
            G4LogicalVolume* scintillatorLogical = new G4LogicalVolume(scintillatorSolid, scintillatorMaterial, baseName+"_barrel_Scintillator_log");

            scintillatorLogical->SetVisAttributes(scintillatorVisAtt);
            
            // Place scintillator layer
            new G4PVPlacement(nullptr, G4ThreeVector(), scintillatorLogical, baseName+"_barrel_Scintillator_phys", segmentLogical, false, layerIndex, true);
            
            // Store scintillator logical volume for later use (sensitive detector attachment, etc.)
            if (layerIndex == 0) {
                *outScintillatorVolume = scintillatorLogical;
            }

        }
    }

}

void DetectorConstruction::ConstructSamplingEndcap(G4String baseName,
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
                                                   G4Colour visColor)
{

    // Create endcaps
    G4Tubs* endcapSolid = new G4Tubs(baseName+"_endcap", 0, endcapRadius, totalThickness/2, 0, twopi);
    G4LogicalVolume* endcapLogical = new G4LogicalVolume(endcapSolid, fWorldMaterial, baseName+"_endcap_log");

    // Visualization attributes
    G4VisAttributes* endcapVisAtt = new G4VisAttributes(visColor);
    endcapVisAtt->SetVisibility(true);
    endcapLogical->SetVisAttributes(endcapVisAtt);
    
    // Place endcaps
    G4double endcapZ = endcapStart + totalThickness/2;
    new G4PVPlacement(0, G4ThreeVector(0, 0, endcapZ), endcapLogical, baseName+"_endcap_pos_phys", parentVolume, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -endcapZ), endcapLogical, baseName+"_endcap_neg_phys", parentVolume, false, 1);

    *outVolume = endcapLogical;

    G4double layerPosition = -totalThickness/2;

    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    scintillatorVisAtt->SetVisibility(false);

    G4VisAttributes* absorberVisAtt = new G4VisAttributes(visColor);
    absorberVisAtt->SetVisibility(false);

    // Create layer structure
    for (G4int layerIndex = 0; layerIndex < numLayers; layerIndex++) {

        // Calculate z distance for this layer
        G4double layerAbsorberZ = layerPosition + layerAbsorberThickness/2;

        // Create absorber layer
        G4Tubs* absorberSolid = new G4Tubs(baseName+"_endcap_Absorber", 0, endcapRadius, layerAbsorberThickness/2, 0, twopi);
        G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid, absorberMaterial, baseName+"_endcap_Absorber_log");

        absorberLogical->SetVisAttributes(absorberVisAtt);
        
        // Place absorber layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerAbsorberZ), absorberLogical, baseName+"_endcap_Absorber_phys", endcapLogical, false, layerIndex, true);
        
        layerPosition += layerAbsorberThickness;

        G4double layerScintillatorZ = layerPosition + (layerThickness - layerAbsorberThickness)/2;

        // Create scintillator layer
        G4Tubs* scintillatorSolid = new G4Tubs(baseName+"_endcap_Scintillator", 0, endcapRadius, (layerThickness - layerAbsorberThickness)/2, 0, twopi);
        G4LogicalVolume* scintillatorLogical = new G4LogicalVolume(scintillatorSolid, scintillatorMaterial, baseName+"_endcap_Scintillator_log");

        scintillatorLogical->SetVisAttributes(scintillatorVisAtt);
        
        // Place scintillator layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerScintillatorZ), scintillatorLogical, baseName+"_endcap_Scintillator_phys", endcapLogical, false, layerIndex, true);
        
        layerPosition += (layerThickness - layerAbsorberThickness);

        // Store scintillator logical volume for later use (sensitive detector attachment, etc.)
        if (layerIndex == 0) {
            *outScintillatorVolume = scintillatorLogical;
        }

    }

}

void DetectorConstruction::ConstructSDandField()
{

    // Create a uniform magnetic field
    G4ThreeVector fieldValue = G4ThreeVector(0., 0., -fMagneticFieldStrength);
    fMagneticField = new G4UniformMagField(fieldValue);

    // Create a field manager
    G4FieldManager* fieldManager = new G4FieldManager(fMagneticField);

    // Assign local magnetic field to logical volume
    fFieldLogical->SetFieldManager(fieldManager, true);

}

// Messenger methods implementation
void DetectorConstruction::SetTPCRadius(G4double radius)
{
    fTPCRadius = radius;
    G4cout << "TPC radius set to " << fTPCRadius/cm << " cm" << G4endl;

    UpdateGeometry();

}

void DetectorConstruction::SetTPCLength(G4double length)
{
    fTPCLength = length;
    G4cout << "TPC length set to " << fTPCLength/cm << " cm" << G4endl;

    UpdateGeometry();

}

void DetectorConstruction::SetECalAbsorberThickness(G4double thickness)
{
    fECalAbsorberThickness = thickness;
    G4cout << "ECal absorber thickness set to " << fECalAbsorberThickness/mm << " mm" << G4endl;

    UpdateGeometry();
}

void DetectorConstruction::SetECalScintillatorThickness(G4double thickness)
{
    fECalScintillatorThickness = thickness;
    G4cout << "ECal scintillator thickness set to " << fECalScintillatorThickness/mm << " mm" << G4endl;

    UpdateGeometry();
}

void DetectorConstruction::SetECalLayers(G4int layers)
{
    fECalLayers = layers;
    G4cout << "ECal layers set to " << fECalLayers << G4endl;

    UpdateGeometry();
}

void DetectorConstruction::SetMuIDAbsorberThickness(G4double thickness)
{
    fMuIDAbsorberThickness = thickness;
    G4cout << "MuID absorber thickness set to " << fMuIDAbsorberThickness/cm << " cm" << G4endl;

    UpdateGeometry();
}

void DetectorConstruction::SetMuIDScintillatorThickness(G4double thickness)
{
    fMuIDScintillatorThickness = thickness;
    G4cout << "MuID scintillator thickness set to " << fMuIDScintillatorThickness/cm << " cm" << G4endl;

    UpdateGeometry();
}

void DetectorConstruction::SetMuIDLayers(G4int layers)
{
    fMuIDLayers = layers;
    G4cout << "MuID layers set to " << fMuIDLayers << G4endl;

    UpdateGeometry();
}

void DetectorConstruction::DefineCommands()
{
    // Create a new messenger
    fMessenger = new DetectorMessenger(this);
}