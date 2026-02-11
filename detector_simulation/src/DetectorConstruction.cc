#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "AnalysisManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
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
  fGeometryType(kGArLike),
  fWorldLogical(nullptr),
  fFieldLogical(nullptr), fTPCLogical(nullptr),
  fECalBarrelLogical(nullptr), fECalEndcapsLogical(nullptr), fECalScintillatorLogical(nullptr),
  fMuIDLogical(nullptr), fMuIDScintillatorLogical(nullptr),
  fLArTPCLogical(nullptr),
  fWorldPhysical(nullptr),
  fFieldPhysical(nullptr), fTPCPhysical(nullptr),
  fLArTPCPhysical(nullptr),
  fMagneticField(nullptr), fMagneticFieldStrength(0.5*tesla),
  fTPCRadius(250.0*cm), fTPCLength(500.0*cm), fTPCPCBThickness(0.5*cm),
  fECalBarrelGap(20.0*cm), fECalEndcapGap(80.0*cm),
  fECalNumSides(12),
  fECalHGAbsorberThickness(0.7*mm), fECalHGScintillatorThickness(5.0*mm), fECalHGBoardThickness(1.0*mm),
  fECalBarrelHGLayers(8), fECalEndcapHGLayers(6),
  fECalLGAbsorberThickness(1.4*mm), fECalLGScintillatorThickness(10.0*mm),
  fECalBarrelLGLayers(34), fECalEndcapLGLayers(36),
  fMuIDBarrelGap(50.0*cm),
  fMuIDAbsorberThickness(10.0*cm), fMuIDScintillatorThickness(2.0*cm), fMuIDNumSides(21), fMuIDLayers(3),
  fLArNModulesX(5), fLArNModulesY(1), fLArNModulesZ(7),
  fLArModuleLength(100.0*cm), fLArModuleWidth(300.0*cm), fLArModuleDepth(100.0*cm),
  fLArModuleGap(1.0*cm), fLArInsulationThickness(1.0*mm), fLArCryostatThickness(5.0*cm),
  fLArEnableMuonWindow(true), fLArMuonWindowThickness(4.5*cm),
  fPressure(10.13*bar), fRefPressure(1.01*bar), fTemperature(294.26*kelvin),
  fGasDensity(0.001677*g/cm3), // at atm pressure (used in GArSoft -- need to check)
  //fGasDensity(0.001677*g/cm3), // at atm pressure (from some random table)
  fFoamDensity(0.2*g/cm3), // low density foam for muon window (200 kg/m³)
  fMessenger(nullptr),
  fGeometryInitialized(false)
{

    // Define command interface
    DefineCommands();

    // Define materials
    DefineMaterials();
    
    fGeometryInitialized = true;
    
}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

void DetectorConstruction::ComputeDerivedQuantities()
{
    // Define some numerical quantities that will be used extensively
    fTPCTotalLength = fTPCLength + 2*fTPCPCBThickness;
    fECalHGLayerThickness = fECalHGAbsorberThickness + fECalHGScintillatorThickness + fECalHGBoardThickness;
    fECalLGLayerThickness = fECalLGAbsorberThickness + fECalLGScintillatorThickness;
    fMuIDLayerThickness = fMuIDAbsorberThickness + fMuIDScintillatorThickness;
    fECalBarrelTotalThickness = fECalHGLayerThickness * fECalBarrelHGLayers + fECalLGLayerThickness * fECalBarrelLGLayers;
    fECalEndcapTotalThickness = fECalHGLayerThickness * fECalEndcapHGLayers + fECalLGLayerThickness * fECalEndcapLGLayers;
    fMuIDTotalThickness = fMuIDLayerThickness * fMuIDLayers;

    fLArTotalLength = fLArNModulesX * fLArModuleLength + (fLArNModulesX + 1) * fLArModuleGap;
    fLArTotalWidth  = fLArNModulesY * fLArModuleWidth  + (fLArNModulesY + 1) * fLArModuleGap;
    fLArTotalDepth  = fLArNModulesZ * fLArModuleDepth  + (fLArNModulesZ + 1) * fLArModuleGap;
}

G4bool DetectorConstruction::UpdateGeometry()
{
    // Tell G4RunManager to rebuild geometry
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    return true;
}

void DetectorConstruction::DefineMaterials()
{
    // Define materials for different detector components
	G4NistManager *nistManager = G4NistManager::Instance();

    // Get needed elements
    G4double z, fractionmass;
	G4int nel, natoms;
    G4Element* H  = nistManager->FindOrBuildElement("H");
    G4Element* C  = nistManager->FindOrBuildElement("C");
    G4Element* N  = nistManager->FindOrBuildElement("N");
    G4Element* O  = nistManager->FindOrBuildElement("O");
    G4Element* Ar = nistManager->FindOrBuildElement("Ar");
    G4Element* B  = nistManager->FindOrBuildElement("B");
    G4Element* Na = nistManager->FindOrBuildElement("Na");
    G4Element* Al = nistManager->FindOrBuildElement("Al");
    G4Element* Si = nistManager->FindOrBuildElement("Si");
    G4Element* K  = nistManager->FindOrBuildElement("K");
    G4Element* Ti = nistManager->FindOrBuildElement("Ti");
    G4Element* Zn = nistManager->FindOrBuildElement("Zn");

    // Define epoxy material
    G4Material* Epoxy = new G4Material("Epoxy", 1.25*g/cm3, 3);
    Epoxy->AddElement(C, 6);
    Epoxy->AddElement(H, 6);
    Epoxy->AddElement(O, 1);

    // Define glass material
    G4Material* Glass = new G4Material("Glass", 2.7*g/cm3, 8);
    Glass->AddElement(Si, 0.2743);  // Silicon: 27.43%
    Glass->AddElement(B,  0.0166);  // Boron: 1.66%
    Glass->AddElement(Al, 0.0207);  // Aluminum: 2.07%
    Glass->AddElement(Na, 0.0449);  // Sodium: 4.49%
    Glass->AddElement(K,  0.0821);  // Potassium: 8.21%
    Glass->AddElement(Zn, 0.0882);  // Zinc: 8.82%
    Glass->AddElement(Ti, 0.0292);  // Titanium: 2.92%
    Glass->AddElement(O,  0.444);   // Oxygen: 44.4%

    G4Material* fr4 = new G4Material("FR4", 1.85*g/cm3, 2);
    fr4->AddMaterial(Epoxy, 0.6);   // Epoxy: 60% (prev 0.206)
    fr4->AddMaterial(Glass, 0.4);   // Glass: 40% (prev 0.794)
	
    // World material definition -- just air
    fWorldMaterial = nistManager->FindOrBuildMaterial("G4_AIR");

    // Define gas material
    fGArTPCMaterial = new G4Material("GasMixture", fGasDensity*fPressure/fRefPressure, nel = 3, kStateGas, fTemperature, fPressure);
    // Add elements to material
    fGArTPCMaterial->AddElement(H,  fractionmass = 0.011);
    fGArTPCMaterial->AddElement(C,  fractionmass = 0.032);
    fGArTPCMaterial->AddElement(Ar, fractionmass = 0.957);

    // Materials for drift chamber
    fTPCPCBMaterial = fr4; 

    // Materials for ECal
    fECalAbsorberMaterial = nistManager->FindOrBuildMaterial("G4_Pb");
    fECalScintillatorMaterial = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    fECalPCBMaterial = fr4;
    
    // Materials for MuID
    fMuIDAbsorberMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    fMuIDScintillatorMaterial = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // Materials for LAr TPC
    fLArTPCMaterial = nistManager->FindOrBuildMaterial("G4_lAr"); // use pre-defined LAr material
    fLArCryostatMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    fLArInsulationMaterial = nistManager->FindOrBuildMaterial("G4_Al"); // Al for module optical isolation

    // Material for low-density foam muon window
    // Typical composition for polyurethane foam: (C3H8N2O)n
    fLArMuonWindowMaterial = new G4Material("FiberglassFoam", fFoamDensity, nel = 4, kStateSolid);
    fLArMuonWindowMaterial->AddElement(C, 3);
    fLArMuonWindowMaterial->AddElement(H, 8);
    fLArMuonWindowMaterial->AddElement(N, 2);
    fLArMuonWindowMaterial->AddElement(O, 1);
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    
    // Recalculate derived quantities
    ComputeDerivedQuantities();

    // Construct world volume
    fWorldPhysical = ConstructWorld();

    // Build detector based on selected geometry type
    switch(fGeometryType) {
        case kGArLike:
            ConstructGArDetector();
            break;
        case kLArLike:
            ConstructLArDetector();
            break;
        default:
            G4cerr << "Unknown geometry type!" << G4endl;
            break;
    }

    return fWorldPhysical;
}

// Add Geometry information to output file
void DetectorConstruction::RecordGeometry()
{
    auto analysisManager = AnalysisManager::GetInstance();

    // Book geometry tree
    analysisManager->BookGeometry();

    // Fill geometry parameters
    GeometryInfo geoInfo;

    // Set geometry type
    geoInfo.geometry_type = static_cast<int>(fGeometryType);

    // GAr TPC parameters
    geoInfo.gar_tpc_radius = fTPCRadius/cm;
    geoInfo.gar_tpc_length = fTPCLength/cm;
    geoInfo.gar_magnetic_field = fMagneticFieldStrength/tesla;
    geoInfo.gar_pressure = fPressure/bar;

    // ECal parameters
    geoInfo.ecal_barrel_gap = fECalBarrelGap/cm;
    geoInfo.ecal_endcap_gap = fECalEndcapGap/cm;
    geoInfo.ecal_num_sides = fECalNumSides;
    geoInfo.ecal_hg_absorber_thickness = fECalHGAbsorberThickness/cm;
    geoInfo.ecal_hg_scintillator_thickness = fECalHGScintillatorThickness/cm;
    geoInfo.ecal_hg_board_thickness = fECalHGBoardThickness/cm;
    geoInfo.ecal_barrel_hg_layers = fECalBarrelHGLayers;
    geoInfo.ecal_endcap_hg_layers = fECalEndcapHGLayers;
    geoInfo.ecal_lg_absorber_thickness = fECalLGAbsorberThickness/cm;
    geoInfo.ecal_lg_scintillator_thickness = fECalLGScintillatorThickness/cm;
    geoInfo.ecal_barrel_lg_layers = fECalBarrelLGLayers;
    geoInfo.ecal_endcap_lg_layers = fECalEndcapLGLayers;

    // MuID parameters
    geoInfo.muid_barrel_gap = fMuIDBarrelGap/cm;
    geoInfo.muid_absorber_thickness = fMuIDAbsorberThickness/cm;
    geoInfo.muid_scintillator_thickness = fMuIDScintillatorThickness/cm;
    geoInfo.muid_num_sides = fMuIDNumSides;
    geoInfo.muid_layers = fMuIDLayers;

    // LAr TPC parameters
    geoInfo.lar_n_modules_x = fLArNModulesX;
    geoInfo.lar_n_modules_y = fLArNModulesY;
    geoInfo.lar_n_modules_z = fLArNModulesZ;
    geoInfo.lar_module_length = fLArModuleLength/cm;
    geoInfo.lar_module_width = fLArModuleWidth/cm;
    geoInfo.lar_module_depth = fLArModuleDepth/cm;
    geoInfo.lar_module_gap = fLArModuleGap/cm;
    geoInfo.lar_insulation_thickness = fLArInsulationThickness/cm;
    geoInfo.lar_cryostat_thickness = fLArCryostatThickness/cm;
    geoInfo.lar_enable_muon_window = fLArEnableMuonWindow;
    geoInfo.lar_muon_window_thickness = fLArMuonWindowThickness/cm;

    analysisManager->FillGeometryInfo(geoInfo);
}

void DetectorConstruction::ConstructGArDetector()
{
    // Construct detector components
    ConstructFieldEnclosure();
    ConstructTPC();
    ConstructECal();
    ConstructMuID();

    G4UserLimits* limitsTPCSteps = new G4UserLimits(5.0*mm);
    G4UserLimits* limitsCaloSteps = new G4UserLimits(5.0*mm);

    fTPCLogical->SetUserLimits(limitsTPCSteps);
    fECalBarrelLogical->SetUserLimits(limitsCaloSteps);
    fECalEndcapsLogical->SetUserLimits(limitsCaloSteps);
    fMuIDLogical->SetUserLimits(limitsCaloSteps);
}

void DetectorConstruction::ConstructLArDetector()
{
    // Check that module boxes don't overlap
    if (fLArInsulationThickness*2 > fLArModuleGap) {
        G4cerr << "LAr TPC modules will overlap!" << G4endl;
    }

    // Create cryostat solid
    G4Box* cryostatOuterSolid = new G4Box("Cryostat_outer",
                                          (fLArTotalLength + 2*fLArCryostatThickness) / 2,
                                          (fLArTotalWidth  + 2*fLArCryostatThickness) / 2,
                                          (fLArTotalDepth  + 2*fLArCryostatThickness) / 2);

    G4Box* cryostatInnerSolid = new G4Box("Cryostat_inner", fLArTotalLength/2, fLArTotalWidth/2, fLArTotalDepth/2);

    G4SubtractionSolid* cryostatSolid = new G4SubtractionSolid("Cryostat", cryostatOuterSolid, cryostatInnerSolid);

    // If muon window enabled cut hole in cryostat solid
    if (fLArEnableMuonWindow) {
        // Check that muon window is thiner than cryostat
        if (fLArMuonWindowThickness > fLArCryostatThickness) {
            G4cerr << "Muon window thicker than cryostat!" << G4endl;
        }
        
        // Remove boc in +X face to make room for muon window
        G4Box* cutoutBox = new G4Box("Cutout",
                                     fLArMuonWindowThickness/2,
                                     fLArTotalWidth/2,
                                     fLArTotalDepth/2);
        
        G4ThreeVector cutoutPos(fLArTotalLength/2 + fLArCryostatThickness/2, 0, 0);
        
        cryostatSolid = new G4SubtractionSolid("Cryostat",
                                               cryostatSolid,
                                               cutoutBox,
                                               0,
                                               cutoutPos);
    }

    // Place cryostat in world
    G4LogicalVolume* cryostatLogical = new G4LogicalVolume(cryostatSolid, fLArCryostatMaterial, "Cryostat_log");
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), cryostatLogical, "Cryostat_phys", fWorldLogical, false, 0, true);

    // Set visualization attributes
    G4VisAttributes* cryostatVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.1));
    cryostatVisAtt->SetVisibility(true);
    cryostatLogical->SetVisAttributes(cryostatVisAtt);

    // If enabled create and place muon window
    if (fLArEnableMuonWindow) {
        
        G4Box* muonWindowSolid = new G4Box("MuonWindow",
                                           fLArMuonWindowThickness/2,
                                           fLArTotalWidth/2,
                                           fLArTotalDepth/2);
        
        G4LogicalVolume* muonWindowLogical = new G4LogicalVolume(muonWindowSolid,
                                                                 fLArMuonWindowMaterial,
                                                                 "MuonWindow_log");
        
        G4ThreeVector windowPos(fLArTotalLength/2 + fLArCryostatThickness/2, 0, 0);
        
        new G4PVPlacement(0, windowPos, muonWindowLogical, "MuonWindow_phys", fWorldLogical, false, 0, true);

        // Set visualization attributes
        G4VisAttributes* muonWindowVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.3));
        muonWindowVisAtt->SetVisibility(true);
        muonWindowLogical->SetVisAttributes(muonWindowVisAtt);
    }

    // Create mother volume for all TPC modules
    G4Box* motherSolid = new G4Box("TPC", fLArTotalLength/2, fLArTotalWidth/2, fLArTotalDepth/2);
    G4LogicalVolume* motherLogical = new G4LogicalVolume(motherSolid, fWorldMaterial, "TPC_log");
    // Place mother volume inside cryostat
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), motherLogical, "TPC_phys", fWorldLogical, false, 0, true);

    // Set visualization attributes
    G4VisAttributes* motherVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    motherVisAtt->SetVisibility(false);
    motherLogical->SetVisAttributes(motherVisAtt);

    // Create an insulation box for TPC modules (will be replicated)
    G4Box* insulationOuterSolid = new G4Box("Insulation_outer",
                                            (fLArModuleLength + 2*fLArInsulationThickness) / 2,
                                            (fLArModuleWidth  + 2*fLArInsulationThickness) / 2,
                                            (fLArModuleDepth  + 2*fLArInsulationThickness) / 2);

    G4Box* insulationInnerSolid = new G4Box("Insulation_inner", fLArModuleLength/2, fLArModuleWidth/2, fLArModuleDepth/2);

    G4SubtractionSolid* insulationSolid = new G4SubtractionSolid("Insulation", insulationOuterSolid, insulationInnerSolid);
    G4LogicalVolume* insulationLogical = new G4LogicalVolume(insulationSolid, fLArInsulationMaterial, "Insulation_log");

    // Create a single TPC module (will be replicated)
    G4Box* moduleSolid = new G4Box("TPC_module", fLArModuleLength/2, fLArModuleWidth/2, fLArModuleDepth/2);
    fLArTPCLogical = new G4LogicalVolume(moduleSolid, fLArTPCMaterial, "TPC_module_log");

    // Place modules and insulation boxes
    G4int moduleID = 0;
    for (G4int ix = 0; ix < fLArNModulesX; ix++) {
        for (G4int iy = 0; iy < fLArNModulesY; iy++) {
            for (G4int iz = 0; iz < fLArNModulesZ; iz++) {

                // Calculate module position
                G4double xPos = -fLArTotalLength/2 + fLArModuleGap + fLArModuleLength/2 + ix * (fLArModuleLength + fLArModuleGap);
                G4double yPos = 0; // Centered in Y
                G4double zPos = -fLArTotalDepth/2  + fLArModuleGap + fLArModuleDepth/2  + iz * (fLArModuleDepth  + fLArModuleGap);

                // Place insulation box
                new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos),
                                  insulationLogical, "Insulation_phys",
                                  motherLogical, false, moduleID, true);

                // Place module inside insulation box
                new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos),
                                  fLArTPCLogical, "TPC_module_phys",
                                  motherLogical, false, moduleID, true);

                moduleID++;

            }
        }
    }

    // Set visualization attributes
    G4VisAttributes* insulationVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.0));
    insulationVisAtt->SetVisibility(false);
    insulationLogical->SetVisAttributes(insulationVisAtt);
    
    G4VisAttributes* moduleVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    moduleVisAtt->SetVisibility(true);
    fLArTPCLogical->SetVisAttributes(moduleVisAtt);
    
    // Set user limits for tracking
    G4UserLimits* limitsTPCSteps = new G4UserLimits(1.0*mm);
    fLArTPCLogical->SetUserLimits(limitsTPCSteps);
}

G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{
    G4double worldSizeX, worldSizeY, worldSizeZ;

    if (fGeometryType == kGArLike) {
        // Calculate world size based on GAr detector dimensions
        worldSizeX = 2.5 * (fTPCRadius + fECalBarrelTotalThickness + fMuIDTotalThickness);
        worldSizeY = worldSizeX;
        worldSizeZ = 2.5 * (fTPCTotalLength/2 + fECalEndcapTotalThickness + fMuIDTotalThickness);
    } else if (fGeometryType == kLArLike) {
        // Calculate world size based on LAr detector dimensions
        worldSizeX = 2.0 * fLArTotalLength;
        worldSizeY = 2.0 * fLArTotalWidth;
        worldSizeZ = 2.0 * fLArTotalDepth;
    } else {
        G4cerr << "Unknown geometry type!" << G4endl;
    }
    
    // Create world box
    G4Box* worldSolid = new G4Box("World", worldSizeX, worldSizeY, worldSizeZ);
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
    G4double fieldEnclosureLength = fTPCTotalLength + 2*fECalEndcapGap + 2*fECalEndcapTotalThickness;
    G4double fieldEnclosureRadius = fTPCRadius + fECalBarrelGap + fECalBarrelTotalThickness + fMuIDBarrelGap;

    // Create field enclosure cylinder
    G4Tubs* fieldSolid = new G4Tubs("FieldEnclosure", 0, fieldEnclosureRadius, fieldEnclosureLength/2, 0, twopi);
    fFieldLogical = new G4LogicalVolume(fieldSolid, fGArTPCMaterial, "FieldEnclosure_log");

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
    G4Tubs* tpcSolid = new G4Tubs("GArTPC", 0, fTPCRadius, fTPCTotalLength/2, 0, twopi);
    fTPCLogical = new G4LogicalVolume(tpcSolid, fGArTPCMaterial, "TPC_log");

    // Place TPC in world
    fTPCPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fTPCLogical, "TPC_phys", fFieldLogical, false, 0);

    // Create and place TPC gas volume
    G4Tubs* tpcGasSolid = new G4Tubs("TPCGas", 0, fTPCRadius, fTPCLength/2, 0, twopi);
    G4LogicalVolume* fTPCGasLogical = new G4LogicalVolume(tpcGasSolid, fGArTPCMaterial, "TPCGas_log");
    G4VisAttributes* tpcGasVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    tpcGasVisAtt->SetVisibility(true);
    fTPCGasLogical->SetVisAttributes(tpcGasVisAtt);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fTPCGasLogical, "TPCGas_phys", fTPCLogical, false, 0, true);

    // Create PCB volume
    G4Tubs* tpcPCBSolid = new G4Tubs("TPCPCB", 0, fTPCRadius, fTPCPCBThickness/2, 0, twopi);
    G4LogicalVolume* tpcPCBLogical = new G4LogicalVolume(tpcPCBSolid, fTPCPCBMaterial, "TPCPCB_log");
    G4VisAttributes* tpcPCBVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.3));
    tpcPCBVisAtt->SetVisibility(true);
    tpcPCBLogical->SetVisAttributes(tpcPCBVisAtt);

    // Rotation matrix to flip the negative PCB
    G4RotationMatrix* flip = new G4RotationMatrix();
    flip->rotateY(180*deg);
    
    // Place PCB volume
    G4double pcbZ = fTPCLength/2 + fTPCPCBThickness/2;
    new G4PVPlacement(0, G4ThreeVector(0, 0, pcbZ), tpcPCBLogical, "TPCPCB_pos_phys", fTPCLogical, false, 0, true);
    new G4PVPlacement(flip, G4ThreeVector(0, 0, -pcbZ), tpcPCBLogical, "TPCPCB_neg_phys", fTPCLogical, false, 0, true);
}

void DetectorConstruction::ConstructECal()
{

    G4double ecalBarrelTotalLength = fTPCTotalLength + 2*fECalEndcapGap + 2*fECalEndcapTotalThickness;
    G4double ecalBarrelInnerDistance = fTPCRadius + fECalBarrelGap;

    ConstructSamplingBarrel("ECal",
                            ecalBarrelTotalLength, ecalBarrelInnerDistance, fECalNumSides, fECalBarrelTotalThickness,
                            fECalHGLayerThickness, fECalHGAbsorberThickness, fECalHGScintillatorThickness,
                            fECalLGLayerThickness, fECalLGAbsorberThickness,
                            fECalBarrelHGLayers, fECalBarrelLGLayers,
                            fECalAbsorberMaterial, fECalScintillatorMaterial, fECalPCBMaterial,
                            fFieldLogical,
                            &fECalBarrelLogical, &fECalScintillatorLogical,
                            G4Colour(0.0, 1.0, 0.0, 0.3));

    G4double ecalEndcapRadius = ecalBarrelInnerDistance;
    G4double ecalEndcapStart = fTPCTotalLength/2 + fECalEndcapGap;

    ConstructSamplingEndcap("ECal",
                            ecalEndcapStart, ecalEndcapRadius, fECalEndcapTotalThickness,
                            fECalHGLayerThickness, fECalHGAbsorberThickness, fECalHGScintillatorThickness,
                            fECalLGLayerThickness, fECalLGAbsorberThickness,
                            fECalEndcapHGLayers, fECalEndcapLGLayers,
                            fECalAbsorberMaterial, fECalScintillatorMaterial, fECalPCBMaterial,
                            fFieldLogical,
                            &fECalEndcapsLogical, &fECalScintillatorLogical,
                            G4Colour(0.0, 1.0, 0.0, 0.3));
}

void DetectorConstruction::ConstructMuID()
{
    G4double muidBarrelTotalLength = fTPCTotalLength + 2*fECalEndcapGap + 2*fECalEndcapTotalThickness;
    G4double muidBarrelInnerDistance = fTPCRadius + fECalBarrelGap + fECalBarrelTotalThickness + fMuIDBarrelGap;

    ConstructSamplingBarrel("MuID",
                            muidBarrelTotalLength, muidBarrelInnerDistance, fMuIDNumSides, fMuIDTotalThickness,
                            0.0, 0.0, 0.0,
                            fMuIDLayerThickness, fMuIDAbsorberThickness,
                            0, fMuIDLayers,
                            fMuIDAbsorberMaterial, fMuIDScintillatorMaterial, fECalPCBMaterial,
                            fWorldLogical,
                            &fMuIDLogical, &fMuIDScintillatorLogical,
                            G4Colour(1.0, 0.0, 0.0, 0.3));
}

void DetectorConstruction::ConstructSamplingBarrel(G4String baseName,
                                                   G4double barrelLength,
                                                   G4double barrelInnerDistance,
                                                   G4int numSides,
                                                   G4double totalThickness,
                                                   G4double layerHGThickness,
                                                   G4double layerHGAbsorberThickness,
                                                   G4double layerHGScintillatorThickness,
                                                   G4double layerLGThickness,
                                                   G4double layerLGAbsorberThickness,
                                                   G4int numHGLayers,
                                                   G4int numLGLayers,
                                                   G4Material* absorberMaterial,
                                                   G4Material* scintillatorMaterial,
                                                   G4Material* boardMaterial,
                                                   G4LogicalVolume* parentVolume,
                                                   G4LogicalVolume** outVolume,
                                                   G4LogicalVolume** outScintillatorVolume,
                                                   G4Colour visColor)
{
    // Define barrel dimensions
    G4double cosineAngle = std::cos(pi / numSides);
    G4double innerApothem = barrelInnerDistance;
    G4double outerApothem = innerApothem + totalThickness;
    G4double halfLength = barrelLength / 2;

    // Create barrel
    G4double phiStart = 0;
    G4double phiTotal = twopi;
    G4double zPlanes[] = {-halfLength, halfLength};
    G4double rInner[]  = {innerApothem, innerApothem};
    G4double rOuter[]  = {outerApothem, outerApothem};
    
    G4Polyhedra* barrelSolid = new G4Polyhedra(baseName+"_barrel", phiStart, phiTotal, numSides, 2, zPlanes, rInner, rOuter);
    G4LogicalVolume* barrelLogical = new G4LogicalVolume(barrelSolid, fWorldMaterial, baseName+"_barrel_log");
    
    // Visualization attributes
    G4VisAttributes* barrelVisAtt = new G4VisAttributes(visColor);
    barrelVisAtt->SetVisibility(true);
    barrelLogical->SetVisAttributes(barrelVisAtt);
    
    // Place barrel
    G4double rotAngle = 90.*deg - 180.*deg / numSides;
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateZ(rotAngle);
    new G4PVPlacement(rotation, G4ThreeVector(0, 0, 0), barrelLogical, baseName+"_barrel_phys", parentVolume, false, 0);

    *outVolume = barrelLogical;

    // Calculate segment angle
    G4double segmentAngle = phiTotal / numSides;

    // Set visualisation attributes
    G4VisAttributes* segmentVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    segmentVisAtt->SetVisibility(false);

    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    scintillatorVisAtt->SetVisibility(false);

    G4VisAttributes* absorberVisAtt = new G4VisAttributes(visColor);
    absorberVisAtt->SetVisibility(false);

    G4VisAttributes* boardVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    boardVisAtt->SetVisibility(false);
    
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
        
        // Create HG layers within each segment
        for (G4int layerIndex = 0; layerIndex < numHGLayers; layerIndex++) {

            // Define unique copy number in terms of segment and layer index
            G4int layerCopyNo = 1000 * segmentIndex + layerIndex;

            // Calculate inner and outer radii for this layer
            G4double layerInnerApothem = innerApothem + layerIndex * layerHGThickness;
            G4double absorberOuterApothem = layerInnerApothem + layerHGAbsorberThickness;
            G4double scintillatorOuterApothem = absorberOuterApothem + layerHGScintillatorThickness;
            G4double layerOuterApothem = layerInnerApothem + layerHGThickness;

            // Initialize layer radii arrays
            G4double layerRInner[]    = {layerInnerApothem, layerInnerApothem};
            G4double absorberROuter[] = {absorberOuterApothem, absorberOuterApothem};
            G4double scintillatorROuter[] = {scintillatorOuterApothem, scintillatorOuterApothem};
            G4double layerROuter[]    = {layerOuterApothem, layerOuterApothem};
            
            // Create absorber layer
            G4Polyhedra* absorberSolid = new G4Polyhedra(baseName+"_barrel_HG_Absorber", segmentPhiStart, segmentAngle, 1, 2, zPlanes, layerRInner, absorberROuter);
            G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid, absorberMaterial, baseName+"_barrel_HG_Absorber_log");

            absorberLogical->SetVisAttributes(absorberVisAtt);
            
            // Place absorber layer
            new G4PVPlacement(nullptr, G4ThreeVector(), absorberLogical, baseName+"_barrel_HG_Absorber_phys", segmentLogical, false, layerCopyNo, true);
            
            // Create scintillator layer
            G4Polyhedra* scintillatorSolid = new G4Polyhedra(baseName+"_barrel_HG_Scintillator", segmentPhiStart, segmentAngle, 1, 2, zPlanes, absorberROuter, scintillatorROuter);
            G4LogicalVolume* scintillatorLogical = new G4LogicalVolume(scintillatorSolid, scintillatorMaterial, baseName+"_barrel_HG_Scintillator_log");

            scintillatorLogical->SetVisAttributes(scintillatorVisAtt);
            
            // Place scintillator layer
            new G4PVPlacement(nullptr, G4ThreeVector(), scintillatorLogical, baseName+"_barrel_HG_Scintillator_phys", segmentLogical, false, layerCopyNo, true);
            
            // Store scintillator logical volume for later use (sensitive detector attachment, etc.)
            if (layerIndex == 0) {
                *outScintillatorVolume = scintillatorLogical;
            }

            // Create PCB layer
            G4Polyhedra* boardSolid = new G4Polyhedra(baseName+"_barrel_HG_PCB", segmentPhiStart, segmentAngle, 1, 2, zPlanes, scintillatorROuter, layerROuter);
            G4LogicalVolume* boardLogical = new G4LogicalVolume(boardSolid, boardMaterial, baseName+"_barrel_HG_PCB_log");

            boardLogical->SetVisAttributes(boardVisAtt);
            
            // Place PCB layer
            new G4PVPlacement(nullptr, G4ThreeVector(), boardLogical, baseName+"_barrel_HG_PCB_phys", segmentLogical, false, layerCopyNo, true);

        }

        // LG layers start at an offset inner radius
        G4double innerLGApothem = innerApothem + numHGLayers * layerHGThickness;

        // Create LG layers within each segment
        for (G4int layerIndex = 0; layerIndex < numLGLayers; layerIndex++) {

            // Define unique copy number in terms of segment and layer index
            G4int layerCopyNo = 1000 * segmentIndex + layerIndex + numHGLayers;

            // Calculate inner and outer radii for this layer
            G4double layerInnerApothem = innerLGApothem + layerIndex * layerLGThickness;
            G4double absorberOuterApothem = layerInnerApothem + layerLGAbsorberThickness;
            G4double layerOuterApothem = layerInnerApothem + layerLGThickness;

            // Initialize layer radii arrays
            G4double layerRInner[]    = {layerInnerApothem, layerInnerApothem};
            G4double absorberROuter[] = {absorberOuterApothem, absorberOuterApothem};
            G4double layerROuter[]    = {layerOuterApothem, layerOuterApothem};
            
            // Create absorber layer
            G4Polyhedra* absorberSolid = new G4Polyhedra(baseName+"_barrel_LG_Absorber", segmentPhiStart, segmentAngle, 1, 2, zPlanes, layerRInner, absorberROuter);
            G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid, absorberMaterial, baseName+"_barrel_LG_Absorber_log");

            absorberLogical->SetVisAttributes(absorberVisAtt);
            
            // Place absorber layer
            new G4PVPlacement(nullptr, G4ThreeVector(), absorberLogical, baseName+"_barrel_LG_Absorber_phys", segmentLogical, false, layerCopyNo, true);
            
            // Create scintillator layer
            G4Polyhedra* scintillatorSolid = new G4Polyhedra(baseName+"_barrel_LG_Scintillator", segmentPhiStart, segmentAngle, 1, 2, zPlanes, absorberROuter, layerROuter);
            G4LogicalVolume* scintillatorLogical = new G4LogicalVolume(scintillatorSolid, scintillatorMaterial, baseName+"_barrel_LG_Scintillator_log");

            scintillatorLogical->SetVisAttributes(scintillatorVisAtt);
            
            // Place scintillator layer
            new G4PVPlacement(nullptr, G4ThreeVector(), scintillatorLogical, baseName+"_barrel_LG_Scintillator_phys", segmentLogical, false, layerCopyNo, true);

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
                                                   G4double layerHGThickness,
                                                   G4double layerHGAbsorberThickness,
                                                   G4double layerHGScintillatorThickness,
                                                   G4double layerLGThickness,
                                                   G4double layerLGAbsorberThickness,
                                                   G4int numHGLayers,
                                                   G4int numLGLayers,
                                                   G4Material* absorberMaterial,
                                                   G4Material* scintillatorMaterial,
                                                   G4Material* boardMaterial,
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

    // Rotation matrix to flip the negative endcap
    G4RotationMatrix* flip = new G4RotationMatrix();
    flip->rotateY(180*deg);
    
    // Place endcaps
    G4double endcapZ = endcapStart + totalThickness/2;
    new G4PVPlacement(0, G4ThreeVector(0, 0, endcapZ), endcapLogical, baseName+"_endcap_pos_phys", parentVolume, false, 0);
    new G4PVPlacement(flip, G4ThreeVector(0, 0, -endcapZ), endcapLogical, baseName+"_endcap_neg_phys", parentVolume, false, 1);

    *outVolume = endcapLogical;

    G4double layerPosition = -totalThickness/2;

    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    scintillatorVisAtt->SetVisibility(false);

    G4VisAttributes* absorberVisAtt = new G4VisAttributes(visColor);
    absorberVisAtt->SetVisibility(false);

    G4VisAttributes* boardVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    boardVisAtt->SetVisibility(false);

    // Create HG layer structure
    for (G4int layerIndex = 0; layerIndex < numHGLayers; layerIndex++) {

        // Calculate z distance for this layer
        G4double layerAbsorberZ = layerPosition + layerHGAbsorberThickness/2;

        // Create absorber layer
        G4Tubs* absorberSolid = new G4Tubs(baseName+"_endcap_HG_Absorber", 0, endcapRadius, layerHGAbsorberThickness/2, 0, twopi);
        G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid, absorberMaterial, baseName+"_endcap_HG_Absorber_log");

        absorberLogical->SetVisAttributes(absorberVisAtt);
        
        // Place absorber layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerAbsorberZ), absorberLogical, baseName+"_endcap_HG_Absorber_phys", endcapLogical, false, layerIndex, true);
        
        layerPosition += layerHGAbsorberThickness;

        G4double layerScintillatorZ = layerPosition + layerHGScintillatorThickness/2;

        // Create scintillator layer
        G4Tubs* scintillatorSolid = new G4Tubs(baseName+"_endcap_HG_Scintillator", 0, endcapRadius, layerHGScintillatorThickness/2, 0, twopi);
        G4LogicalVolume* scintillatorLogical = new G4LogicalVolume(scintillatorSolid, scintillatorMaterial, baseName+"_endcap_HG_Scintillator_log");

        scintillatorLogical->SetVisAttributes(scintillatorVisAtt);
        
        // Place scintillator layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerScintillatorZ), scintillatorLogical, baseName+"_endcap_HG_Scintillator_phys", endcapLogical, false, layerIndex, true);
        
        layerPosition += layerHGScintillatorThickness;

        // Store scintillator logical volume for later use (sensitive detector attachment, etc.)
        if (layerIndex == 0) {
            *outScintillatorVolume = scintillatorLogical;
        }

        G4double layerBoardZ = layerPosition + (layerHGThickness - layerHGAbsorberThickness - layerHGScintillatorThickness)/2;

        // Create PCB layer
        G4Tubs* boardSolid = new G4Tubs(baseName+"_endcap_HG_PCB", 0, endcapRadius, (layerHGThickness - layerHGAbsorberThickness - layerHGScintillatorThickness)/2, 0, twopi);
        G4LogicalVolume* boardLogical = new G4LogicalVolume(boardSolid, boardMaterial, baseName+"_endcap_HG_PCB_log");

        boardLogical->SetVisAttributes(boardVisAtt);
        
        // Place PCB layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerBoardZ), boardLogical, baseName+"_endcap_HG_PCB_phys", endcapLogical, false, layerIndex, true);
        
        layerPosition += (layerHGThickness - layerHGAbsorberThickness - layerHGScintillatorThickness);

    }

    // Create LG layer structure
    for (G4int layerIndex = 0; layerIndex < numLGLayers; layerIndex++) {

        // Calculate z distance for this layer
        G4double layerAbsorberZ = layerPosition + layerLGAbsorberThickness/2;

        // Create absorber layer
        G4Tubs* absorberSolid = new G4Tubs(baseName+"_endcap_LG_Absorber", 0, endcapRadius, layerLGAbsorberThickness/2, 0, twopi);
        G4LogicalVolume* absorberLogical = new G4LogicalVolume(absorberSolid, absorberMaterial, baseName+"_endcap_LG_Absorber_log");

        absorberLogical->SetVisAttributes(absorberVisAtt);
        
        // Place absorber layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerAbsorberZ), absorberLogical, baseName+"_endcap_LG_Absorber_phys", endcapLogical, false, layerIndex + numHGLayers, true);
        
        layerPosition += layerLGAbsorberThickness;

        G4double layerScintillatorZ = layerPosition + (layerLGThickness - layerLGAbsorberThickness)/2;

        // Create scintillator layer
        G4Tubs* scintillatorSolid = new G4Tubs(baseName+"_endcap_LG_Scintillator", 0, endcapRadius, (layerLGThickness - layerLGAbsorberThickness)/2, 0, twopi);
        G4LogicalVolume* scintillatorLogical = new G4LogicalVolume(scintillatorSolid, scintillatorMaterial, baseName+"_endcap_LG_Scintillator_log");

        scintillatorLogical->SetVisAttributes(scintillatorVisAtt);
        
        // Place scintillator layer
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, layerScintillatorZ), scintillatorLogical, baseName+"_endcap_LG_Scintillator_phys", endcapLogical, false, layerIndex + numHGLayers, true);
        
        layerPosition += (layerLGThickness - layerLGAbsorberThickness);

        // Store scintillator logical volume for later use (sensitive detector attachment, etc.)
        if (layerIndex == 0) {
            *outScintillatorVolume = scintillatorLogical;
        }

    }

}

void DetectorConstruction::ConstructSDandField()
{
    if (fGeometryType == kGArLike) {
        // Create a uniform magnetic field
        G4ThreeVector fieldValue = G4ThreeVector(0., 0., -fMagneticFieldStrength);
        fMagneticField = new G4UniformMagField(fieldValue);

        // Create a field manager
        G4FieldManager* fieldManager = new G4FieldManager(fMagneticField);

        // Assign local magnetic field to logical volume
        fFieldLogical->SetFieldManager(fieldManager, true);
    }
    // No magnetic field for LAr TPC
}

// Messenger methods implementation
void DetectorConstruction::SetGeometryType(G4String type)
{
    if (type == "gar") {
        fGeometryType = kGArLike;
        G4cout << "Geometry set to ND-GAr-like detector (magentised GAr TPC + ECal + MuID)" << G4endl;
    } else if (type == "lar") {
        fGeometryType = kLArLike;
        G4cout << "Geometry set to ND-LAr-like detector" << G4endl;
    } else {
        G4cerr << "Unknown geometry type: " << type << G4endl;
        G4cerr << "Valid options are: 'gar' or 'lar'" << G4endl;
        return;
    }

    if (fGeometryInitialized) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetMagneticFieldStrength(G4double bfield)
{
    fMagneticFieldStrength = bfield;
    G4cout << "Magnetic field strength set to " << fMagneticFieldStrength/tesla << " T" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetTPCRadius(G4double radius)
{
    fTPCRadius = radius;
    G4cout << "TPC radius set to " << fTPCRadius/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetTPCLength(G4double length)
{
    fTPCLength = length;
    G4cout << "TPC length set to " << fTPCLength/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalAbsorberMaterial(G4String material)
{
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (!mat) {
        G4cerr << "Unknown ECal absorber material: " << material << G4endl;
        G4cerr << "Use a valid NIST name" << G4endl;   
        return;
    }
    fECalAbsorberMaterial = mat;
    G4cout << "ECal absorber material set to " << material << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalScintillatorMaterial(G4String material)
{
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (!mat) {
        G4cerr << "Unknown ECal scintillator material: " << material << G4endl;
        G4cerr << "Use a valid NIST name" << G4endl;   
        return;
    }
    fECalScintillatorMaterial = mat;
    G4cout << "ECal scintillator material set to " << material << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalHGAbsorberThickness(G4double thickness)
{
    fECalHGAbsorberThickness = thickness;
    G4cout << "ECal HG absorber thickness set to " << fECalHGAbsorberThickness/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalHGScintillatorThickness(G4double thickness)
{
    fECalHGScintillatorThickness = thickness;
    G4cout << "ECal HG scintillator thickness set to " << fECalHGScintillatorThickness/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalHGBoardThickness(G4double thickness)
{
    fECalHGBoardThickness = thickness;
    G4cout << "ECal HG PCB thickness set to " << fECalHGBoardThickness/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalLGAbsorberThickness(G4double thickness)
{
    fECalLGAbsorberThickness = thickness;
    G4cout << "ECal LG absorber thickness set to " << fECalLGAbsorberThickness/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalLGScintillatorThickness(G4double thickness)
{
    fECalLGScintillatorThickness = thickness;
    G4cout << "ECal LG scintillator thickness set to " << fECalLGScintillatorThickness/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalBarrelHGLayers(G4int layers)
{
    fECalBarrelHGLayers = layers;
    G4cout << "ECal barrel HG layers set to " << fECalBarrelHGLayers << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalBarrelLGLayers(G4int layers)
{
    fECalBarrelLGLayers = layers;
    G4cout << "ECal barrel LG layers set to " << fECalBarrelLGLayers << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalEndcapHGLayers(G4int layers)
{
    fECalEndcapHGLayers = layers;
    G4cout << "ECal end cap HG layers set to " << fECalEndcapHGLayers << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetECalEndcapLGLayers(G4int layers)
{
    fECalEndcapLGLayers = layers;
    G4cout << "ECal end cap LG layers set to " << fECalEndcapLGLayers << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetMuIDAbsorberThickness(G4double thickness)
{
    fMuIDAbsorberThickness = thickness;
    G4cout << "MuID absorber thickness set to " << fMuIDAbsorberThickness/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetMuIDScintillatorThickness(G4double thickness)
{
    fMuIDScintillatorThickness = thickness;
    G4cout << "MuID scintillator thickness set to " << fMuIDScintillatorThickness/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetMuIDLayers(G4int layers)
{
    fMuIDLayers = layers;
    G4cout << "MuID layers set to " << fMuIDLayers << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArNModulesX(G4int nx)
{
    fLArNModulesX = nx;
    G4cout << "Number of LAr modules in the X direction set to " << fLArNModulesX << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArNModulesY(G4int ny)
{
    fLArNModulesY = ny;
    G4cout << "Number of LAr modules in the Y direction set to " << fLArNModulesY << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArNModulesZ(G4int nz)
{
    fLArNModulesZ = nz;
    G4cout << "Number of LAr modules in the Z direction set to " << fLArNModulesZ << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArModuleLength(G4double length)
{
    fLArModuleLength = length;
    G4cout << "LAr module length set to " << fLArModuleLength/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArModuleWidth(G4double width)
{
    fLArModuleWidth = width;
    G4cout << "LAr module width set to " << fLArModuleWidth/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArModuleDepth(G4double depth)
{
    fLArModuleDepth = depth;
    G4cout << "LAr module depth set to " << fLArModuleDepth/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArModuleGap(G4double gap)
{
    fLArModuleGap = gap;
    G4cout << "LAr module gap set to " << fLArModuleGap/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArInsulationThickness(G4double thickness)
{
    fLArInsulationThickness = thickness;
    G4cout << "LAr module insulation thickness set to " << fLArInsulationThickness/mm << " mm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArCryostatThickness(G4double thickness)
{
    fLArCryostatThickness = thickness;
    G4cout << "LAr TPC cryostat thickness set to " << fLArCryostatThickness/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArEnableMuonWindow(G4bool enable)
{
    fLArEnableMuonWindow = enable;
    G4cout << "LAr muon window set to " << fLArEnableMuonWindow << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetLArMuonWindowThickness(G4double thickness)
{
    fLArMuonWindowThickness = thickness;
    G4cout << "LAr muon window thickness set to " << fLArMuonWindowThickness/cm << " cm" << G4endl;

    if (fGeometryInitialized && fGeometryType == kLArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::SetPressure(G4double pressure)
{
    fPressure = pressure;
    G4cout << "Gas pressure set to " << fPressure/bar << " bar" << G4endl;

    if (fGeometryInitialized && fGeometryType == kGArLike) {
        UpdateGeometry();
    }
}

void DetectorConstruction::DefineCommands()
{
    // Create a new messenger
    fMessenger = new DetectorMessenger(this);
}
