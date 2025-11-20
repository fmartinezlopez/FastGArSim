#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "DataTypes.hh"
#include "ROOTDataTypes.hh"

#include "G4AnalysisManager.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// Geometry info structure
struct GeometryInfo {
  // Geometry type
  Int_t geometry_type;  // 0 = GArLike, 1 = LArLike
  
  // GAr TPC parameters
  Double_t gar_tpc_radius;
  Double_t gar_tpc_length;
  Double_t gar_magnetic_field;
  Double_t gar_pressure;
  
  // ECal parameters
  Double_t ecal_barrel_gap;
  Double_t ecal_endcap_gap;
  Int_t ecal_num_sides;
  Double_t ecal_hg_absorber_thickness;
  Double_t ecal_hg_scintillator_thickness;
  Double_t ecal_hg_board_thickness;
  Int_t ecal_barrel_hg_layers;
  Int_t ecal_endcap_hg_layers;
  Double_t ecal_lg_absorber_thickness;
  Double_t ecal_lg_scintillator_thickness;
  Int_t ecal_barrel_lg_layers;
  Int_t ecal_endcap_lg_layers;
  
  // MuID parameters
  Double_t muid_barrel_gap;
  Double_t muid_absorber_thickness;
  Double_t muid_scintillator_thickness;
  Int_t muid_num_sides;
  Int_t muid_layers;
  
  // LAr TPC parameters
  Int_t lar_n_modules_x;
  Int_t lar_n_modules_y;
  Int_t lar_n_modules_z;
  Double_t lar_module_length;
  Double_t lar_module_width;
  Double_t lar_module_depth;
  Double_t lar_module_gap;
  Double_t lar_insulation_thickness;
  Double_t lar_cryostat_thickness;
  Bool_t lar_enable_muon_window;
  Double_t lar_muon_window_thickness;
};

class AnalysisManager
{
public:
  static AnalysisManager* GetInstance();
  virtual ~AnalysisManager();
  
  // Initialization methods for the beginning and end of run
  void Book();
  void BookGeometry();
  void Save();
  void Close();

  // Fill Geometry TTree
  void FillGeometryInfo(const GeometryInfo& geoInfo);

  // Method to update eventID when reading external input file(s)
  void UpdateEventID(const G4int id);
  
  // Methods to record data
  void RecordEnergyDeposit(const G4Step* step);
  void RecordTrackInfo(const G4Track* track, G4bool end=false);

  // Methods for Event/Particle-level data
  void BeginEvent(G4int eventID);
  void EndEvent();
  void AddTPCHit(const G4Track* track, const G4ThreeVector& pos, G4double edep, G4double step);
  void AddECalHit(const G4Track* track, const G4ThreeVector& pos, G4double time, G4double edep, G4int segment, G4int layer, G4int detID);
  void AddMuIDHit(const G4Track* track, const G4ThreeVector& pos, G4double time, G4double edep, G4int segment, G4int layer, G4int detID);
  Particle* GetParticle(const G4Track* track);

  // Setter methods
  void SetOutputFileName(const G4String& name) { fOutputFileName = name; }
  void SetEnergyCut(G4double cut);

  // Getter methods
  G4bool GetWriteTrajectory() const { return fWriteTrajectory; }

private:
  AnalysisManager();
  static G4ThreadLocal AnalysisManager* fInstance;
  
  // Helper methods
  std::map<G4int, Particle*> SimplifyParticleCollection();
  void WriteEvent(const Event& g4Event);

  void DefineCommands();
  
  // Data members
  G4String fOutputFileName;
  G4double fEnergyCut;

  // Current event data
  Event* fCurrentEvent;
  G4int fCurrentEventID;
  std::map<G4int, Particle*> fTrackMap;  // Maps trackID to Particle

  // ROOT file and tree
  TFile* fRootFile;
  TTree* fEventTree;
  root::Event* fStoredEvent;  // ROOT Event object that will be stored in the tree
  TTree* fGeometryTree;
  GeometryInfo fGeometryInfo;
  
  // Configuration parameters
  G4bool fWriteTrajectory;
  G4bool fWriteTPCHits;
  G4bool fWriteECalHits;
  G4bool fWriteMuIDHits;

  // Messenger for macro commands
  G4GenericMessenger* fMessenger;
};

#endif