 /***************************************************************************
 * AnalysisEvent.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Data holders for the three input trees consumed by an analysis:
 *
 *     - GenieEvent   -> GENIE "gst" tree  (truth-level interaction record)
 *     - SimEvent     -> FastGArSim "AnaTree" (flat per-event particle/hit ntuple)
 *     - GeometryInfo -> FastGArSim "GeoTree" (detector configuration, 1 entry)
 *
 *   Each holder knows how to attach itself to a TTree. Branches that are
 *   absent from a given file are reported once and left at their default
 *   values, so the same reader works with older files, gun-only samples, etc.
 *
 ***************************************************************************/

#ifndef AnalysisEvent_hh
#define AnalysisEvent_hh

#include <string>
#include <unordered_map>
#include <vector>

#include "Rtypes.h"
#include "TVector3.h"

class TTree;

namespace ana {

/* -------------------------------------------------------------------------- */
/*                              Useful constants                              */
/* -------------------------------------------------------------------------- */

// Fixed max dimension to read arrays from the gst TTree
constexpr Int_t kNPmax = 1000;

// Sub-detector identifiers stored in ecalHitDetID / muidHitDetID
// (see common/include/SimDataTypes.hh)
enum EDetID {
    kEndCapMinus = 1,
    kBarrel      = 2,
    kEndCapPlus  = 3
};

// Geometry types stored in GeometryInfo::geometry_type
enum EGeometryType {
    kGArLike = 0,
    kLArLike = 1
};

/* -------------------------------------------------------------------------- */
/*                        GENIE gst truth-level record                        */
/* -------------------------------------------------------------------------- */

struct GenieEvent {
    Int_t iev = 0;               // event number
    Int_t neu = 0;               // neutrino PDG code
    Bool_t qel = kFALSE;         // is QE event?
    Bool_t mec = kFALSE;         // is MEC event?
    Bool_t res = kFALSE;         // is RES event?
    Bool_t dis = kFALSE;         // is DIS event?
    Bool_t coh = kFALSE;         // is COH event?
    Int_t resid = 0;             // baryon resonance ID
    Bool_t cc = kFALSE;          // is CC?
    Bool_t nc = kFALSE;          // is NC?
    Double_t Ev = 0.;            // neutrino energy
    Double_t pxv = 0., pyv = 0., pzv = 0.;   // neutrino momentum
    Double_t Q2 = 0.;            // 4-momentum transfer
    Double_t W = 0.;             // hadronic invariant mass
    Double_t x = 0.;             // Bjorken x
    Double_t y = 0.;             // inelasticity
    Double_t El = 0.;            // primary lepton energy
    Double_t pxl = 0., pyl = 0., pzl = 0.;   // primary lepton momentum
    Int_t nf = 0;                // number of final state hadrons
    Int_t pdgf[kNPmax] = {};     // i-th hadron PDG code
    Double_t Ef[kNPmax] = {};    // i-th hadron energy
    Double_t pxf[kNPmax] = {};   // i-th hadron px
    Double_t pyf[kNPmax] = {};   // i-th hadron py
    Double_t pzf[kNPmax] = {};   // i-th hadron pz

    // Attach every branch above to `tree`
    void Connect(TTree* tree);

    // Convenience accessors
    TVector3 NuMomentum() const { return TVector3(pxv, pyv, pzv); }
    TVector3 LeptonMomentum() const { return TVector3(pxl, pyl, pzl); }
    TVector3 HadronMomentum(Int_t i) const { return TVector3(pxf[i], pyf[i], pzf[i]); }
};

/* -------------------------------------------------------------------------- */
/*                     FastGArSim flat ntuple (AnaTree)                       */
/* -------------------------------------------------------------------------- */

struct SimEvent {
    Int_t eventID = 0;

    // Particle properties
    std::vector<Int_t> *trackID = nullptr;
    std::vector<Int_t> *pdgCode = nullptr;
    std::vector<Int_t> *motherID = nullptr;
    std::vector<std::string> *creatorProcess = nullptr;
    std::vector<std::string> *endProcess = nullptr;

    // Start and end trajectory points
    std::vector<Float_t> *startX = nullptr, *startY = nullptr, *startZ = nullptr;
    std::vector<Float_t> *endX = nullptr, *endY = nullptr, *endZ = nullptr;

    // Initial and final momenta
    std::vector<Float_t> *startPX = nullptr, *startPY = nullptr, *startPZ = nullptr;
    std::vector<Float_t> *endPX = nullptr, *endPY = nullptr, *endPZ = nullptr;

    // TPC hits
    std::vector<Int_t> *tpcHitTrackID = nullptr;
    std::vector<Bool_t> *tpcHitIsSec = nullptr;
    std::vector<Float_t> *tpcHitX = nullptr, *tpcHitY = nullptr, *tpcHitZ = nullptr;
    std::vector<Float_t> *tpcHitEdep = nullptr;
    std::vector<Float_t> *tpcHitStepSize = nullptr;

    // ECal hits
    std::vector<Int_t> *ecalHitTrackID = nullptr;
    std::vector<Bool_t> *ecalHitIsSec = nullptr;
    std::vector<Float_t> *ecalHitX = nullptr, *ecalHitY = nullptr, *ecalHitZ = nullptr;
    std::vector<Float_t> *ecalHitTime = nullptr;
    std::vector<Float_t> *ecalHitEdep = nullptr;
    std::vector<Int_t> *ecalHitSegment = nullptr;
    std::vector<Int_t> *ecalHitLayer = nullptr;
    std::vector<Int_t> *ecalHitDetID = nullptr;

    // MuID hits
    std::vector<Int_t> *muidHitTrackID = nullptr;
    std::vector<Bool_t> *muidHitIsSec = nullptr;
    std::vector<Float_t> *muidHitX = nullptr, *muidHitY = nullptr, *muidHitZ = nullptr;
    std::vector<Float_t> *muidHitTime = nullptr;
    std::vector<Float_t> *muidHitEdep = nullptr;
    std::vector<Int_t> *muidHitSegment = nullptr;
    std::vector<Int_t> *muidHitLayer = nullptr;
    std::vector<Int_t> *muidHitDetID = nullptr;

    // Attach every branch above to `tree`
    void Connect(TTree* tree);

    // Invalidate the per-event lookup caches. AnalysisBase calls this after
    // every GetEntry(), so analyses never need to.
    void Update();

    /* ------------------------------ Collection sizes ---------------------- */

    size_t NParticles() const { return trackID ? trackID->size() : 0; }
    size_t NTPCHits()   const { return tpcHitTrackID ? tpcHitTrackID->size() : 0; }
    size_t NECalHits()  const { return ecalHitTrackID ? ecalHitTrackID->size() : 0; }
    size_t NMuIDHits()  const { return muidHitTrackID ? muidHitTrackID->size() : 0; }

    /* ------------------------------ Vector accessors ---------------------- */

    TVector3 StartPosition(size_t i) const;
    TVector3 EndPosition(size_t i) const;
    TVector3 StartMomentum(size_t i) const;
    TVector3 EndMomentum(size_t i) const;

    TVector3 TPCHitPosition(size_t i) const;
    TVector3 ECalHitPosition(size_t i) const;
    TVector3 MuIDHitPosition(size_t i) const;

    /* ------------------------------ Track lookups ------------------------- */

    // Index of the particle with the given Geant4 track ID, or -1 if absent.
    Int_t IndexOfTrack(Int_t id) const;

    // Indices of the hits belonging to a given track ID. The maps are built
    // once per event on first use, so repeated calls are cheap.
    const std::vector<size_t>& TPCHitsOfTrack(Int_t id) const;
    const std::vector<size_t>& ECalHitsOfTrack(Int_t id) const;
    const std::vector<size_t>& MuIDHitsOfTrack(Int_t id) const;

    // Summed energy deposit of a track in each sub-detector
    Double_t TPCEdepOfTrack(Int_t id) const;
    Double_t ECalEdepOfTrack(Int_t id) const;
    Double_t MuIDEdepOfTrack(Int_t id) const;

private:
    using HitIndexMap = std::unordered_map<Int_t, std::vector<size_t>>;

    static const std::vector<size_t>& Lookup(const HitIndexMap& map, Int_t id);
    static void BuildHitMap(const std::vector<Int_t>* ids, HitIndexMap& map);

    mutable std::unordered_map<Int_t, Int_t> fTrackIndex;
    mutable HitIndexMap fTPCHitsByTrack;
    mutable HitIndexMap fECalHitsByTrack;
    mutable HitIndexMap fMuIDHitsByTrack;

    mutable Bool_t fTrackIndexValid = kFALSE;
    mutable Bool_t fTPCMapValid = kFALSE;
    mutable Bool_t fECalMapValid = kFALSE;
    mutable Bool_t fMuIDMapValid = kFALSE;
};

/* -------------------------------------------------------------------------- */
/*                    FastGArSim geometry record (GeoTree)                    */
/* -------------------------------------------------------------------------- */

struct GeometryInfo {
    // Geometry type -- 0 = GArLike, 1 = LArLike
    Int_t geometry_type = 0;

    // GAr TPC parameters
    Double_t gar_tpc_radius = 0.;                 // TPC radius (cm)
    Double_t gar_tpc_length = 0.;                 // TPC length (cm)
    Double_t gar_magnetic_field = 0.;             // magnetic field strength (T)
    Double_t gar_pressure = 0.;                   // gas pressure (bar)

    // ECal parameters
    Double_t ecal_barrel_gap = 0.;                // distance between TPC and ECal barrel (cm)
    Double_t ecal_endcap_gap = 0.;                // distance between TPC and ECal end caps (cm)
    Int_t ecal_num_sides = 0;                     // number of sides of the ECal barrel
    Double_t ecal_hg_absorber_thickness = 0.;     // absorber thickness ECal HG layers
    Double_t ecal_hg_scintillator_thickness = 0.; // scintillator thickness ECal HG layers
    Double_t ecal_hg_board_thickness = 0.;        // PCB thickness ECal HG layers
    Int_t ecal_barrel_hg_layers = 0;              // ECal barrel number of HG layers
    Int_t ecal_endcap_hg_layers = 0;              // ECal end caps number of HG layers
    Double_t ecal_lg_absorber_thickness = 0.;     // absorber thickness ECal LG layers
    Double_t ecal_lg_scintillator_thickness = 0.; // scintillator thickness ECal LG layers
    Int_t ecal_barrel_lg_layers = 0;              // ECal barrel number of LG layers
    Int_t ecal_endcap_lg_layers = 0;              // ECal end caps number of LG layers

    // MuID parameters
    Double_t muid_barrel_gap = 0.;                // distance between ECal and MuID barrels (cm)
    Double_t muid_absorber_thickness = 0.;        // absorber thickness MuID layers
    Double_t muid_scintillator_thickness = 0.;    // scintillator thickness MuID layers
    Int_t muid_num_sides = 0;                     // number of sides of the MuID barrel
    Int_t muid_layers = 0;                        // MuID barrel number of layers

    // LAr TPC parameters
    Int_t lar_n_modules_x = 0;
    Int_t lar_n_modules_y = 0;
    Int_t lar_n_modules_z = 0;
    Double_t lar_module_length = 0.;
    Double_t lar_module_width = 0.;
    Double_t lar_module_depth = 0.;
    Double_t lar_module_gap = 0.;
    Double_t lar_insulation_thickness = 0.;
    Double_t lar_cryostat_thickness = 0.;
    Bool_t lar_enable_muon_window = kFALSE;
    Double_t lar_muon_window_thickness = 0.;

    // Attach every branch above to `tree`
    void Connect(TTree* tree);

    // Dump the configuration to stdout
    void Print() const;

    // Total number of ECal layers in barrel / end caps
    Int_t ECalBarrelLayers() const { return ecal_barrel_hg_layers + ecal_barrel_lg_layers; }
    Int_t ECalEndCapLayers() const { return ecal_endcap_hg_layers + ecal_endcap_lg_layers; }
};

} // namespace ana

#endif
