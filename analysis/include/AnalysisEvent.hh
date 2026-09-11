 /***************************************************************************
 * AnalysisEvent.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Data holders for the three input trees consumed by an analysis:
 *
 *     - GenieEvent   -> GENIE "gst" tree (truth-level interaction record)
 *     - SimEvent     -> FastGArSim "Events" tree (root::Event objects)
 *     - GeometryInfo -> FastGArSim "Geometry" tree (detector configuration)
 *
 *   Each holder knows how to attach itself to a TTree. GenieEvent and
 *   GeometryInfo read one branch per field, and a branch a given file does
 *   not have is reported once and left at its default value, so the same
 *   reader works with older files, gun-only samples, and files with no
 *   geometry record.
 *
 *   SimEvent reads the simulation's objects as they were written: no flat
 *   ntuple stage in between. The event itself is reachable as sim.event, and
 *   on top of it SimEvent offers the two views an analysis usually wants --
 *   particles by index, and every hit of a sub-detector in one flat list,
 *   regardless of which particle produced it. The flat list is built lazily,
 *   once per event, and only if something asks for it.
 *
 ***************************************************************************/

#ifndef AnalysisEvent_hh
#define AnalysisEvent_hh

#include <string>
#include <unordered_map>
#include <vector>

#include "Rtypes.h"
#include "TVector3.h"

#include "SimDataTypes.hh"

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
/*                    FastGArSim simulation output (Events)                   */
/* -------------------------------------------------------------------------- */

// A flat view of one sub-detector's hits, gathered from every particle in the
// event. Built on first use and thrown away when the next event is loaded.
template <class Hit>
struct HitView {
    std::vector<const Hit*> hits;
    std::vector<Int_t> trackIDs;      // particle each hit belongs to
    std::vector<Bool_t> secondary;    // was it in the particle's sec_ list?
    std::unordered_map<Int_t, std::vector<size_t>> byTrack;
    Bool_t built = kFALSE;
};

struct SimEvent {

    // The event exactly as the simulation wrote it: ROOT refills it on every
    // GetEntry(), so treat it as read-only. Null until Connect(). This is the
    // thing to reach for when the views below do not cover what you want.
    root::Event* event = nullptr;

    Int_t eventID = 0;   // refreshed every entry

    // Attach to the "Event" branch of `tree`
    void Connect(TTree* tree, const char* branchName = "Event");

    // Drop the per-event views. AnalysisBase calls this after every
    // GetEntry(), so analyses never need to.
    void Update();

    Bool_t IsValid() const { return event != nullptr; }

    /* ------------------------------- Particles ---------------------------- */

    size_t NParticles() const { return event ? event->particles.size() : 0; }

    // Undefined for i >= NParticles(); check first, as with any vector
    const root::Particle& Particle(size_t i) const { return event->particles[i]; }

    Int_t TrackID(size_t i) const;
    Int_t PdgCode(size_t i) const;
    Int_t MotherID(size_t i) const;
    std::string CreatorProcess(size_t i) const;
    std::string EndProcess(size_t i) const;

    // First and last point of the stored trajectory. A particle with no
    // trajectory points gives the null vector.
    TVector3 StartPosition(size_t i) const;
    TVector3 EndPosition(size_t i) const;
    TVector3 StartMomentum(size_t i) const;
    TVector3 EndMomentum(size_t i) const;

    // Whether there is a trajectory to take those from
    Bool_t HasTrajectory(size_t i) const;

    /* --------------------------------- Hits ------------------------------- */
    /* One flat list per sub-detector, over every particle in the event, in    */
    /* particle order with each particle's own hits before its secondaries'.   */

    size_t NTPCHits()  const;
    size_t NECalHits() const;
    size_t NMuIDHits() const;

    const root::TPCHit&  TPCHit(size_t k) const;
    const root::ECalHit& ECalHit(size_t k) const;
    const root::MuIDHit& MuIDHit(size_t k) const;

    TVector3 TPCHitPosition(size_t k) const;
    TVector3 ECalHitPosition(size_t k) const;
    TVector3 MuIDHitPosition(size_t k) const;

    // The particle a hit belongs to, and whether it was booked as one of that
    // particle's secondaries (a delta ray folded back onto its parent)
    Int_t TPCHitTrackID(size_t k) const;
    Int_t ECalHitTrackID(size_t k) const;
    Int_t MuIDHitTrackID(size_t k) const;

    Bool_t TPCHitIsSecondary(size_t k) const;
    Bool_t ECalHitIsSecondary(size_t k) const;
    Bool_t MuIDHitIsSecondary(size_t k) const;

    /* ------------------------------ Track lookups ------------------------- */

    // Index of the particle with the given Geant4 track ID, or -1 if absent.
    Int_t IndexOfTrack(Int_t id) const;

    // Indices into the flat hit lists above. The maps are built once per
    // event on first use, so repeated calls are cheap.
    const std::vector<size_t>& TPCHitsOfTrack(Int_t id) const;
    const std::vector<size_t>& ECalHitsOfTrack(Int_t id) const;
    const std::vector<size_t>& MuIDHitsOfTrack(Int_t id) const;

    // Summed energy deposit of a track in each sub-detector
    Double_t TPCEdepOfTrack(Int_t id) const;
    Double_t ECalEdepOfTrack(Int_t id) const;
    Double_t MuIDEdepOfTrack(Int_t id) const;

private:
    const HitView<root::TPCHit>&  TPCView() const;
    const HitView<root::ECalHit>& ECalView() const;
    const HitView<root::MuIDHit>& MuIDView() const;

    static const std::vector<size_t>& Lookup(
        const std::unordered_map<Int_t, std::vector<size_t>>& map, Int_t id);

    mutable HitView<root::TPCHit>  fTPC;
    mutable HitView<root::ECalHit> fECal;
    mutable HitView<root::MuIDHit> fMuID;

    mutable std::unordered_map<Int_t, Int_t> fTrackIndex;
    mutable Bool_t fTrackIndexValid = kFALSE;
};

/* -------------------------------------------------------------------------- */
/*                   FastGArSim geometry record (Geometry)                    */
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
