 /***************************************************************************
 * AnalysisEvent.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Branch wiring and per-event lookup helpers for GenieEvent, SimEvent and
 *   GeometryInfo.
 *
 ***************************************************************************/

#include "AnalysisEvent.hh"

#include <iostream>

#include "TTree.h"

namespace ana {

namespace {

/* -------------------------------------------------------------------------- */
/*                            Branch wiring helpers                           */
/* -------------------------------------------------------------------------- */

// Collects the names of branches a file does not provide, so that a single
// warning can be issued instead of one message per branch.
class BranchConnector {
public:
    BranchConnector(TTree* tree, const char* label) : fTree(tree), fLabel(label) {}

    template <typename T>
    void operator()(const char* name, T* address)
    {
        if (!fTree->GetBranch(name)) {
            fMissing.push_back(name);
            return;
        }
        fTree->SetBranchAddress(name, address);
    }

    ~BranchConnector()
    {
        if (fMissing.empty()) return;
        std::cout << "Warning: tree '" << fTree->GetName() << "' (" << fLabel
                  << ") is missing " << fMissing.size() << " expected branch(es):";
        for (const auto& name : fMissing) std::cout << " " << name;
        std::cout << "\n         these variables keep their default values."
                  << std::endl;
    }

private:
    TTree* fTree;
    const char* fLabel;
    std::vector<std::string> fMissing;
};

} // anonymous namespace

/* -------------------------------------------------------------------------- */
/*                                 GenieEvent                                 */
/* -------------------------------------------------------------------------- */

void GenieEvent::Connect(TTree* tree)
{
    if (!tree) return;

    BranchConnector connect(tree, "GENIE gst");

    connect("iev",   &iev);
    connect("neu",   &neu);
    connect("qel",   &qel);
    connect("mec",   &mec);
    connect("res",   &res);
    connect("dis",   &dis);
    connect("coh",   &coh);
    connect("resid", &resid);
    connect("cc",    &cc);
    connect("nc",    &nc);
    connect("Ev",    &Ev);
    connect("pxv",   &pxv);
    connect("pyv",   &pyv);
    connect("pzv",   &pzv);
    connect("Q2",    &Q2);
    connect("W",     &W);
    connect("x",     &x);
    connect("y",     &y);
    connect("El",    &El);
    connect("pxl",   &pxl);
    connect("pyl",   &pyl);
    connect("pzl",   &pzl);
    connect("nf",    &nf);
    connect("pdgf",  pdgf);
    connect("Ef",    Ef);
    connect("pxf",   pxf);
    connect("pyf",   pyf);
    connect("pzf",   pzf);
}

/* -------------------------------------------------------------------------- */
/*                                  SimEvent                                  */
/* -------------------------------------------------------------------------- */

namespace {

// Gather one sub-detector's hits from every particle into a flat list. The
// order is the one an analysis expects when it walks the particles itself:
// particle by particle, each particle's own hits before its secondaries'.
template <class Hit>
void FillView(const root::Event* event,
              const std::vector<Hit> root::Particle::* primary,
              const std::vector<Hit> root::Particle::* secondary,
              HitView<Hit>& view)
{
    view.hits.clear();
    view.trackIDs.clear();
    view.secondary.clear();
    view.byTrack.clear();
    view.built = kTRUE;

    if (!event) return;

    // One pass to size the vectors, so the pointers gathered below are not
    // moved about while the second pass is taking them
    size_t total = 0;
    for (const root::Particle& particle : event->particles) {
        total += (particle.*primary).size() + (particle.*secondary).size();
    }
    view.hits.reserve(total);
    view.trackIDs.reserve(total);
    view.secondary.reserve(total);

    for (const root::Particle& particle : event->particles) {
        for (const Hit& hit : particle.*primary) {
            view.byTrack[particle.trackID].push_back(view.hits.size());
            view.hits.push_back(&hit);
            view.trackIDs.push_back(particle.trackID);
            view.secondary.push_back(kFALSE);
        }
        for (const Hit& hit : particle.*secondary) {
            view.byTrack[particle.trackID].push_back(view.hits.size());
            view.hits.push_back(&hit);
            view.trackIDs.push_back(particle.trackID);
            view.secondary.push_back(kTRUE);
        }
    }
}

// Element access that returns a harmless default rather than reading past the
// end, matching how the flat-ntuple reader behaved
template <class T>
T At(const std::vector<T>& values, size_t i, T fallback = T())
{
    return (i < values.size()) ? values[i] : fallback;
}

TVector3 PointToVector(const root::TrajectoryPoint& point)
{
    return TVector3(point.x, point.y, point.z);
}

TVector3 PointToVector(const root::MomentumPoint& point)
{
    return TVector3(point.x, point.y, point.z);
}

template <class Hit>
Double_t SumEdep(const HitView<Hit>& view, const std::vector<size_t>& indices)
{
    Double_t sum = 0.;
    for (const size_t i : indices) {
        if (i < view.hits.size()) sum += view.hits[i]->energyDeposit;
    }
    return sum;
}

} // anonymous namespace

void SimEvent::Connect(TTree* tree, const char* branchName)
{
    if (!tree) return;

    if (!tree->GetBranch(branchName)) {
        std::cout << "Warning: tree '" << tree->GetName() << "' has no '"
                  << branchName << "' branch; the simulation record will be empty."
                  << std::endl;
        return;
    }

    tree->SetBranchAddress(branchName, &event);
}

void SimEvent::Update()
{
    eventID = event ? event->eventID : 0;

    fTPC.built = kFALSE;
    fECal.built = kFALSE;
    fMuID.built = kFALSE;
    fTrackIndexValid = kFALSE;
}

/* ------------------------------- Particles ------------------------------- */

Int_t SimEvent::TrackID(size_t i) const
{
    return (i < NParticles()) ? event->particles[i].trackID : -1;
}

Int_t SimEvent::PdgCode(size_t i) const
{
    return (i < NParticles()) ? event->particles[i].pdgCode : 0;
}

Int_t SimEvent::MotherID(size_t i) const
{
    return (i < NParticles()) ? event->particles[i].motherID : -1;
}

std::string SimEvent::CreatorProcess(size_t i) const
{
    return (i < NParticles()) ? event->particles[i].creatorProcess.Data() : std::string();
}

std::string SimEvent::EndProcess(size_t i) const
{
    return (i < NParticles()) ? event->particles[i].endProcess.Data() : std::string();
}

Bool_t SimEvent::HasTrajectory(size_t i) const
{
    return i < NParticles() && !event->particles[i].trajectory.points.empty();
}

TVector3 SimEvent::StartPosition(size_t i) const
{
    if (!HasTrajectory(i)) return TVector3();
    return PointToVector(event->particles[i].trajectory.points.front());
}

TVector3 SimEvent::EndPosition(size_t i) const
{
    if (!HasTrajectory(i)) return TVector3();
    return PointToVector(event->particles[i].trajectory.points.back());
}

TVector3 SimEvent::StartMomentum(size_t i) const
{
    if (i >= NParticles()) return TVector3();
    const auto& points = event->particles[i].trajectory.mom_points;
    if (points.empty()) return TVector3();
    return PointToVector(points.front());
}

TVector3 SimEvent::EndMomentum(size_t i) const
{
    if (i >= NParticles()) return TVector3();
    const auto& points = event->particles[i].trajectory.mom_points;
    if (points.empty()) return TVector3();
    return PointToVector(points.back());
}

/* ---------------------------------- Hits --------------------------------- */

const HitView<root::TPCHit>& SimEvent::TPCView() const
{
    if (!fTPC.built) {
        FillView(event, &root::Particle::tpcHits, &root::Particle::sec_tpcHits, fTPC);
    }
    return fTPC;
}

const HitView<root::ECalHit>& SimEvent::ECalView() const
{
    if (!fECal.built) {
        FillView(event, &root::Particle::ecalHits, &root::Particle::sec_ecalHits, fECal);
    }
    return fECal;
}

const HitView<root::MuIDHit>& SimEvent::MuIDView() const
{
    if (!fMuID.built) {
        FillView(event, &root::Particle::muidHits, &root::Particle::sec_muidHits, fMuID);
    }
    return fMuID;
}

size_t SimEvent::NTPCHits()  const { return TPCView().hits.size(); }
size_t SimEvent::NECalHits() const { return ECalView().hits.size(); }
size_t SimEvent::NMuIDHits() const { return MuIDView().hits.size(); }

const root::TPCHit&  SimEvent::TPCHit(size_t k)  const { return *TPCView().hits[k]; }
const root::ECalHit& SimEvent::ECalHit(size_t k) const { return *ECalView().hits[k]; }
const root::MuIDHit& SimEvent::MuIDHit(size_t k) const { return *MuIDView().hits[k]; }

TVector3 SimEvent::TPCHitPosition(size_t k) const
{
    const HitView<root::TPCHit>& view = TPCView();
    if (k >= view.hits.size()) return TVector3();
    return TVector3(view.hits[k]->x, view.hits[k]->y, view.hits[k]->z);
}

TVector3 SimEvent::ECalHitPosition(size_t k) const
{
    const HitView<root::ECalHit>& view = ECalView();
    if (k >= view.hits.size()) return TVector3();
    return TVector3(view.hits[k]->x, view.hits[k]->y, view.hits[k]->z);
}

TVector3 SimEvent::MuIDHitPosition(size_t k) const
{
    const HitView<root::MuIDHit>& view = MuIDView();
    if (k >= view.hits.size()) return TVector3();
    return TVector3(view.hits[k]->x, view.hits[k]->y, view.hits[k]->z);
}

Int_t SimEvent::TPCHitTrackID(size_t k)  const { return At(TPCView().trackIDs,  k, Int_t(-1)); }
Int_t SimEvent::ECalHitTrackID(size_t k) const { return At(ECalView().trackIDs, k, Int_t(-1)); }
Int_t SimEvent::MuIDHitTrackID(size_t k) const { return At(MuIDView().trackIDs, k, Int_t(-1)); }

Bool_t SimEvent::TPCHitIsSecondary(size_t k)  const { return At(TPCView().secondary,  k, kFALSE); }
Bool_t SimEvent::ECalHitIsSecondary(size_t k) const { return At(ECalView().secondary, k, kFALSE); }
Bool_t SimEvent::MuIDHitIsSecondary(size_t k) const { return At(MuIDView().secondary, k, kFALSE); }

/* ------------------------------ Track lookups ---------------------------- */

Int_t SimEvent::IndexOfTrack(Int_t id) const
{
    if (!fTrackIndexValid) {
        fTrackIndex.clear();
        for (size_t i = 0; i < NParticles(); ++i) {
            fTrackIndex.emplace(event->particles[i].trackID, static_cast<Int_t>(i));
        }
        fTrackIndexValid = kTRUE;
    }

    const auto it = fTrackIndex.find(id);
    return (it == fTrackIndex.end()) ? -1 : it->second;
}

const std::vector<size_t>& SimEvent::Lookup(
    const std::unordered_map<Int_t, std::vector<size_t>>& map, Int_t id)
{
    static const std::vector<size_t> kEmpty;
    const auto it = map.find(id);
    return (it == map.end()) ? kEmpty : it->second;
}

const std::vector<size_t>& SimEvent::TPCHitsOfTrack(Int_t id) const
{ return Lookup(TPCView().byTrack, id); }

const std::vector<size_t>& SimEvent::ECalHitsOfTrack(Int_t id) const
{ return Lookup(ECalView().byTrack, id); }

const std::vector<size_t>& SimEvent::MuIDHitsOfTrack(Int_t id) const
{ return Lookup(MuIDView().byTrack, id); }

Double_t SimEvent::TPCEdepOfTrack(Int_t id) const
{ return SumEdep(TPCView(), TPCHitsOfTrack(id)); }

Double_t SimEvent::ECalEdepOfTrack(Int_t id) const
{ return SumEdep(ECalView(), ECalHitsOfTrack(id)); }

Double_t SimEvent::MuIDEdepOfTrack(Int_t id) const
{ return SumEdep(MuIDView(), MuIDHitsOfTrack(id)); }

/* -------------------------------------------------------------------------- */
/*                                GeometryInfo                                */
/* -------------------------------------------------------------------------- */

void GeometryInfo::Connect(TTree* tree)
{
    if (!tree) return;

    BranchConnector connect(tree, "FastGArSim Geometry");

    connect("geometry_type", &geometry_type);

    // GAr TPC
    connect("gar_tpc_radius",     &gar_tpc_radius);
    connect("gar_tpc_length",     &gar_tpc_length);
    connect("gar_magnetic_field", &gar_magnetic_field);
    connect("gar_pressure",       &gar_pressure);

    // ECal
    connect("ecal_barrel_gap",                &ecal_barrel_gap);
    connect("ecal_endcap_gap",                &ecal_endcap_gap);
    connect("ecal_num_sides",                 &ecal_num_sides);
    connect("ecal_hg_absorber_thickness",     &ecal_hg_absorber_thickness);
    connect("ecal_hg_scintillator_thickness", &ecal_hg_scintillator_thickness);
    connect("ecal_hg_board_thickness",        &ecal_hg_board_thickness);
    connect("ecal_barrel_hg_layers",          &ecal_barrel_hg_layers);
    connect("ecal_endcap_hg_layers",          &ecal_endcap_hg_layers);
    connect("ecal_lg_absorber_thickness",     &ecal_lg_absorber_thickness);
    connect("ecal_lg_scintillator_thickness", &ecal_lg_scintillator_thickness);
    connect("ecal_barrel_lg_layers",          &ecal_barrel_lg_layers);
    connect("ecal_endcap_lg_layers",          &ecal_endcap_lg_layers);

    // MuID
    connect("muid_barrel_gap",             &muid_barrel_gap);
    connect("muid_absorber_thickness",     &muid_absorber_thickness);
    connect("muid_scintillator_thickness", &muid_scintillator_thickness);
    connect("muid_num_sides",              &muid_num_sides);
    connect("muid_layers",                 &muid_layers);

    // LAr TPC
    connect("lar_n_modules_x",           &lar_n_modules_x);
    connect("lar_n_modules_y",           &lar_n_modules_y);
    connect("lar_n_modules_z",           &lar_n_modules_z);
    connect("lar_module_length",         &lar_module_length);
    connect("lar_module_width",          &lar_module_width);
    connect("lar_module_depth",          &lar_module_depth);
    connect("lar_module_gap",            &lar_module_gap);
    connect("lar_insulation_thickness",  &lar_insulation_thickness);
    connect("lar_cryostat_thickness",    &lar_cryostat_thickness);
    connect("lar_enable_muon_window",    &lar_enable_muon_window);
    connect("lar_muon_window_thickness", &lar_muon_window_thickness);
}

void GeometryInfo::Print() const
{
    std::cout << "Detector geometry\n"
              << "    type:           " << (geometry_type == kLArLike ? "LArLike" : "GArLike") << "\n";

    if (geometry_type == kLArLike) {
        std::cout << "    LAr modules:    " << lar_n_modules_x << " x "
                  << lar_n_modules_y << " x " << lar_n_modules_z << "\n"
                  << "    module size:    " << lar_module_length << " x "
                  << lar_module_width << " x " << lar_module_depth << " cm\n";
    } else {
        std::cout << "    TPC radius:     " << gar_tpc_radius << " cm\n"
                  << "    TPC length:     " << gar_tpc_length << " cm\n"
                  << "    magnetic field: " << gar_magnetic_field << " T\n"
                  << "    gas pressure:   " << gar_pressure << " bar\n";
    }

    std::cout << "    ECal sides:     " << ecal_num_sides << "\n"
              << "    ECal layers:    " << ECalBarrelLayers() << " (barrel), "
              << ECalEndCapLayers() << " (end caps)\n"
              << "    MuID layers:    " << muid_layers << "\n"
              << std::endl;
}

} // namespace ana
