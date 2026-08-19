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

// Safe element access for the optional vector branches
template <typename T>
T Element(const std::vector<T>* v, size_t i, T fallback = T())
{
    return (v && i < v->size()) ? v->at(i) : fallback;
}

TVector3 MakeVector(const std::vector<Float_t>* vx,
                    const std::vector<Float_t>* vy,
                    const std::vector<Float_t>* vz,
                    size_t i)
{
    return TVector3(Element<Float_t>(vx, i), Element<Float_t>(vy, i), Element<Float_t>(vz, i));
}

Double_t SumOver(const std::vector<size_t>& indices, const std::vector<Float_t>* values)
{
    if (!values) return 0.;
    Double_t sum = 0.;
    for (size_t i : indices) {
        if (i < values->size()) sum += values->at(i);
    }
    return sum;
}

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

void SimEvent::Connect(TTree* tree)
{
    if (!tree) return;

    BranchConnector connect(tree, "FastGArSim AnaTree");

    connect("eventID", &eventID);

    // Particle properties
    connect("trackID",        &trackID);
    connect("pdgCode",        &pdgCode);
    connect("motherID",       &motherID);
    connect("creatorProcess", &creatorProcess);
    connect("endProcess",     &endProcess);

    // Start and end trajectory points
    connect("startX", &startX);
    connect("startY", &startY);
    connect("startZ", &startZ);
    connect("endX",   &endX);
    connect("endY",   &endY);
    connect("endZ",   &endZ);

    // Initial and final momenta
    connect("startPX", &startPX);
    connect("startPY", &startPY);
    connect("startPZ", &startPZ);
    connect("endPX",   &endPX);
    connect("endPY",   &endPY);
    connect("endPZ",   &endPZ);

    // TPC hits
    connect("tpcHitTrackID",  &tpcHitTrackID);
    connect("tpcHitIsSec",    &tpcHitIsSec);
    connect("tpcHitX",        &tpcHitX);
    connect("tpcHitY",        &tpcHitY);
    connect("tpcHitZ",        &tpcHitZ);
    connect("tpcHitEdep",     &tpcHitEdep);
    connect("tpcHitStepSize", &tpcHitStepSize);

    // ECal hits
    connect("ecalHitTrackID", &ecalHitTrackID);
    connect("ecalHitIsSec",   &ecalHitIsSec);
    connect("ecalHitX",       &ecalHitX);
    connect("ecalHitY",       &ecalHitY);
    connect("ecalHitZ",       &ecalHitZ);
    connect("ecalHitTime",    &ecalHitTime);
    connect("ecalHitEdep",    &ecalHitEdep);
    connect("ecalHitSegment", &ecalHitSegment);
    connect("ecalHitLayer",   &ecalHitLayer);
    connect("ecalHitDetID",   &ecalHitDetID);

    // MuID hits
    connect("muidHitTrackID", &muidHitTrackID);
    connect("muidHitIsSec",   &muidHitIsSec);
    connect("muidHitX",       &muidHitX);
    connect("muidHitY",       &muidHitY);
    connect("muidHitZ",       &muidHitZ);
    connect("muidHitTime",    &muidHitTime);
    connect("muidHitEdep",    &muidHitEdep);
    connect("muidHitSegment", &muidHitSegment);
    connect("muidHitLayer",   &muidHitLayer);
    connect("muidHitDetID",   &muidHitDetID);
}

void SimEvent::Update()
{
    fTrackIndexValid = kFALSE;
    fTPCMapValid = kFALSE;
    fECalMapValid = kFALSE;
    fMuIDMapValid = kFALSE;
}

TVector3 SimEvent::StartPosition(size_t i) const { return MakeVector(startX, startY, startZ, i); }
TVector3 SimEvent::EndPosition(size_t i) const { return MakeVector(endX, endY, endZ, i); }
TVector3 SimEvent::StartMomentum(size_t i) const { return MakeVector(startPX, startPY, startPZ, i); }
TVector3 SimEvent::EndMomentum(size_t i) const { return MakeVector(endPX, endPY, endPZ, i); }

TVector3 SimEvent::TPCHitPosition(size_t i) const { return MakeVector(tpcHitX, tpcHitY, tpcHitZ, i); }
TVector3 SimEvent::ECalHitPosition(size_t i) const { return MakeVector(ecalHitX, ecalHitY, ecalHitZ, i); }
TVector3 SimEvent::MuIDHitPosition(size_t i) const { return MakeVector(muidHitX, muidHitY, muidHitZ, i); }

Int_t SimEvent::IndexOfTrack(Int_t id) const
{
    if (!fTrackIndexValid) {
        fTrackIndex.clear();
        if (trackID) {
            for (size_t i = 0; i < trackID->size(); ++i) {
                fTrackIndex.emplace(trackID->at(i), static_cast<Int_t>(i));
            }
        }
        fTrackIndexValid = kTRUE;
    }

    auto it = fTrackIndex.find(id);
    return (it == fTrackIndex.end()) ? -1 : it->second;
}

void SimEvent::BuildHitMap(const std::vector<Int_t>* ids, HitIndexMap& map)
{
    map.clear();
    if (!ids) return;
    for (size_t i = 0; i < ids->size(); ++i) {
        map[ids->at(i)].push_back(i);
    }
}

const std::vector<size_t>& SimEvent::Lookup(const HitIndexMap& map, Int_t id)
{
    static const std::vector<size_t> kEmpty;
    auto it = map.find(id);
    return (it == map.end()) ? kEmpty : it->second;
}

const std::vector<size_t>& SimEvent::TPCHitsOfTrack(Int_t id) const
{
    if (!fTPCMapValid) {
        BuildHitMap(tpcHitTrackID, fTPCHitsByTrack);
        fTPCMapValid = kTRUE;
    }
    return Lookup(fTPCHitsByTrack, id);
}

const std::vector<size_t>& SimEvent::ECalHitsOfTrack(Int_t id) const
{
    if (!fECalMapValid) {
        BuildHitMap(ecalHitTrackID, fECalHitsByTrack);
        fECalMapValid = kTRUE;
    }
    return Lookup(fECalHitsByTrack, id);
}

const std::vector<size_t>& SimEvent::MuIDHitsOfTrack(Int_t id) const
{
    if (!fMuIDMapValid) {
        BuildHitMap(muidHitTrackID, fMuIDHitsByTrack);
        fMuIDMapValid = kTRUE;
    }
    return Lookup(fMuIDHitsByTrack, id);
}

Double_t SimEvent::TPCEdepOfTrack(Int_t id) const { return SumOver(TPCHitsOfTrack(id), tpcHitEdep); }
Double_t SimEvent::ECalEdepOfTrack(Int_t id) const { return SumOver(ECalHitsOfTrack(id), ecalHitEdep); }
Double_t SimEvent::MuIDEdepOfTrack(Int_t id) const { return SumOver(MuIDHitsOfTrack(id), muidHitEdep); }

/* -------------------------------------------------------------------------- */
/*                                GeometryInfo                                */
/* -------------------------------------------------------------------------- */

void GeometryInfo::Connect(TTree* tree)
{
    if (!tree) return;

    BranchConnector connect(tree, "FastGArSim GeoTree");

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
