//
// SimDataTypes.hh - Geant4 simulation output data classes
//
// These are ROOT-serializable (TObject-based) types written to the "Events"
// TTree by the detector_simulation. They are shared between the simulation
// and reconstruction to avoid duplicating definitions.
//
// Namespace: root  (kept for backwards compatibility with existing files)
//

#ifndef SimDataTypes_hh
#define SimDataTypes_hh

#include <vector>
#include "TObject.h"
#include "TString.h"

namespace root {

    class TrajectoryPoint : public TObject
    {
    public:
        TrajectoryPoint() : x(0), y(0), z(0), time(0) {}
        TrajectoryPoint(Double_t px, Double_t py, Double_t pz, Double_t t)
            : x(px), y(py), z(pz), time(t) {}

        Double_t x, y, z;
        Double_t time;

        ClassDef(TrajectoryPoint, 1)
    };

    class MomentumPoint : public TObject
    {
    public:
        MomentumPoint() : x(0), y(0), z(0) {}
        MomentumPoint(Double_t px, Double_t py, Double_t pz)
            : x(px), y(py), z(pz) {}

        Double_t x, y, z;

        ClassDef(MomentumPoint, 1)
    };

    class Trajectory : public TObject
    {
    public:
        Trajectory() = default;

        std::vector<TrajectoryPoint> points;
        std::vector<MomentumPoint>   mom_points;
        TString start_vol_name;
        TString stop_vol_name;

        ClassDef(Trajectory, 1)
    };

    // Hit in the Time Projection Chamber
    class TPCHit : public TObject
    {
    public:
        TPCHit() : x(0), y(0), z(0), energyDeposit(0), stepSize(0) {}
        TPCHit(Double_t px, Double_t py, Double_t pz, Double_t edep, Double_t step)
            : x(px), y(py), z(pz), energyDeposit(edep), stepSize(step) {}

        Double_t x, y, z;          // Position [cm]
        Double_t energyDeposit;    // Energy deposit [MeV]
        Double_t stepSize;         // Step length [cm]

        ClassDef(TPCHit, 1)
    };

    // Hit in the Electromagnetic Calorimeter
    class ECalHit : public TObject
    {
    public:
        ECalHit() : x(0), y(0), z(0), time(0), energyDeposit(0),
                    segment(0), layer(0), detID(0) {}
        ECalHit(Double_t px, Double_t py, Double_t pz, Double_t t,
                Double_t edep, Int_t s, Int_t l, Int_t id)
            : x(px), y(py), z(pz), time(t), energyDeposit(edep),
              segment(s), layer(l), detID(id) {}

        Double_t x, y, z;          // Position [cm]
        Double_t time;             // Hit time [ns]
        Double_t energyDeposit;    // Energy deposit [MeV]
        Int_t    segment;          // Barrel segment or endcap copy number / 1000
        Int_t    layer;            // Layer index (copy number % 1000)
        Int_t    detID;            // 1 = left endcap, 2 = barrel, 3 = right endcap

        ClassDef(ECalHit, 1)
    };

    // Hit in the Muon ID system
    class MuIDHit : public TObject
    {
    public:
        MuIDHit() : x(0), y(0), z(0), time(0), energyDeposit(0),
                    segment(0), layer(0), detID(0) {}
        MuIDHit(Double_t px, Double_t py, Double_t pz, Double_t t,
                Double_t edep, Int_t s, Int_t l, Int_t id)
            : x(px), y(py), z(pz), time(t), energyDeposit(edep),
              segment(s), layer(l), detID(id) {}

        Double_t x, y, z;
        Double_t time;
        Double_t energyDeposit;
        Int_t    segment;
        Int_t    layer;
        Int_t    detID;

        ClassDef(MuIDHit, 1)
    };

    class Particle : public TObject
    {
    public:
        Particle() : trackID(0), pdgCode(0), creatorProcess(""),
                     endProcess(""), motherID(0) {}
        Particle(Int_t id, Int_t pdg, const TString& process,
                 const TString& end_process, Int_t momid)
            : trackID(id), pdgCode(pdg), creatorProcess(process),
              endProcess(end_process), motherID(momid) {}

        Int_t    trackID;
        Int_t    pdgCode;
        TString  creatorProcess;
        TString  endProcess;
        Int_t    motherID;
        Trajectory           trajectory;
        std::vector<TPCHit>  tpcHits;
        std::vector<ECalHit> ecalHits;
        std::vector<MuIDHit> muidHits;
        std::vector<TPCHit>  sec_tpcHits;
        std::vector<ECalHit> sec_ecalHits;
        std::vector<MuIDHit> sec_muidHits;

        ClassDef(Particle, 1)
    };

    class Event : public TObject
    {
    public:
        Event() : eventID(0) {}
        Event(Int_t id) : eventID(id) {}

        Int_t                  eventID;
        std::vector<Particle>  particles;

        ClassDef(Event, 1)
    };

} // namespace root

#endif
