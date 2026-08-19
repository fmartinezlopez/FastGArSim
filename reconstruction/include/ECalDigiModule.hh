//
// ECalDigiModule.hh - ECal digitization reconstruction module
//
// Reads raw root::ECalHit objects directly from the detector_simulation
// "Events" TTree (no ntuple conversion needed) and applies a realistic SiPM
// response model:
//
//   HG layers (tiles)  - single SiPM per tile
//   LG layers (strips) - dual SiPM with propagation-time position reconstruction
//
// Outputs (registered in the RecoStore and written as TTree branches):
//   "ECalTileDigiHits"  -> std::vector<digi::TileDigiHit>
//   "ECalStripDigiHits" -> std::vector<digi::StripDigiHit>
//

#ifndef ECalDigiModule_h
#define ECalDigiModule_h 1

#include "RecoModule.hh"
#include "DigiDataTypes.hh"

#include <vector>
#include <unordered_map>
#include <random>
#include <tuple>

namespace root { class Event; }

// ---- Internal helper types (not persisted) --------------------------------

enum LayerType { kHG, kLG };

struct RawHit {
    float   energy;
    float   time;
    int     trackID;
};

struct HitKey {
    int segment, layer, row, col;
    bool operator==(const HitKey& o) const {
        return segment==o.segment && layer==o.layer && row==o.row && col==o.col;
    }
};

struct HitKeyHash {
    std::size_t operator()(const HitKey& k) const {
        std::size_t h = 0;
        auto mix = [&](int v) {
            h ^= std::hash<int>{}(v) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
        };
        mix(k.segment); mix(k.layer); mix(k.row); mix(k.col);
        return h;
    }
};

struct StripKey {
    int segment, layer, row;
    bool operator==(const StripKey& o) const {
        return segment==o.segment && layer==o.layer && row==o.row;
    }
};

struct StripKeyHash {
    std::size_t operator()(const StripKey& k) const {
        std::size_t h = 0;
        auto mix = [&](int v) {
            h ^= std::hash<int>{}(v) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
        };
        mix(k.segment); mix(k.layer); mix(k.row);
        return h;
    }
};

struct StripInfo {
    std::vector<RawHit> hitsLeft;
    std::vector<RawHit> hitsRight;
    std::vector<RawHit> hitsTrue;   // unattenuated hits for MC truth accounting
    float stripLength = 0;
    bool  isEndcap    = false;
};

using HitMap   = std::unordered_map<HitKey,  std::vector<RawHit>, HitKeyHash>;
using StripMap = std::unordered_map<StripKey, StripInfo,           StripKeyHash>;

// ---- Module ---------------------------------------------------------------

class ECalDigiModule : public RecoModule {
public:
    ECalDigiModule();
    virtual ~ECalDigiModule();

    virtual void Initialize() override;
    virtual void Execute()    override;
    virtual void Finalize()   override;

    virtual std::vector<ObjectSpec> GetInputSpec()  const override;
    virtual std::vector<ObjectSpec> GetOutputSpec() const override;

private:
    // --- Digitization pipeline ---
    LayerType GetLayerType(int detID, int layer) const;
    float     ComputeStripLength(int detID, int layer, int row, float localRadius) const;
    void      MapHitsToReadout(const root::Event* event,
                               HitMap& tileMap, StripMap& stripMap) const;

    // Waveform digitization
    struct TimeBin {
        float energy = 0, time = 0;
        std::unordered_map<int,float> trackContrib;
        void add(const RawHit& h) { energy += h.energy; trackContrib[h.trackID] += h.energy; }
    };
    struct Waveform {
        std::vector<int>     adc;
        std::vector<TimeBin> bins;
        float tMin = 0;
        int   nBins = 0;
    };
    Waveform DigitizeHits(const std::vector<RawHit>& hits);

    struct Cluster {
        float adcSum = 0, time = 0, trueEnergy = 0;
        std::unordered_map<int,float> trackEnergy;
    };
    std::vector<Cluster> FindClusters(const Waveform& w) const;

    digi::TileDigiHit  MakeTileHit (const HitKey&  key, const Cluster& c) const;
    digi::StripDigiHit MakeStripHit(const StripKey& key, float stripLen,
                                    const Cluster* cL, const Cluster* cR,
                                    float totalE,
                                    const std::unordered_map<int,float>& trackE) const;

    void DigitizeTiles (const HitMap&   tileMap);
    void DigitizeStrips(const StripMap& stripMap);

    // --- I/O ---
    root::Event*                    fSimEvent;      // branch address on input tree
    std::vector<digi::TileDigiHit>*  fTileHits;
    std::vector<digi::StripDigiHit>* fStripHits;

    // --- Geometry (loaded from Geometry TTree) ---
    double fGarTpcRadius, fGarTpcLength;
    double fEcalBarrelGap, fEcalEndcapGap;
    int    fEcalNumSides;
    double fEcalHgAbsorberThickness, fEcalHgScintillatorThickness, fEcalHgBoardThickness;
    int    fEcalBarrelHgLayers, fEcalEndcapHgLayers;
    double fEcalLgAbsorberThickness, fEcalLgScintillatorThickness;
    int    fEcalBarrelLgLayers, fEcalEndcapLgLayers;

    // Derived geometry
    float fInnerApothem;
    float fHgLayerThickness, fLgLayerThickness;
    float fBarrelLength, fEndcapRadius;
    float fBaseHgScintillatorRadius, fBaseLgScintillatorRadius;
    std::vector<float> fRotAngles;   // rotation angles for each barrel segment

    // --- Digitization parameters ---
    float fTileSize, fStripSize;
    float fTimeRes, fPropVelocity;
    float fLightYield, fEffectivePE, fSiPMNoise, fSiPMGain;
    float fAttLength;   // Scintillator attenuation length [cm]; 0 = no attenuation
    int   fADCRange, fADCThreshold, fMinClusterBins;
    float fMaxPropTime;

    std::mt19937 fRng;
};

#endif
