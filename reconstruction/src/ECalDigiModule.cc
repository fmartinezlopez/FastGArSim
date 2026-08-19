//
// ECalDigiModule.cc - ECal digitization reconstruction module
//

#include "ECalDigiModule.hh"
#include "ModuleFactory.hh"
#include "RecoStore.hh"
#include "SimDataTypes.hh"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <cmath>
#include <algorithm>
#include <iostream>

// Self-registration
namespace {
    bool kRegistered = ModuleFactory::Instance().Register(
        "ECalDigiModule",
        []() -> RecoModule* { return new ECalDigiModule(); });
}

// ---------------------------------------------------------------------------
// Constructor / Destructor
// ---------------------------------------------------------------------------

ECalDigiModule::ECalDigiModule()
    : RecoModule("ECalDigi"),
      fSimEvent(nullptr),
      fTileHits(nullptr),
      fStripHits(nullptr),
      // Geometry (filled from tree)
      fGarTpcRadius(0), fGarTpcLength(0),
      fEcalBarrelGap(0), fEcalEndcapGap(0),
      fEcalNumSides(0),
      fEcalHgAbsorberThickness(0), fEcalHgScintillatorThickness(0), fEcalHgBoardThickness(0),
      fEcalBarrelHgLayers(0), fEcalEndcapHgLayers(0),
      fEcalLgAbsorberThickness(0), fEcalLgScintillatorThickness(0),
      fEcalBarrelLgLayers(0), fEcalEndcapLgLayers(0),
      // Derived geometry (computed in Initialize)
      fInnerApothem(0), fHgLayerThickness(0), fLgLayerThickness(0),
      fBarrelLength(0), fEndcapRadius(0),
      fBaseHgScintillatorRadius(0), fBaseLgScintillatorRadius(0),
      // Digitization defaults (overridden by macro parameters)
      fTileSize(0.5f), fStripSize(0.5f),
      fTimeRes(0.1f), fPropVelocity(10.0f),
      fLightYield(1000.0f), fEffectivePE(5000.0f),
      fSiPMNoise(1.0f), fSiPMGain(10.0f),
      fAttLength(150.0f),
      fADCRange(4096), fADCThreshold(10),
      fMinClusterBins(2), fMaxPropTime(50.0f)
{
}

ECalDigiModule::~ECalDigiModule()
{
}

// ---------------------------------------------------------------------------
// Module spec
// ---------------------------------------------------------------------------

std::vector<ObjectSpec> ECalDigiModule::GetInputSpec() const
{
    // Reads directly from the "Event" branch of the simulation Events TTree
    return {{"Event", "root::Event"}};
}

std::vector<ObjectSpec> ECalDigiModule::GetOutputSpec() const
{
    return {
        {"ECalTileDigiHits",  "std::vector<digi::TileDigiHit>"},
        {"ECalStripDigiHits", "std::vector<digi::StripDigiHit>"}
    };
}

// ---------------------------------------------------------------------------
// Initialize
// ---------------------------------------------------------------------------

void ECalDigiModule::Initialize()
{
    Print("Initializing ECal digitization module");

    // ---- Read digitization parameters from macro ----
    fTileSize        = GetParameterDouble("tileSize",       0.5);
    fStripSize       = GetParameterDouble("stripSize",      0.5);
    fTimeRes         = GetParameterDouble("timeRes",        0.1);
    fPropVelocity    = GetParameterDouble("propVelocity",   10.0);
    fLightYield      = GetParameterDouble("lightYield",     1000.0);
    fEffectivePE     = GetParameterDouble("effectivePE",    5000.0);
    fSiPMNoise       = GetParameterDouble("sipmNoise",      1.0);
    fSiPMGain        = GetParameterDouble("sipmGain",       10.0);
    fAttLength       = GetParameterDouble("attLength",      150.0);
    fADCRange        = GetParameterInt   ("adcRange",       4096);
    fADCThreshold    = GetParameterInt   ("adcThreshold",   10);
    fMinClusterBins  = GetParameterInt   ("minClusterBins", 2);
    fMaxPropTime     = GetParameterDouble("maxPropTime",    50.0);

    Print("Digitization parameters:");
    std::cout << "   Tile size [cm]: "      << fTileSize      << "\n"
              << "   Strip size [cm]: "     << fStripSize     << "\n"
              << "   Time resolution [ns]: "<< fTimeRes       << "\n"
              << "   Light yield [pe/MeV]: "<< fLightYield    << "\n"
              << "   SiPM effective PE: "   << fEffectivePE   << "\n"
              << "   ADC range: "           << fADCRange      << "\n"
              << "   ADC threshold: "       << fADCThreshold  << "\n"
              << "   Attenuation length [cm]: " << fAttLength << std::endl;

    // ---- Load detector geometry from the Geometry TTree ----
    if (!fInputFile) {
        std::cerr << "[ECalDigi] ERROR: no input file available" << std::endl;
        return;
    }

    TTree* geoTree = (TTree*)fInputFile->Get("Geometry");
    if (!geoTree) {
        std::cerr << "[ECalDigi] ERROR: cannot find 'Geometry' tree in input file" << std::endl;
        return;
    }

    geoTree->SetBranchAddress("gar_tpc_radius",                &fGarTpcRadius);
    geoTree->SetBranchAddress("gar_tpc_length",                &fGarTpcLength);
    geoTree->SetBranchAddress("ecal_barrel_gap",               &fEcalBarrelGap);
    geoTree->SetBranchAddress("ecal_endcap_gap",               &fEcalEndcapGap);
    geoTree->SetBranchAddress("ecal_num_sides",                &fEcalNumSides);
    geoTree->SetBranchAddress("ecal_hg_absorber_thickness",    &fEcalHgAbsorberThickness);
    geoTree->SetBranchAddress("ecal_hg_scintillator_thickness",&fEcalHgScintillatorThickness);
    geoTree->SetBranchAddress("ecal_hg_board_thickness",       &fEcalHgBoardThickness);
    geoTree->SetBranchAddress("ecal_barrel_hg_layers",         &fEcalBarrelHgLayers);
    geoTree->SetBranchAddress("ecal_endcap_hg_layers",         &fEcalEndcapHgLayers);
    geoTree->SetBranchAddress("ecal_lg_absorber_thickness",    &fEcalLgAbsorberThickness);
    geoTree->SetBranchAddress("ecal_lg_scintillator_thickness",&fEcalLgScintillatorThickness);
    geoTree->SetBranchAddress("ecal_barrel_lg_layers",         &fEcalBarrelLgLayers);
    geoTree->SetBranchAddress("ecal_endcap_lg_layers",         &fEcalEndcapLgLayers);
    geoTree->GetEntry(0);

    // ---- Compute derived geometry ----
    fInnerApothem  = fGarTpcRadius + fEcalBarrelGap;
    fHgLayerThickness = fEcalHgAbsorberThickness
                      + fEcalHgScintillatorThickness
                      + fEcalHgBoardThickness;
    fLgLayerThickness = fEcalLgAbsorberThickness
                      + fEcalLgScintillatorThickness;
    fBarrelLength  = fGarTpcLength / 2.0
                   + fEcalEndcapGap
                   + fHgLayerThickness * fEcalEndcapHgLayers
                   + fLgLayerThickness * fEcalEndcapLgLayers;
    fEndcapRadius  = fGarTpcRadius + fEcalBarrelGap;

    double innerRadius = fInnerApothem / std::cos(TMath::Pi() / 12.0);
    fBaseHgScintillatorRadius = innerRadius
        + fEcalHgAbsorberThickness + fEcalHgScintillatorThickness * 0.5;
    fBaseLgScintillatorRadius = fBaseHgScintillatorRadius
        + fEcalLgAbsorberThickness + fEcalLgScintillatorThickness * 0.5;

    fRotAngles.clear();
    for (int i = 0; i < fEcalNumSides; ++i) {
        fRotAngles.push_back(TMath::Pi() - (i+1) * (2.0*TMath::Pi() / fEcalNumSides));
    }

    Print("Geometry loaded:");
    std::cout << "   TPC radius [cm]: " << fGarTpcRadius    << "\n"
              << "   TPC length [cm]: " << fGarTpcLength     << "\n"
              << "   ECal sides: "      << fEcalNumSides     << "\n"
              << "   HG layers (barrel/endcap): "
              << fEcalBarrelHgLayers << "/" << fEcalEndcapHgLayers << "\n"
              << "   LG layers (barrel/endcap): "
              << fEcalBarrelLgLayers << "/" << fEcalEndcapLgLayers << std::endl;

    // ---- Set up input branch (Events TTree) ----
    fSimEvent = nullptr;
    fInputTree->SetBranchAddress("Event", &fSimEvent);

    // ---- Create output containers, register in store, create TTree branches ----
    fTileHits  = new std::vector<digi::TileDigiHit>();
    fStripHits = new std::vector<digi::StripDigiHit>();

    fStore->Register("ECalTileDigiHits",  fTileHits);
    fStore->Register("ECalStripDigiHits", fStripHits);

    if (fOutputTree) {
        fOutputTree->Branch("ECalTileDigiHits",  &fTileHits);
        fOutputTree->Branch("ECalStripDigiHits", &fStripHits);
    }

    // Seed RNG
    std::random_device rd;
    fRng.seed(rd());
}

// ---------------------------------------------------------------------------
// Execute
// ---------------------------------------------------------------------------

void ECalDigiModule::Execute()
{
    fTileHits->clear();
    fStripHits->clear();

    if (!fSimEvent) return;

    HitMap   tileMap;
    StripMap stripMap;
    tileMap.reserve(10000);

    MapHitsToReadout(fSimEvent, tileMap, stripMap);
    DigitizeTiles(tileMap);
    DigitizeStrips(stripMap);
}

// ---------------------------------------------------------------------------
// Finalize
// ---------------------------------------------------------------------------

void ECalDigiModule::Finalize()
{
    Print("ECal digitization summary:");
    std::cout << "   Total events processed" << std::endl;

    delete fTileHits;
    delete fStripHits;
    fTileHits  = nullptr;
    fStripHits = nullptr;
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

LayerType ECalDigiModule::GetLayerType(int detID, int layer) const
{
    if (detID == 2) {  // barrel
        return (layer < fEcalBarrelHgLayers) ? kHG : kLG;
    } else {           // endcaps (detID 1 or 3)
        return (layer < fEcalEndcapHgLayers) ? kHG : kLG;
    }
}

float ECalDigiModule::ComputeStripLength(int detID, int layer, int row,
                                         float localRadius) const
{
    bool isEndcap   = (detID == 1 || detID == 3);
    bool stripAlongX = (layer % 2 == 0);

    if (isEndcap) {
        float xCircle = row * fStripSize;
        float rSq = fEndcapRadius * fEndcapRadius - xCircle * xCircle;
        return (rSq > 0) ? 2.0f * std::sqrt(rSq) : 0.0f;
    } else {
        if (stripAlongX) {
            // Strip spans from -fBarrelLength to +fBarrelLength (fBarrelLength is the half-length)
            return 2.0f * fBarrelLength;
        } else {
            return 2.0f * localRadius
                 * std::tan(TMath::Pi() / fEcalNumSides);
        }
    }
}

// ---------------------------------------------------------------------------
// Map raw ECal hits to readout channels
// ---------------------------------------------------------------------------

void ECalDigiModule::MapHitsToReadout(const root::Event* event,
                                      HitMap& tileMap, StripMap& stripMap) const
{
    for (const auto& particle : event->particles) {
        // Process both primary and secondary ECal hits
        auto processHits = [&](const std::vector<root::ECalHit>& hits) {
            for (const auto& hit : hits) {
                int   detID   = hit.detID;
                int   segment = hit.segment;
                int   layer   = hit.layer;
                float x       = hit.x;
                float y       = hit.y;
                float z       = hit.z;
                float t       = hit.time;
                float energy  = hit.energyDeposit;
                int   trackID = particle.trackID;

                LayerType lType = GetLayerType(detID, layer);
                bool  isEndcap  = (detID == 1 || detID == 3);

                // Compute local (row/col) coordinates
                float xlocal, ylocal, localRadius = 0.0f;

                if (isEndcap) {
                    xlocal = x;
                    ylocal = y;
                } else {
                    // Rotate into barrel segment frame
                    float theta = fRotAngles[segment];
                    float c = std::cos(theta - TMath::PiOver2());
                    float s = std::sin(theta - TMath::PiOver2());
                    float yr = s * x + c * y;
                    float xr = c * x - s * y;

                    xlocal = z;
                    if (lType == kHG) {
                        localRadius = fBaseHgScintillatorRadius
                                    + layer * fHgLayerThickness;
                    } else {
                        localRadius = fBaseLgScintillatorRadius
                                    + (layer - fEcalBarrelHgLayers) * fLgLayerThickness;
                    }
                    ylocal = (xr != 0) ? localRadius * (yr / xr) : 0.0f;
                }

                if (lType == kHG) {
                    int row = static_cast<int>(xlocal / fTileSize);
                    int col = static_cast<int>(ylocal / fTileSize);
                    tileMap[{segment, layer, row, col}].push_back({energy, t, trackID});

                } else {
                    bool stripAlongX = (layer % 2 == 0);
                    int  row         = static_cast<int>(
                        (stripAlongX ? ylocal : xlocal) / fStripSize);
                    float posAlong   = stripAlongX ? xlocal : ylocal;
                    float stripLen   = ComputeStripLength(detID, layer, row, localRadius);

                    // Distance from hit to each strip end
                    float dLeft  = stripLen/2.0f - posAlong;
                    float dRight = stripLen/2.0f + posAlong;

                    float t1 = t + dLeft  / fPropVelocity;
                    float t2 = t + dRight / fPropVelocity;

                    // Clamp distances to [0, stripLen] to guard against floating-point
                    // edge cases at strip boundaries
                    dLeft  = std::max(0.0f, std::min(dLeft,  stripLen));
                    dRight = std::max(0.0f, std::min(dRight, stripLen));

                    // Exponential light attenuation; fall back to 0.5 if
                    // fAttLength <= 0 (backwards-compatible / no-attenuation mode)
                    float attLeft, attRight;
                    if (fAttLength > 0) {
                        attLeft  = std::exp(-dLeft  / fAttLength);
                        attRight = std::exp(-dRight / fAttLength);
                    } else {
                        attLeft = attRight = 0.5f;
                    }

                    StripInfo& info = stripMap[{segment, layer, row}];
                    info.stripLength = stripLen;
                    info.isEndcap    = isEndcap;
                    info.hitsLeft .push_back({energy * attLeft,  t1, trackID});
                    info.hitsRight.push_back({energy * attRight, t2, trackID});
                    info.hitsTrue .push_back({energy,            t,  trackID});
                }
            }
        };

        processHits(particle.ecalHits);
        processHits(particle.sec_ecalHits);
    }
}

// ---------------------------------------------------------------------------
// Waveform digitization
// ---------------------------------------------------------------------------

ECalDigiModule::Waveform
ECalDigiModule::DigitizeHits(const std::vector<RawHit>& hits)
{
    Waveform w;
    if (hits.empty()) return w;

    float tMin = hits[0].time, tMax = hits[0].time;
    for (const auto& h : hits) {
        if (h.time < tMin) tMin = h.time;
        if (h.time > tMax) tMax = h.time;
    }
    tMin -= fTimeRes;
    tMax += fTimeRes;

    w.nBins = std::max(1, static_cast<int>(std::ceil((tMax - tMin) / fTimeRes)));
    w.tMin  = tMin;
    w.bins.resize(w.nBins);
    for (int i = 0; i < w.nBins; ++i) {
        w.bins[i].time = tMin + (i + 0.5f) * fTimeRes;
    }
    for (const auto& h : hits) {
        int idx = static_cast<int>((h.time - tMin) / fTimeRes);
        if (idx < 0) idx = 0;
        if (idx >= w.nBins) idx = w.nBins - 1;
        w.bins[idx].add(h);
    }

    // Photon statistics: energy → mean PE → SiPM saturation → binomial smearing
    w.adc.resize(w.nBins);
    std::normal_distribution<float> noise(0.0f, fSiPMNoise);
    for (int i = 0; i < w.nBins; ++i) {
        float meanPE = w.bins[i].energy * fLightYield;
        float satPE  = fEffectivePE * (1.0f - std::exp(-meanPE / fEffectivePE));
        float p = std::min(std::max(satPE / fEffectivePE, 0.0f), 1.0f);

        std::binomial_distribution<int> binom(static_cast<int>(fEffectivePE), p);
        float signal = binom(fRng) + noise(fRng);
        float raw    = signal * fSiPMGain;
        raw = std::min(std::max(raw, 0.0f), float(fADCRange - 1));
        w.adc[i] = static_cast<int>(std::floor(raw));
    }
    return w;
}

std::vector<ECalDigiModule::Cluster>
ECalDigiModule::FindClusters(const Waveform& w) const
{
    std::vector<Cluster> clusters;
    if (w.nBins == 0) return clusters;

    int i = 0;
    while (i < w.nBins) {
        if (w.adc[i] < fADCThreshold) { ++i; continue; }

        Cluster c;
        int startBin = i;
        while (i < w.nBins && w.adc[i] >= fADCThreshold) {
            c.adcSum     += w.adc[i];
            c.time       += w.bins[i].time * w.adc[i];
            c.trueEnergy += w.bins[i].energy;
            for (const auto& kv : w.bins[i].trackContrib) {
                c.trackEnergy[kv.first] += kv.second;
            }
            ++i;
        }
        int endBin = i;
        if (c.adcSum > 0) c.time /= c.adcSum;
        if (endBin - startBin >= fMinClusterBins) {
            clusters.push_back(std::move(c));
        }
    }
    return clusters;
}

// ---------------------------------------------------------------------------
// Output hit builders
// ---------------------------------------------------------------------------

digi::TileDigiHit ECalDigiModule::MakeTileHit(const HitKey& key,
                                               const Cluster& c) const
{
    digi::TileDigiHit dh;
    dh.segment    = key.segment;
    dh.layer      = key.layer;
    dh.row        = key.row;
    dh.col        = key.col;
    dh.adcSum     = c.adcSum;
    dh.time       = c.time;
    dh.trueEnergy = c.trueEnergy;
    for (const auto& kv : c.trackEnergy) {
        dh.trackIDs.push_back(kv.first);
        dh.trackFractions.push_back(c.trueEnergy > 0 ? kv.second / c.trueEnergy : 0);
    }
    return dh;
}

digi::StripDigiHit ECalDigiModule::MakeStripHit(
    const StripKey& key, float stripLen,
    const Cluster* cL, const Cluster* cR,
    float totalE,
    const std::unordered_map<int,float>& trackE) const
{
    digi::StripDigiHit sdh;
    sdh.segment     = key.segment;
    sdh.layer       = key.layer;
    sdh.row         = key.row;
    sdh.stripLength = stripLen;
    sdh.adcLeft     = cL ? cL->adcSum : 0;
    sdh.adcRight    = cR ? cR->adcSum : 0;
    sdh.adcCombined = std::sqrt(sdh.adcLeft * sdh.adcRight);
    sdh.timeLeft    = cL ? cL->time   : 0;
    sdh.timeRight   = cR ? cR->time   : 0;
    // Reconstructed position and time from propagation time difference
    sdh.recoPosition = (sdh.timeRight - sdh.timeLeft) * fPropVelocity / 2.0f;
    sdh.recoTime     = (sdh.timeLeft + sdh.timeRight) / 2.0f
                     - stripLen / (2.0f * fPropVelocity);
    sdh.trueEnergy  = totalE;
    for (const auto& kv : trackE) {
        sdh.trackIDs.push_back(kv.first);
        sdh.trackFractions.push_back(totalE > 0 ? kv.second / totalE : 0);
    }
    return sdh;
}

// ---------------------------------------------------------------------------
// Per-event digitization
// ---------------------------------------------------------------------------

void ECalDigiModule::DigitizeTiles(const HitMap& tileMap)
{
    for (const auto& kv : tileMap) {
        std::vector<Cluster> clusters = FindClusters(DigitizeHits(kv.second));
        for (const auto& c : clusters) {
            fTileHits->push_back(MakeTileHit(kv.first, c));
        }
    }
}

void ECalDigiModule::DigitizeStrips(const StripMap& stripMap)
{
    for (const auto& kv : stripMap) {
        const StripKey&  key  = kv.first;
        const StripInfo& info = kv.second;

        // Collect true (unattenuated) energy and track contributions
        float totalE = 0;
        std::unordered_map<int,float> trackE;
        for (const auto& h : info.hitsTrue) {
            totalE            += h.energy;
            trackE[h.trackID] += h.energy;
        }

        auto clustL = FindClusters(DigitizeHits(info.hitsLeft));
        auto clustR = FindClusters(DigitizeHits(info.hitsRight));

        std::vector<bool> rMatched(clustR.size(), false);

        for (const auto& cL : clustL) {
            bool matched = false;
            for (size_t iR = 0; iR < clustR.size(); ++iR) {
                const auto& cR = clustR[iR];
                // Time coincidence check (allow for propagation time)
                if (std::abs(cL.time - cR.time) <= info.stripLength/fPropVelocity + fMaxPropTime) {
                    fStripHits->push_back(MakeStripHit(key, info.stripLength,
                                                       &cL, &cR, totalE, trackE));
                    rMatched[iR] = true;
                    matched      = true;
                }
            }
            if (!matched)
                fStripHits->push_back(MakeStripHit(key, info.stripLength,
                                                   &cL, nullptr, totalE, trackE));
        }
        for (size_t iR = 0; iR < clustR.size(); ++iR) {
            if (!rMatched[iR])
                fStripHits->push_back(MakeStripHit(key, info.stripLength,
                                                   nullptr, &clustR[iR], totalE, trackE));
        }
    }
}
