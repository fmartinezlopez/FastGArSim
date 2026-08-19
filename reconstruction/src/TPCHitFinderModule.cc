//
// TPCHitFinderModule.cc - Pulse finding on digitized TPC pad waveforms
//

#include "TPCHitFinderModule.hh"
#include "ModuleFactory.hh"
#include "RecoStore.hh"
#include "TPCConditions.hh"

#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>

// Self-registration
namespace {
    bool kRegistered = ModuleFactory::Instance().Register(
        "TPCHitFinderModule",
        []() -> RecoModule* { return new TPCHitFinderModule(); });
}

// ---------------------------------------------------------------------------
// Constructor / Destructor
// ---------------------------------------------------------------------------

TPCHitFinderModule::TPCHitFinderModule()
    : RecoModule("TPCHitFinder"),
      fWaveforms(nullptr), fHits(nullptr), fConditions(nullptr),
      fThreshold(12.0), fMinPulseTicks(2), fMinPulseADC(0.0),
      fValleyFraction(0.7), fChargeCalibration(1.0), fLifetimeCorrection(true),
      fNEvents(0), fNHits(0), fTotalCharge(0), fTotalEnergy(0)
{
}

TPCHitFinderModule::~TPCHitFinderModule()
{
}

// ---------------------------------------------------------------------------
// Module spec
// ---------------------------------------------------------------------------

std::vector<ObjectSpec> TPCHitFinderModule::GetInputSpec() const
{
    return {
        {"TPCWaveforms",  "std::vector<digi::TPCWaveform>"},
        {"TPCConditions", "TPCConditions"}
    };
}

std::vector<ObjectSpec> TPCHitFinderModule::GetOutputSpec() const
{
    return {{"TPCHits", "std::vector<digi::TPCHit>"}};
}

// ---------------------------------------------------------------------------
// Initialize
// ---------------------------------------------------------------------------

void TPCHitFinderModule::Initialize()
{
    Print("Initializing TPC hit finder module");

    fThreshold          = GetParameterDouble("threshold",         12.0);
    fMinPulseTicks      = GetParameterInt   ("minPulseTicks",     2);
    fMinPulseADC        = GetParameterDouble("minPulseADC",       0.0);
    fValleyFraction     = GetParameterDouble("valleyFraction",    0.7);
    fChargeCalibration  = GetParameterDouble("chargeCalibration", 1.0);
    fLifetimeCorrection = GetParameterBool  ("lifetimeCorrection", true);

    fWaveforms  = fStore->Get<std::vector<digi::TPCWaveform>>("TPCWaveforms");
    fConditions = fStore->Get<TPCConditions>("TPCConditions");

    Print("Hit finding parameters:");
    std::cout << "   Threshold above pedestal [ADC]: " << fThreshold      << "\n"
              << "   Minimum pulse length [ticks]: "   << fMinPulseTicks  << "\n"
              << "   Minimum pulse integral [ADC]: "   << fMinPulseADC    << "\n"
              << "   Peak splitting valley fraction: " << fValleyFraction << "\n"
              << "   Charge calibration: "             << fChargeCalibration << "\n"
              << "   Lifetime correction: "
              << (fLifetimeCorrection ? "on" : "off") << std::endl;

    fHits = new std::vector<digi::TPCHit>();
    fStore->Register("TPCHits", fHits);

    if (fOutputTree) {
        fOutputTree->Branch("TPCHits", &fHits);
    }
}

// ---------------------------------------------------------------------------
// Execute
// ---------------------------------------------------------------------------

void TPCHitFinderModule::Execute()
{
    fHits->clear();
    if (!fWaveforms) return;

    std::vector<float> signal;

    for (const auto& wf : *fWaveforms) {
        // Work on the pedestal-subtracted waveform
        signal.resize(wf.adc.size());
        for (size_t i = 0; i < wf.adc.size(); ++i) {
            signal[i] = wf.adc[i] - wf.pedestal;
        }

        for (const auto& pulse : FindPulses(signal)) {
            fHits->push_back(MakeHit(wf, signal, pulse));
            fTotalCharge += fHits->back().charge;
            fTotalEnergy += fHits->back().energy;
        }
    }

    fNEvents++;
    fNHits += fHits->size();
}

// ---------------------------------------------------------------------------
// Finalize
// ---------------------------------------------------------------------------

void TPCHitFinderModule::Finalize()
{
    Print("TPC hit finding summary:");
    std::cout << "   Events processed: " << fNEvents << "\n"
              << "   Hits found: "       << fNHits   << "\n"
              << "   Reconstructed charge [electrons]: " << fTotalCharge << "\n"
              << "   Reconstructed energy [MeV]: "       << fTotalEnergy << std::endl;
    if (fNEvents > 0) {
        std::cout << "   Hits per event: "
                  << double(fNHits) / fNEvents << std::endl;
    }

    delete fHits;
    fHits = nullptr;
}

// ---------------------------------------------------------------------------
// Pulse finding
// ---------------------------------------------------------------------------

std::vector<TPCHitFinderModule::Pulse>
TPCHitFinderModule::FindPulses(const std::vector<float>& signal) const
{
    std::vector<Pulse> pulses;
    const int n = static_cast<int>(signal.size());

    int i = 0;
    while (i < n) {
        if (signal[i] < fThreshold) { ++i; continue; }

        // One contiguous region above threshold
        const int begin = i;
        while (i < n && signal[i] >= fThreshold) ++i;
        const int end = i;

        // Local maxima inside it are the candidate pulses
        std::vector<int> peaks;
        for (int k = begin; k < end; ++k) {
            const bool risingInto  = (k == begin)   || signal[k] >= signal[k-1];
            const bool fallingAway = (k == end - 1) || signal[k] >  signal[k+1];
            if (risingInto && fallingAway) peaks.push_back(k);
        }
        if (peaks.empty()) peaks.push_back(begin);

        // Merge peaks that are not separated by a deep enough valley, dropping
        // the smaller one each time, until every survivor is well isolated
        bool merged = true;
        while (merged && peaks.size() > 1) {
            merged = false;
            for (size_t k = 0; k + 1 < peaks.size(); ++k) {
                const int valley = static_cast<int>(
                    std::min_element(signal.begin() + peaks[k],
                                     signal.begin() + peaks[k+1] + 1)
                    - signal.begin());
                const float lower = std::min(signal[peaks[k]], signal[peaks[k+1]]);

                if (signal[valley] > fValleyFraction * lower) {
                    const size_t drop =
                        (signal[peaks[k]] < signal[peaks[k+1]]) ? k : k + 1;
                    peaks.erase(peaks.begin() + drop);
                    merged = true;
                    break;
                }
            }
        }

        // Split the region at the valleys between the surviving peaks
        for (size_t k = 0; k < peaks.size(); ++k) {
            Pulse pulse;
            pulse.peak  = peaks[k];
            pulse.begin = (k == 0) ? begin
                        : static_cast<int>(
                              std::min_element(signal.begin() + peaks[k-1],
                                               signal.begin() + peaks[k] + 1)
                              - signal.begin());
            pulse.end   = (k + 1 == peaks.size()) ? end
                        : static_cast<int>(
                              std::min_element(signal.begin() + peaks[k],
                                               signal.begin() + peaks[k+1] + 1)
                              - signal.begin());

            if (pulse.end - pulse.begin < fMinPulseTicks) continue;

            double integral = 0;
            for (int s = pulse.begin; s < pulse.end; ++s) integral += signal[s];
            if (integral < fMinPulseADC) continue;

            pulses.push_back(pulse);
        }
    }

    return pulses;
}

// ---------------------------------------------------------------------------
// Hit building
// ---------------------------------------------------------------------------

digi::TPCHit TPCHitFinderModule::MakeHit(const digi::TPCWaveform& wf,
                                         const std::vector<float>& signal,
                                         const Pulse& pulse) const
{
    const double dt = fConditions->samplePeriod;

    // Charge-weighted moments of the pulse, in absolute ticks
    double sum = 0, sumT = 0, sumT2 = 0;
    for (int i = pulse.begin; i < pulse.end; ++i) {
        const double t = (wf.tickStart + i + 0.5) * dt;
        sum   += signal[i];
        sumT  += signal[i] * t;
        sumT2 += signal[i] * t * t;
    }

    digi::TPCHit hit;
    hit.channel   = wf.channel;
    hit.plane     = wf.plane;
    hit.row       = wf.row;
    hit.col       = wf.col;
    hit.x         = wf.padX;
    hit.y         = wf.padY;
    hit.adcSum    = sum;
    hit.peakAmplitude = signal[pulse.peak];
    hit.startTick = wf.tickStart + pulse.begin;
    hit.endTick   = wf.tickStart + pulse.end - 1;

    // Shaping delays the pulse and broadens it; take both back out
    const double meanTime = (sum > 0) ? sumT / sum : 0.0;
    hit.driftTime = std::max(0.0, meanTime - fConditions->responseCentroid);

    const double variance = (sum > 0) ? sumT2 / sum - meanTime * meanTime : 0.0;
    hit.sigmaT = std::sqrt(std::max(0.0,
        variance - fConditions->responseSigma * fConditions->responseSigma));
    hit.sigmaZ = hit.sigmaT * fConditions->driftVelocity;

    hit.z = fConditions->DriftToZ(wf.plane, hit.driftTime);

    // Integral back to collected electrons, then to the energy that made them
    double charge = (fConditions->adcPerElectron > 0)
                  ? sum / fConditions->adcPerElectron : 0.0;
    charge *= fChargeCalibration;
    if (fLifetimeCorrection && fConditions->electronLifetime > 0) {
        charge *= std::exp(hit.driftTime / fConditions->electronLifetime);
    }
    hit.charge = charge;

    const double efficiency = (fConditions->collectionEff > 0)
                            ? fConditions->collectionEff : 1.0;
    hit.energy = charge * fConditions->wValue * 1.0e-6 / efficiency;

    // MC truth carried along by the waveform
    if (!wf.trueEnergy.empty()) {
        for (int i = pulse.begin; i < pulse.end && i < (int)wf.trueEnergy.size(); ++i) {
            hit.trueEnergy += wf.trueEnergy[i];
        }
        hit.trackIDs       = wf.trackIDs;
        hit.trackFractions = wf.trackFractions;
    }

    return hit;
}
