//
// TPCDigiModule.cc - TPC readout simulation
//

#include "TPCDigiModule.hh"
#include "ModuleFactory.hh"
#include "RecoStore.hh"
#include "SimDataTypes.hh"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>

// Self-registration
namespace {
    bool kRegistered = ModuleFactory::Instance().Register(
        "TPCDigiModule",
        []() -> RecoModule* { return new TPCDigiModule(); });

    // Shaping function is truncated once it has decayed away
    constexpr double kResponseLengthInTau = 10.0;
}

// ---------------------------------------------------------------------------
// Constructor / Destructor
// ---------------------------------------------------------------------------

TPCDigiModule::TPCDigiModule()
    : RecoModule("TPCDigi"),
      fSimEvent(nullptr),
      fWaveforms(nullptr),
      fReadout(nullptr),
      // Readout defaults (overridden by macro parameters)
      fPadPitch(0.6),
      fResponseWidth(2.5), fResponseShape(4.0), fResponseRange(1),
      fNormalizeResponse(true),
      fElectronsPerCluster(20.0), fMinClusters(10),
      fNoiseLevel(5.0), fADCThreshold(12.0), fROIPadding(5),
      fMaxTicks(0), fSaveTruth(true), fWriteWaveforms(false),
      fNEvents(0), fNDeposits(0), fNWaveforms(0), fTotalElectrons(0)
{
}

TPCDigiModule::~TPCDigiModule()
{
    delete fReadout;
}

// ---------------------------------------------------------------------------
// Module spec
// ---------------------------------------------------------------------------

std::vector<ObjectSpec> TPCDigiModule::GetInputSpec() const
{
    return {{"Event", "root::Event"}};
}

std::vector<ObjectSpec> TPCDigiModule::GetOutputSpec() const
{
    return {
        {"TPCWaveforms",       "std::vector<digi::TPCWaveform>"},
        {"TPCReadoutGeometry", "TPCReadoutGeometry"},
        {"TPCConditions",      "TPCConditions"}
    };
}

// ---------------------------------------------------------------------------
// Initialize
// ---------------------------------------------------------------------------

void TPCDigiModule::Initialize()
{
    Print("Initializing TPC digitization module");

    // ---- Readout parameters ----
    fPadPitch            = GetParameterDouble("padPitch",           0.6);
    fResponseWidth       = GetParameterDouble("padResponseWidth",   2.5);
    fResponseShape       = GetParameterDouble("padResponseShape",   4.0);
    fResponseRange       = GetParameterInt   ("padResponseRange",   1);
    fNormalizeResponse   = GetParameterBool  ("padResponseNormalize", true);
    fElectronsPerCluster = GetParameterDouble("electronsPerCluster", 20.0);
    fMinClusters         = GetParameterInt   ("minClusters",        10);

    // ---- Drift and gas parameters ----
    fConditions.driftVelocity    = GetParameterDouble("driftVelocity",    3.011);
    fConditions.electronLifetime = GetParameterDouble("electronLifetime", 3.0e6);
    fConditions.diffusionT       = GetParameterDouble("diffusionT",       0.0160);
    fConditions.diffusionL       = GetParameterDouble("diffusionL",       0.0201);
    fConditions.wValue           = GetParameterDouble("wValue",           26.4);
    fConditions.fanoFactor       = GetParameterDouble("fanoFactor",       0.16);
    fConditions.collectionEff    = GetParameterDouble("collectionEfficiency", 0.7);

    // ---- Electronics parameters ----
    const double samplingRate = GetParameterDouble("samplingRate", 20.0);  // [MHz]
    fConditions.samplePeriod  = (samplingRate > 0) ? 1.0 / samplingRate : 0.05;
    fConditions.shapingTime   = GetParameterDouble("shapingTime", 0.1);
    fConditions.gain          = GetParameterDouble("gain",        0.5);
    fConditions.pedestal      = GetParameterDouble("adcPedestal", 100.0);
    fConditions.adcRange      = GetParameterDouble("adcRange",    4096.0);
    fNoiseLevel               = GetParameterDouble("noiseLevel",  5.0);

    // ---- Zero suppression and output ----
    fADCThreshold   = GetParameterDouble("adcThreshold",   12.0);
    fROIPadding     = GetParameterInt   ("roiPadding",     5);
    fMaxTicks       = GetParameterInt   ("maxTicks",       0);
    fSaveTruth      = GetParameterBool  ("saveTruth",      true);
    fWriteWaveforms = GetParameterBool  ("writeWaveforms", false);

    const bool doubleSided = GetParameterBool("doubleSided", true);
    fConditions.nPlanes    = doubleSided ? 2 : 1;
    fConditions.readoutSide = GetParameterInt("readoutSide", 1);

    // ---- Detector geometry from the Geometry TTree ----
    double tpcRadius = 250.0, tpcLength = 500.0;
    double pressure = 0.0, magneticField = 0.0;

    if (!fInputFile) {
        std::cerr << "[TPCDigi] ERROR: no input file available" << std::endl;
        return;
    }

    TTree* geoTree = (TTree*)fInputFile->Get("Geometry");
    if (!geoTree) {
        std::cerr << "[TPCDigi] ERROR: cannot find 'Geometry' tree in input file"
                  << std::endl;
        return;
    }

    geoTree->SetBranchAddress("gar_tpc_radius",     &tpcRadius);
    geoTree->SetBranchAddress("gar_tpc_length",     &tpcLength);
    geoTree->SetBranchAddress("gar_pressure",       &pressure);
    geoTree->SetBranchAddress("gar_magnetic_field", &magneticField);
    geoTree->GetEntry(0);

    fConditions.halfLength = tpcLength / 2.0;

    fReadout = new TPCReadoutGeometry(tpcRadius, fPadPitch, fConditions.nPlanes,
                                      fResponseWidth, fResponseShape, fResponseRange);

    const size_t nNeighbours = (2 * fResponseRange + 1) * (2 * fResponseRange + 1);
    fPadWeights.resize(nNeighbours);
    fPadRows.resize(nNeighbours);
    fPadCols.resize(nNeighbours);

    // ---- Shaping function and the charge calibration that follows from it ----
    BuildResponseFunction();

    if (fMaxTicks <= 0) {
        const double maxTime = fConditions.MaxDrift() / fConditions.driftVelocity;
        fMaxTicks = static_cast<int>(std::ceil(maxTime / fConditions.samplePeriod))
                  + static_cast<int>(fResponse.size());
    }

    Print("Geometry loaded:");
    std::cout << "   TPC radius [cm]: "   << tpcRadius     << "\n"
              << "   TPC length [cm]: "   << tpcLength     << "\n"
              << "   Gas pressure [bar]: "<< pressure      << "\n"
              << "   Magnetic field [T]: "<< magneticField << "\n"
              << "   Readout: " << (doubleSided ? "double-sided (cathode at z = 0)"
                                                : "single-sided") << "\n"
              << "   Maximum drift [cm]: " << fConditions.MaxDrift() << std::endl;
    fReadout->Print();

    Print("Response parameters:");
    std::cout << "   W value [eV]: "            << fConditions.wValue          << "\n"
              << "   Fano factor: "             << fConditions.fanoFactor      << "\n"
              << "   Collection efficiency: "   << fConditions.collectionEff   << "\n"
              << "   Drift velocity [cm/us]: "  << fConditions.driftVelocity   << "\n"
              << "   Electron lifetime [us]: "  << fConditions.electronLifetime<< "\n"
              << "   Diffusion T/L [cm/sqrt(cm)]: " << fConditions.diffusionT
              << " / " << fConditions.diffusionL << "\n"
              << "   Sample period [us]: "      << fConditions.samplePeriod    << "\n"
              << "   Shaping time [us]: "       << fConditions.shapingTime     << "\n"
              << "   Gain [ADC/electron]: "     << fConditions.gain            << "\n"
              << "   ADC per electron (integral): " << fConditions.adcPerElectron << "\n"
              << "   Shaping centroid/sigma [us]: " << fConditions.responseCentroid
              << " / " << fConditions.responseSigma << "\n"
              << "   Noise [electrons]: "       << fNoiseLevel                 << "\n"
              << "   Zero suppression [ADC]: "  << fADCThreshold               << "\n"
              << "   Maximum ticks: "           << fMaxTicks                   << std::endl;

    // ---- Input branch ----
    fSimEvent = nullptr;
    fInputTree->SetBranchAddress("Event", &fSimEvent);

    // ---- Outputs ----
    fWaveforms = new std::vector<digi::TPCWaveform>();

    fStore->Register("TPCWaveforms",       fWaveforms);
    fStore->Register("TPCReadoutGeometry", fReadout);
    fStore->Register("TPCConditions",      &fConditions);

    if (fOutputTree && fWriteWaveforms) {
        fOutputTree->Branch("TPCWaveforms", &fWaveforms);
    } else {
        Print("Waveforms are kept in memory only "
              "(set writeWaveforms true to write them out)");
    }

    const int seed = GetParameterInt("seed", 0);
    if (seed > 0) {
        fRng.seed(seed);
    } else {
        std::random_device rd;
        fRng.seed(rd());
    }
}

// ---------------------------------------------------------------------------
// Execute
// ---------------------------------------------------------------------------

void TPCDigiModule::Execute()
{
    if (!fWaveforms) return;   // Initialize bailed out
    fWaveforms->clear();

    if (!fSimEvent || !fReadout) return;

    PadMap pads;
    pads.reserve(4096);

    DriftDeposits(fSimEvent, pads);
    DigitizePads(pads);

    fNEvents++;
    fNWaveforms += fWaveforms->size();
}

// ---------------------------------------------------------------------------
// Finalize
// ---------------------------------------------------------------------------

void TPCDigiModule::Finalize()
{
    Print("TPC digitization summary:");
    std::cout << "   Events digitized: "        << fNEvents    << "\n"
              << "   Energy deposits drifted: " << fNDeposits  << "\n"
              << "   Electrons collected: "     << fTotalElectrons << "\n"
              << "   Waveforms produced: "      << fNWaveforms << std::endl;
    if (fNEvents > 0) {
        std::cout << "   Waveforms per event: "
                  << double(fNWaveforms) / fNEvents << std::endl;
    }

    delete fWaveforms;
    fWaveforms = nullptr;
}

// ---------------------------------------------------------------------------
// Electronics shaping
// ---------------------------------------------------------------------------

void TPCDigiModule::BuildResponseFunction()
{
    // CR-RC^2 shaping, normalized to unit peak height so that `gain` is the
    // ADC amplitude one collected electron produces at the top of the pulse
    const double tau = fConditions.shapingTime;
    const double dt  = fConditions.samplePeriod;
    const int    n   = std::max(2, static_cast<int>(
                           std::ceil(kResponseLengthInTau * tau / dt)));

    fResponse.assign(n, 0.0);
    const double peak = 4.0 * std::exp(-2.0);   // max of (t/tau)^2 exp(-t/tau)

    double integral = 0.0, firstMoment = 0.0, secondMoment = 0.0;
    for (int k = 0; k < n; ++k) {
        const double u = k * dt / tau;
        fResponse[k] = u * u * std::exp(-u) / peak;

        integral     += fResponse[k];
        firstMoment  += fResponse[k] * k * dt;
        secondMoment += fResponse[k] * k * dt * k * dt;
    }

    // A pulse of one electron spreads over the shaping function, so its ADC
    // integral -- what the hit finder measures -- is gain times this sum
    fConditions.adcPerElectron = fConditions.gain * integral;

    // Shaping delays the pulse and widens it; the hit finder takes both out
    fConditions.responseCentroid = firstMoment / integral;
    fConditions.responseSigma    = std::sqrt(
        std::max(0.0, secondMoment / integral
                    - fConditions.responseCentroid * fConditions.responseCentroid));
}

// ---------------------------------------------------------------------------
// Ionization
// ---------------------------------------------------------------------------

double TPCDigiModule::IonizationElectrons(double energyDeposit)
{
    // Mean pairs from the W value, fluctuating with the Fano factor. Below a
    // handful of pairs the Gaussian is a poor description, so use Poisson.
    const double mean = energyDeposit * 1.0e6 / fConditions.wValue;
    if (mean <= 0) return 0.0;

    if (mean < 20.0) {
        std::poisson_distribution<int> poisson(mean);
        return poisson(fRng);
    }

    std::normal_distribution<double> gauss(
        mean, std::sqrt(fConditions.fanoFactor * mean));
    return std::max(0.0, gauss(fRng));
}

// ---------------------------------------------------------------------------
// Drift
// ---------------------------------------------------------------------------

void TPCDigiModule::DriftDeposits(const root::Event* event, PadMap& pads)
{
    std::normal_distribution<double> unitGauss(0.0, 1.0);

    for (const auto& particle : event->particles) {

        auto driftHits = [&](const std::vector<root::TPCHit>& hits) {
            for (const auto& hit : hits) {

                const int    plane = fConditions.PlaneOf(hit.z);
                const double drift = std::max(0.0,
                                              fConditions.DriftDistance(plane, hit.z));

                double electrons = IonizationElectrons(hit.energyDeposit);
                if (electrons <= 0) continue;

                fNDeposits++;

                // Attenuation along the drift and collection at the pad plane
                const double driftTime = drift / fConditions.driftVelocity;
                electrons *= std::exp(-driftTime / fConditions.electronLifetime);
                electrons *= fConditions.collectionEff;
                if (electrons <= 0) continue;

                // Diffusion grows as the square root of the drift distance
                const double sigmaT = fConditions.diffusionT * std::sqrt(drift);
                const double sigmaL = fConditions.diffusionL * std::sqrt(drift);
                const double sigmaTime = sigmaL / fConditions.driftVelocity;

                // Drift macro-electrons rather than single electrons
                int nClusters = static_cast<int>(
                    std::ceil(electrons / fElectronsPerCluster));
                nClusters = std::max(nClusters, fMinClusters);
                nClusters = std::min(nClusters,
                                     std::max(1, static_cast<int>(electrons)));

                const double perCluster = electrons / nClusters;
                const double energyPerCluster = hit.energyDeposit / nClusters;

                fTotalElectrons += electrons;

                for (int i = 0; i < nClusters; ++i) {
                    const double x = hit.x + sigmaT    * unitGauss(fRng);
                    const double y = hit.y + sigmaT    * unitGauss(fRng);
                    const double t = driftTime + sigmaTime * unitGauss(fRng);
                    if (t < 0) continue;

                    CollectElectrons(plane, x, y, t, perCluster,
                                     energyPerCluster, particle.trackID, pads);
                }
            }
        };

        driftHits(particle.tpcHits);
        driftHits(particle.sec_tpcHits);
    }
}

// ---------------------------------------------------------------------------
// Pad response
// ---------------------------------------------------------------------------

void TPCDigiModule::CollectElectrons(int plane, double x, double y, double t,
                                     double electrons, double energy, int trackID,
                                     PadMap& pads)
{
    const int tick = static_cast<int>(t / fConditions.samplePeriod);
    if (tick < 0 || tick >= fMaxTicks) return;

    // Pad the arrival point falls on, and the neighbours that share its charge.
    // The pad itself may not exist -- the point can land just outside the
    // circle and still be collected by a neighbour inside it.
    const int row0  = fReadout->RowOf(x);
    const int col0  = fReadout->ColOf(y);
    const int range = fReadout->ResponseRange();

    int    nPads = 0;
    double sumWeight = 0.0;

    for (int dr = -range; dr <= range; ++dr) {
        for (int dc = -range; dc <= range; ++dc) {
            const int row = row0 + dr;
            const int col = col0 + dc;
            if (!fReadout->IsValidPad(row, col)) continue;

            const double w = fReadout->PadResponse(x - fReadout->PadX(row),
                                                   y - fReadout->PadY(col));
            if (w <= 0) continue;

            fPadWeights[nPads] = w;
            fPadRows[nPads]    = row;
            fPadCols[nPads]    = col;
            sumWeight         += w;
            nPads++;
        }
    }

    if (nPads == 0 || sumWeight <= 0) return;   // Electron missed the readout

    // Normalized, the response only decides how the charge is shared between
    // pads; unnormalized it also rejects charge landing between them, which is
    // the behaviour of the original parametrization
    const double norm = fNormalizeResponse ? 1.0 / sumWeight : 1.0;

    for (int i = 0; i < nPads; ++i) {
        const int channel = fReadout->Channel(plane, fPadRows[i], fPadCols[i]);
        if (channel < 0) continue;

        PadSignal& pad = pads[channel];
        pad.plane = plane;
        pad.row   = fPadRows[i];
        pad.col   = fPadCols[i];
        pad.ticks[tick].Add(electrons * fPadWeights[i] * norm,
                            energy    * fPadWeights[i] * norm, trackID);
    }
}

void TPCDigiModule::TickSignal::Add(float ne, float e, int trackID)
{
    electrons += ne;
    energy    += e;
    for (auto& contribution : trackEnergy) {
        if (contribution.first == trackID) {
            contribution.second += e;
            return;
        }
    }
    trackEnergy.emplace_back(trackID, e);
}

// ---------------------------------------------------------------------------
// Waveform digitization
// ---------------------------------------------------------------------------

void TPCDigiModule::DigitizePads(const PadMap& pads)
{
    const int nResponse = static_cast<int>(fResponse.size());

    std::normal_distribution<double> noise(0.0, fNoiseLevel * fConditions.gain);

    std::vector<float> adc, trueEnergy;

    for (const auto& entry : pads) {
        const PadSignal& pad = entry.second;
        if (pad.ticks.empty()) continue;

        // Window covering the shaped signal, plus room for the ROI padding
        const int firstTick = pad.ticks.begin()->first;
        const int lastTick  = pad.ticks.rbegin()->first;

        const int tickOffset = std::max(0, firstTick - fROIPadding);
        const int nSamples   = (lastTick + nResponse + fROIPadding) - tickOffset;
        if (nSamples <= 0) continue;

        adc.assign(nSamples, 0.0f);
        trueEnergy.assign(nSamples, 0.0f);

        // Convolve the collected charge with the shaping function. The true
        // energy is left at the tick the charge arrived in, so that summing it
        // over a pulse gives back the energy that made the pulse.
        for (const auto& tickEntry : pad.ticks) {
            const int   i0 = tickEntry.first - tickOffset;
            const float ne = tickEntry.second.electrons;

            for (int k = 0; k < nResponse; ++k) {
                adc[i0 + k] += ne * fResponse[k] * fConditions.gain;
            }
            trueEnergy[i0] += tickEntry.second.energy;
        }

        // Pedestal, electronics noise, ADC range and quantization
        for (int i = 0; i < nSamples; ++i) {
            double value = fConditions.pedestal + adc[i] + noise(fRng);
            value = std::min(std::max(value, 0.0), fConditions.adcRange - 1.0);
            adc[i] = std::floor(value);
        }

        ExtractROIs(pad, adc, trueEnergy, tickOffset);
    }
}

// ---------------------------------------------------------------------------
// Zero suppression
// ---------------------------------------------------------------------------

void TPCDigiModule::ExtractROIs(const PadSignal& pad,
                                const std::vector<float>& adc,
                                const std::vector<float>& trueEnergy,
                                int tickOffset)
{
    const int nSamples = static_cast<int>(adc.size());
    const double threshold = fConditions.pedestal + fADCThreshold;

    int i = 0;
    while (i < nSamples) {
        if (adc[i] < threshold) { ++i; continue; }

        // Grow the region until the samples have been below threshold for
        // longer than the padding, so that nearby pulses stay in one waveform
        int begin = i;
        int end   = i;
        int j     = i;
        while (j < nSamples && j - end <= fROIPadding) {
            if (adc[j] >= threshold) end = j;
            ++j;
        }
        i = j;

        begin = std::max(0, begin - fROIPadding);
        end   = std::min(nSamples - 1, end + fROIPadding);

        digi::TPCWaveform wf;
        wf.channel   = fReadout->Channel(pad.plane, pad.row, pad.col);
        wf.plane     = pad.plane;
        wf.row       = pad.row;
        wf.col       = pad.col;
        wf.padX      = fReadout->PadX(pad.row);
        wf.padY      = fReadout->PadY(pad.col);
        wf.tickStart = tickOffset + begin;
        wf.pedestal  = fConditions.pedestal;

        wf.adc.reserve(end - begin + 1);
        for (int k = begin; k <= end; ++k) {
            wf.adc.push_back(static_cast<Short_t>(adc[k]));
        }

        if (fSaveTruth) {
            wf.trueEnergy.assign(trueEnergy.begin() + begin,
                                 trueEnergy.begin() + end + 1);

            // Track contributions from the deposits inside this region
            std::vector<std::pair<int, float>> tracks;
            float total = 0;
            auto first = pad.ticks.lower_bound(wf.tickStart);
            auto last  = pad.ticks.upper_bound(tickOffset + end);
            for (auto it = first; it != last; ++it) {
                for (const auto& contribution : it->second.trackEnergy) {
                    total += contribution.second;
                    auto found = std::find_if(tracks.begin(), tracks.end(),
                        [&](const std::pair<int, float>& t) {
                            return t.first == contribution.first; });
                    if (found != tracks.end()) found->second += contribution.second;
                    else tracks.push_back(contribution);
                }
            }
            for (const auto& track : tracks) {
                wf.trackIDs.push_back(track.first);
                wf.trackFractions.push_back(total > 0 ? track.second / total : 0);
            }
        }

        fWaveforms->push_back(std::move(wf));
    }
}
