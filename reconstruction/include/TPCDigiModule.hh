//
// TPCDigiModule.hh - TPC readout simulation
//
// Reads raw root::TPCHit energy deposits from the detector_simulation "Events"
// TTree and turns them into digitized pad waveforms, following the chain
//
//   energy deposit -> ionization electrons (W value, Fano fluctuations)
//                  -> drift to the readout plane (attenuation by the electron
//                     lifetime, transverse and longitudinal diffusion)
//                  -> collection on pads (generalized-Gaussian pad response)
//                  -> shaped, amplified and digitized waveform (CR-RC^2
//                     shaping, gain, electronics noise, ADC quantization)
//                  -> zero suppression into regions of interest
//
// Electrons are not drifted one by one: the electrons from one deposit are
// grouped into macro-electron clusters, each of which is diffused as a unit.
//
// Outputs:
//   "TPCWaveforms"        -> std::vector<digi::TPCWaveform>  (RecoStore + TTree)
//   "TPCReadoutGeometry"  -> TPCReadoutGeometry              (RecoStore)
//   "TPCConditions"       -> TPCConditions                   (RecoStore)
//

#ifndef TPCDigiModule_h
#define TPCDigiModule_h 1

#include "RecoModule.hh"
#include "DigiDataTypes.hh"
#include "TPCConditions.hh"
#include "TPCReadoutGeometry.hh"

#include <map>
#include <random>
#include <unordered_map>
#include <vector>

namespace root { class Event; }

class TPCDigiModule : public RecoModule {
public:
    TPCDigiModule();
    virtual ~TPCDigiModule();

    virtual void Initialize() override;
    virtual void Execute()    override;
    virtual void Finalize()   override;

    virtual std::vector<ObjectSpec> GetInputSpec()  const override;
    virtual std::vector<ObjectSpec> GetOutputSpec() const override;

private:
    // Charge and truth collected by one pad in one time bin
    struct TickSignal {
        float electrons = 0;
        float energy    = 0;   // True energy behind that charge [MeV]
        std::vector<std::pair<int, float>> trackEnergy;  // Usually one or two entries

        void Add(float ne, float e, int trackID);
    };

    // A pad that saw charge: its ticks, in time order
    struct PadSignal {
        int plane = 0, row = 0, col = 0;
        std::map<int, TickSignal> ticks;
    };

    using PadMap = std::unordered_map<int, PadSignal>;   // channel -> signal

    // --- Pipeline steps ---
    void   BuildResponseFunction();
    double IonizationElectrons(double energyDeposit);
    void   DriftDeposits(const root::Event* event, PadMap& pads);
    void   CollectElectrons(int plane, double x, double y, double t,
                            double electrons, double energy, int trackID,
                            PadMap& pads);
    void   DigitizePads(const PadMap& pads);
    void   ExtractROIs(const PadSignal& pad,
                       const std::vector<float>& adc,
                       const std::vector<float>& trueEnergy,
                       int tickOffset);

    // --- I/O ---
    root::Event*                     fSimEvent;    // Branch address on the input tree
    std::vector<digi::TPCWaveform>*  fWaveforms;

    // --- Detector description, shared with downstream modules ---
    TPCReadoutGeometry* fReadout;
    TPCConditions       fConditions;

    // --- Parameters ---
    double fPadPitch;
    double fResponseWidth, fResponseShape;
    int    fResponseRange;
    bool   fNormalizeResponse;
    double fElectronsPerCluster;
    int    fMinClusters;
    double fNoiseLevel;         // Electronics noise [electrons]
    double fADCThreshold;       // Zero suppression threshold above pedestal [ADC]
    int    fROIPadding;         // Ticks kept either side of a threshold crossing
    int    fMaxTicks;
    bool   fSaveTruth, fWriteWaveforms;

    // --- Shaping function, sampled at the tick spacing ---
    std::vector<double> fResponse;

    // Scratch space for the pads sharing one macro-electron's charge
    std::vector<double> fPadWeights;
    std::vector<int>    fPadRows, fPadCols;

    std::mt19937 fRng;

    // --- Running totals for the end-of-job summary ---
    long   fNEvents, fNDeposits, fNWaveforms;
    double fTotalElectrons;
};

#endif
