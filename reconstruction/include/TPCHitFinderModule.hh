//
// TPCHitFinderModule.hh - Pulse finding on digitized TPC pad waveforms
//
// Reads the zero-suppressed waveforms produced by TPCDigiModule and turns each
// charge pulse into a reconstructed hit. Pulses are found by threshold
// crossing and split at the local minima between merged peaks, so that two
// deposits arriving close in time on the same pad still give two hits.
//
// The pad gives the hit its (x,y); the charge-weighted time of the pulse gives
// its z, through the drift velocity and the readout plane it was seen on. The
// pedestal-subtracted integral is converted back to collected electrons and
// then, correcting for the electron lifetime, to energy.
//
// Inputs:
//   "TPCWaveforms"        <- std::vector<digi::TPCWaveform>
//   "TPCConditions"       <- TPCConditions
// Output:
//   "TPCHits"             -> std::vector<digi::TPCHit>  (RecoStore + TTree)
//

#ifndef TPCHitFinderModule_h
#define TPCHitFinderModule_h 1

#include "RecoModule.hh"
#include "DigiDataTypes.hh"

#include <vector>

struct TPCConditions;

class TPCHitFinderModule : public RecoModule {
public:
    TPCHitFinderModule();
    virtual ~TPCHitFinderModule();

    virtual void Initialize() override;
    virtual void Execute()    override;
    virtual void Finalize()   override;

    virtual std::vector<ObjectSpec> GetInputSpec()  const override;
    virtual std::vector<ObjectSpec> GetOutputSpec() const override;

private:
    // Half-open [begin, end) range of samples making up one pulse
    struct Pulse {
        int begin = 0;
        int end   = 0;
        int peak  = 0;
    };

    std::vector<Pulse> FindPulses(const std::vector<float>& signal) const;
    digi::TPCHit       MakeHit(const digi::TPCWaveform& wf,
                               const std::vector<float>& signal,
                               const Pulse& pulse) const;

    // --- I/O ---
    const std::vector<digi::TPCWaveform>* fWaveforms;
    std::vector<digi::TPCHit>*            fHits;
    const TPCConditions*                  fConditions;

    // --- Parameters ---
    double fThreshold;          // Pulse finding threshold above pedestal [ADC]
    int    fMinPulseTicks;      // Shortest pulse accepted
    double fMinPulseADC;        // Smallest integral accepted [ADC]
    double fValleyFraction;     // Local minimum this far below the lower of two
                                // neighbouring peaks splits them
    double fChargeCalibration;  // Multiplies the reconstructed charge
    bool   fLifetimeCorrection;

    // --- Running totals for the end-of-job summary ---
    long   fNEvents, fNHits;
    double fTotalCharge, fTotalEnergy;
};

#endif
