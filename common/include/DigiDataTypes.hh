//
// DigiDataTypes.hh - Digitized hit data classes
//
// Output types produced by digitization reconstruction modules.
// These are ROOT-serializable and stored in the RecoStore for downstream
// modules (e.g. cluster finding). They are also written as TTree branches.
//
// TileDigiHit  - ECal, single HG tile readout (single SiPM)
// StripDigiHit - ECal, LG strip readout (dual SiPM with position reconstruction)
// TPCWaveform  - TPC, zero-suppressed ADC waveform of one pad channel
// TPCHit       - TPC, reconstructed hit (one pulse on one pad channel)
// TPCCluster   - TPC, group of TPCHits belonging to the same energy deposit
//

#ifndef DigiDataTypes_hh
#define DigiDataTypes_hh

#include <vector>
#include "TObject.h"

namespace digi {

    // Digitized hit from one HG tile cell (single SiPM readout)
    class TileDigiHit : public TObject
    {
    public:
        TileDigiHit()
            : segment(0), layer(0), row(0), col(0),
              adcSum(0), time(0), trueEnergy(0) {}

        Int_t   segment;       // Barrel segment or endcap identifier
        Int_t   layer;         // Layer index
        Int_t   row;           // Tile row index
        Int_t   col;           // Tile column index
        Float_t adcSum;        // Integrated ADC counts
        Float_t time;          // ADC-weighted mean time [ns]
        Float_t trueEnergy;    // True deposited energy [MeV] (MC truth)

        // Track contributions (parallel vectors)
        std::vector<Int_t>   trackIDs;       // Contributing track IDs
        std::vector<Float_t> trackFractions; // Energy fraction per track

        ClassDef(TileDigiHit, 1)
    };

    // Digitized hit from one LG strip (dual SiPM readout)
    class StripDigiHit : public TObject
    {
    public:
        StripDigiHit()
            : segment(0), layer(0), row(0),
              adcLeft(0), adcRight(0), adcCombined(0),
              timeLeft(0), timeRight(0), recoTime(0), recoPosition(0),
              stripLength(0), trueEnergy(0) {}

        Int_t   segment;        // Barrel segment or endcap identifier
        Int_t   layer;          // Layer index
        Int_t   row;            // Strip index

        Float_t adcLeft;        // ADC from left SiPM
        Float_t adcRight;       // ADC from right SiPM
        Float_t adcCombined;    // Geometric mean of left and right ADC
        Float_t timeLeft;       // Hit time at left SiPM [ns]
        Float_t timeRight;      // Hit time at right SiPM [ns]
        Float_t recoTime;       // Reconstructed interaction time [ns]
        Float_t recoPosition;   // Reconstructed position along strip [cm]
        Float_t stripLength;    // Physical strip length [cm]
        Float_t trueEnergy;     // True deposited energy [MeV] (MC truth)

        std::vector<Int_t>   trackIDs;
        std::vector<Float_t> trackFractions;

        ClassDef(StripDigiHit, 1)
    };

    // ------------------------------------------------------------------
    // TPC
    // ------------------------------------------------------------------

    // Zero-suppressed ADC waveform for one pad channel.
    //
    // One object per region of interest, so a channel hit at two well
    // separated drift times yields two waveforms. Sample i of `adc`
    // corresponds to tick (tickStart + i), i.e. to a drift time of
    // (tickStart + i + 0.5) * samplePeriod.
    class TPCWaveform : public TObject
    {
    public:
        TPCWaveform()
            : channel(0), plane(0), row(0), col(0),
              padX(0), padY(0), tickStart(0), pedestal(0) {}

        Int_t   channel;        // Global channel number
        Int_t   plane;          // Readout plane: 0 = +z end, 1 = -z end
        Int_t   row;            // Pad row index (along x)
        Int_t   col;            // Pad column index (along y)
        Float_t padX;           // Pad centre x [cm]
        Float_t padY;           // Pad centre y [cm]

        Int_t   tickStart;      // Tick of the first sample
        Float_t pedestal;       // Baseline the samples sit on [ADC]

        std::vector<Short_t> adc;         // ADC samples, pedestal included

        // MC truth (empty when the module runs with saveTruth off). The energy
        // recorded is that of the deposits whose ionization this pad collected,
        // booked at the tick the charge arrived in rather than spread over the
        // shaped pulse, so summing it over a pulse gives the energy behind it.
        std::vector<Float_t> trueEnergy;  // [MeV] per sample
        std::vector<Int_t>   trackIDs;        // Contributing track IDs
        std::vector<Float_t> trackFractions;  // Energy fraction per track

        ClassDef(TPCWaveform, 1)
    };

    // Reconstructed hit: one charge pulse found on one pad channel
    class TPCHit : public TObject
    {
    public:
        TPCHit()
            : channel(0), plane(0), row(0), col(0),
              x(0), y(0), z(0), driftTime(0),
              peakAmplitude(0), adcSum(0), charge(0), energy(0),
              sigmaT(0), sigmaZ(0), startTick(0), endTick(0),
              trueEnergy(0) {}

        Int_t   channel;        // Global channel number
        Int_t   plane;          // Readout plane: 0 = +z end, 1 = -z end
        Int_t   row;            // Pad row index (along x)
        Int_t   col;            // Pad column index (along y)

        Float_t x, y;           // Pad centre [cm]
        Float_t z;              // Drift coordinate from the pulse time [cm]
        Float_t driftTime;      // Charge-weighted drift time [us]

        Float_t peakAmplitude;  // Pulse height above pedestal [ADC]
        Float_t adcSum;         // Pedestal-subtracted integral [ADC]
        Float_t charge;         // Collected charge [electrons]
        Float_t energy;         // Calibrated energy [MeV]

        Float_t sigmaT;         // Pulse RMS in time [us]
        Float_t sigmaZ;         // Pulse RMS along the drift direction [cm]

        Int_t   startTick;      // First tick of the pulse
        Int_t   endTick;        // Last tick of the pulse (inclusive)

        Float_t trueEnergy;     // True deposited energy [MeV] (MC truth)
        std::vector<Int_t>   trackIDs;
        std::vector<Float_t> trackFractions;

        ClassDef(TPCHit, 1)
    };

    // Group of TPCHits, e.g. the hits left by one track segment
    class TPCCluster : public TObject
    {
    public:
        TPCCluster()
            : clusterID(0), plane(-1), nHits(0),
              x(0), y(0), z(0), charge(0), energy(0),
              rmsX(0), rmsY(0), rmsZ(0),
              dirX(0), dirY(0), dirZ(0), length(0), width(0),
              startX(0), startY(0), startZ(0),
              endX(0), endY(0), endZ(0),
              trueEnergy(0) {}

        Int_t   clusterID;
        Int_t   plane;          // Readout plane, or -1 if the hits span both
        Int_t   nHits;

        Float_t x, y, z;        // Charge-weighted centroid [cm]
        Float_t charge;         // Summed hit charge [electrons]
        Float_t energy;         // Summed hit energy [MeV]
        Float_t rmsX, rmsY, rmsZ;  // Charge-weighted spread about the centroid [cm]

        Float_t dirX, dirY, dirZ;  // Principal axis (unit vector)
        Float_t length;         // Extent along the principal axis [cm]
        Float_t width;          // RMS extent transverse to it [cm]

        Float_t startX, startY, startZ;  // Extreme hits along the principal axis [cm]
        Float_t endX, endY, endZ;

        std::vector<Int_t> hitIndices;   // Indices into the TPCHit collection

        Float_t trueEnergy;
        std::vector<Int_t>   trackIDs;
        std::vector<Float_t> trackFractions;

        ClassDef(TPCCluster, 1)
    };

} // namespace digi

#endif
