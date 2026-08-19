//
// DigiDataTypes.hh - Digitized hit data classes
//
// Output types produced by digitization reconstruction modules.
// These are ROOT-serializable and stored in the RecoStore for downstream
// modules (e.g. cluster finding). They are also written as TTree branches.
//
// TileDigiHit  - single HG tile readout (single SiPM)
// StripDigiHit - LG strip readout (dual SiPM with position reconstruction)
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

} // namespace digi

#endif
