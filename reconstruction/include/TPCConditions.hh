//
// TPCConditions.hh - Drift and electronics constants of the TPC response
//
// TPCDigiModule fills this and registers it in the RecoStore under
// "TPCConditions"; TPCHitFinderModule reads it back, so the two modules
// cannot disagree about the sampling period, the drift velocity or the charge
// calibration without the macro saying so explicitly.
//
// The TPC gas volume spans z in [-halfLength, +halfLength]. A double-sided
// readout has a cathode at z = 0 and anodes at both ends, so the maximum drift
// distance is halfLength; a single-sided one has the cathode at the far end
// and drifts the full length to the anode at readoutSide * halfLength.
//

#ifndef TPCConditions_h
#define TPCConditions_h 1

struct TPCConditions {
    // --- Drift ---
    double driftVelocity    = 3.011;   // [cm/us]
    double electronLifetime = 3.0e6;   // [us]
    double diffusionT       = 0.0160;  // Transverse diffusion   [cm/sqrt(cm)]
    double diffusionL       = 0.0201;  // Longitudinal diffusion [cm/sqrt(cm)]

    // --- Ionization ---
    double wValue           = 26.4;    // [eV] per electron-ion pair
    double fanoFactor       = 0.16;
    double collectionEff    = 0.7;     // Fraction of electrons reaching the pads

    // --- Electronics ---
    double samplePeriod     = 0.05;    // [us] (1 / sampling rate)
    double shapingTime      = 0.1;     // [us]
    double gain             = 0.5;     // Peak ADC per collected electron
    double pedestal         = 100.0;   // [ADC]
    double adcRange         = 4096.0;  // [ADC]

    // Properties of the shaping function the digitizer used, so that the hit
    // finder can undo it: the ADC integral one collected electron produces,
    // and the mean and RMS of the shaping function in time.
    double adcPerElectron   = 1.0;
    double responseCentroid = 0.0;    // [us] after the charge arrives
    double responseSigma    = 0.0;    // [us]

    // --- Drift geometry ---
    double halfLength       = 250.0;   // Half length of the gas volume [cm]
    int    nPlanes          = 2;       // 1 = single-sided, 2 = double-sided
    int    readoutSide      = 1;       // Single-sided: +1 reads out at +z, -1 at -z

    // Which end of the TPC plane `plane` sits at: +1 for the +z end, -1 for -z
    double PlaneSide(int plane) const {
        if (nPlanes == 2) return (plane == 0) ? +1.0 : -1.0;
        return (readoutSide < 0) ? -1.0 : +1.0;
    }

    double PlaneZ(int plane)  const { return PlaneSide(plane) * halfLength; }

    // Longest drift a deposit can have before reaching a readout plane
    double MaxDrift() const { return (nPlanes == 2) ? halfLength : 2.0 * halfLength; }

    // Drift distance of a deposit at z, and the z it is reconstructed back to
    double DriftDistance(int plane, double z) const {
        return PlaneSide(plane) * (PlaneZ(plane) - z);
    }
    double DriftToZ(int plane, double driftTime) const {
        return PlaneZ(plane) - PlaneSide(plane) * driftVelocity * driftTime;
    }

    // Plane a deposit at z drifts to
    int PlaneOf(double z) const {
        if (nPlanes == 2) return (z >= 0.0) ? 0 : 1;
        return 0;
    }
};

#endif
