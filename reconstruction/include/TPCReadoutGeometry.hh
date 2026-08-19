//
// TPCReadoutGeometry.hh - Pad plane description for the HPgTPC readout
//
// The TPC is a cylinder centred on the origin with its axis along z, so the
// readout planes are the two circular end faces at z = +-L/2 and the drift
// direction is z. Each plane is tiled with square pads of side `pitch`, laid
// out on a regular grid in (x,y); pads whose centre falls outside the TPC
// radius do not exist.
//
// Pads are addressed by (row, col), the grid indices along x and y, and are
// numbered channel 0 .. NChannels()-1, plane by plane and row by row. Channel
// numbers are contiguous: the gaps left by the pads outside the circle are
// skipped rather than numbered.
//
// The class also carries the pad response function: the probability that an
// electron arriving at a given point on the plane is collected by a given pad,
// modelled by the generalized ("squared") 2D Gaussian
//
//     PRF(dx, dy) = exp( -(|dx|^p + |dy|^p) / (2 sigma^p) ),
//     sigma       = pitch / responseWidth
//
// with dx, dy the offsets from the pad centre. With the default p = 4 this is
// a flat-topped, steep-sided function roughly the size of one pad.
//

#ifndef TPCReadoutGeometry_h
#define TPCReadoutGeometry_h 1

#include <vector>

class TPCReadoutGeometry {
public:
    // radius, pitch in cm; nPlanes is 1 (single-sided) or 2 (double-sided)
    TPCReadoutGeometry(double radius, double pitch, int nPlanes,
                       double responseWidth = 2.5, double responseShape = 4.0,
                       int responseRange = 1);

    // --- Layout ---
    double Radius()  const { return fRadius; }
    double Pitch()   const { return fPitch;  }
    int    NPlanes() const { return fNPlanes; }
    int    NRows()   const { return fNRows;  }
    int    NPadsPerPlane() const { return fNPadsPerPlane; }
    int    NChannels()     const { return fNPadsPerPlane * fNPlanes; }

    // Grid index range: rows and columns run over [-NRows()/2, +NRows()/2]
    int RowMin() const { return -fHalfRows; }
    int RowMax() const { return  fHalfRows; }

    bool IsValidPad(int row, int col) const;

    // Pad centre [cm]. Undefined for invalid pads; check with IsValidPad first.
    double PadX(int row) const { return row * fPitch; }
    double PadY(int col) const { return col * fPitch; }

    // Grid indices of the pad a coordinate falls on, whether or not that pad
    // exists -- FindPad checks that it does
    int RowOf(double x) const;
    int ColOf(double y) const;

    // Grid indices of the pad containing (x,y). Returns false if the point
    // falls outside the instrumented area.
    bool FindPad(double x, double y, int& row, int& col) const;

    // Global channel number, or -1 for an invalid pad
    int Channel(int plane, int row, int col) const;

    // --- Pad response ---
    // Number of pad rows/columns either side of the pad the electron landed on
    // that are considered as charge collectors
    int ResponseRange() const { return fResponseRange; }

    // Collection weight of a pad whose centre is (dx,dy) away from the
    // electron's arrival point [cm]
    double PadResponse(double dx, double dy) const;

    void Print() const;

private:
    double fRadius, fPitch;
    int    fNPlanes;
    int    fNRows, fHalfRows;
    int    fNPadsPerPlane;

    double fResponseSigma, fResponseShape;
    int    fResponseRange;

    // Per row: valid column range and the channel number of its first pad.
    // Indexed by (row - RowMin()).
    std::vector<int> fColMin, fColMax, fRowOffset;
};

#endif
