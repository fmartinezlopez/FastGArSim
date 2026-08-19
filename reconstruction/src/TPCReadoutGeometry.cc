//
// TPCReadoutGeometry.cc - Pad plane description for the HPgTPC readout
//

#include "TPCReadoutGeometry.hh"

#include <cmath>
#include <iostream>

TPCReadoutGeometry::TPCReadoutGeometry(double radius, double pitch, int nPlanes,
                                       double responseWidth, double responseShape,
                                       int responseRange)
    : fRadius(radius), fPitch(pitch), fNPlanes(nPlanes),
      fNRows(0), fHalfRows(0), fNPadsPerPlane(0),
      fResponseSigma(pitch / responseWidth), fResponseShape(responseShape),
      fResponseRange(responseRange)
{
    // Odd number of rows, so that one pad is centred on the beam axis
    fHalfRows = static_cast<int>(std::floor(fRadius / fPitch));
    fNRows    = 2 * fHalfRows + 1;

    fColMin.resize(fNRows);
    fColMax.resize(fNRows);
    fRowOffset.resize(fNRows);

    // A pad exists if its centre is inside the circle, so for a given row the
    // valid columns are a contiguous band symmetric about zero
    for (int i = 0; i < fNRows; ++i) {
        const int    row = i - fHalfRows;
        const double x   = PadX(row);
        const double y2  = fRadius * fRadius - x * x;

        fRowOffset[i] = fNPadsPerPlane;
        if (y2 <= 0) {                 // row lies outside the circle entirely
            fColMin[i] = 1;
            fColMax[i] = 0;
            continue;
        }
        fColMax[i] = static_cast<int>(std::floor(std::sqrt(y2) / fPitch));
        fColMin[i] = -fColMax[i];
        fNPadsPerPlane += fColMax[i] - fColMin[i] + 1;
    }
}

bool TPCReadoutGeometry::IsValidPad(int row, int col) const
{
    if (row < RowMin() || row > RowMax()) return false;
    const int i = row - RowMin();
    return col >= fColMin[i] && col <= fColMax[i];
}

int TPCReadoutGeometry::RowOf(double x) const
{
    return static_cast<int>(std::floor(x / fPitch + 0.5));
}

int TPCReadoutGeometry::ColOf(double y) const
{
    return static_cast<int>(std::floor(y / fPitch + 0.5));
}

bool TPCReadoutGeometry::FindPad(double x, double y, int& row, int& col) const
{
    row = RowOf(x);
    col = ColOf(y);
    return IsValidPad(row, col);
}

int TPCReadoutGeometry::Channel(int plane, int row, int col) const
{
    if (plane < 0 || plane >= fNPlanes)  return -1;
    if (!IsValidPad(row, col))           return -1;

    const int i = row - RowMin();
    return plane * fNPadsPerPlane + fRowOffset[i] + (col - fColMin[i]);
}

double TPCReadoutGeometry::PadResponse(double dx, double dy) const
{
    const double s = std::pow(std::fabs(dx), fResponseShape)
                   + std::pow(std::fabs(dy), fResponseShape);
    return std::exp(-s / (2.0 * std::pow(fResponseSigma, fResponseShape)));
}

void TPCReadoutGeometry::Print() const
{
    std::cout << "   Pad plane radius [cm]: " << fRadius << "\n"
              << "   Pad pitch [cm]: "        << fPitch  << "\n"
              << "   Readout planes: "        << fNPlanes << "\n"
              << "   Pad grid: "              << fNRows << " x " << fNRows << "\n"
              << "   Pads per plane: "        << fNPadsPerPlane << "\n"
              << "   Total channels: "        << NChannels() << "\n"
              << "   Pad response sigma [cm]: " << fResponseSigma
              << " (shape " << fResponseShape
              << ", range +-" << fResponseRange << " pads)" << std::endl;
}
