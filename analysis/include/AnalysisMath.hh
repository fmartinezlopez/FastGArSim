 /***************************************************************************
 * AnalysisMath.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Small numerical helpers shared by the analysis macros.
 *
 ***************************************************************************/

#ifndef AnalysisMath_hh
#define AnalysisMath_hh

#include <cstddef>
#include <vector>

#include "Rtypes.h"

namespace ana {

// Mean of the lowest `fraction` of `values` -- the standard truncated mean
// used to tame the Landau tail of ionisation samples.
//
//   fraction -- fraction of the samples to keep, in (0, 1]. 0.6 keeps the
//               lowest 60%.
//   nKept    -- if given, receives the number of samples actually averaged.
//
// The number of samples kept is floor(size * fraction), with a floor of one
// sample so that a non-empty input always yields a value. Returns 0 for an
// empty input.
//
// NOTE: `values` is sorted in place, so pass a copy if the original order
// matters.
Double_t TruncatedMean(std::vector<Double_t>& values,
                       Double_t fraction,
                       size_t* nKept = nullptr);

// Logarithmically spaced bin edges, for use with the TH1/TH2 constructors
// that take an explicit edge array. Returns nBins + 1 edges spanning
// [min, max]; both must be strictly positive.
std::vector<Double_t> LogSpacedBins(Int_t nBins, Double_t min, Double_t max);

} // namespace ana

#endif
