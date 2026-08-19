 /***************************************************************************
 * AnalysisMath.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the shared numerical helpers.
 *
 ***************************************************************************/

#include "AnalysisMath.hh"

#include <algorithm>
#include <cmath>

namespace ana {

Double_t TruncatedMean(std::vector<Double_t>& values, Double_t fraction, size_t* nKept)
{
    if (nKept) *nKept = 0;
    if (values.empty()) return 0.;

    if (fraction <= 0.) fraction = 0.;
    if (fraction > 1.) fraction = 1.;

    // Keep at least one sample, so a non-empty input always gives a value
    size_t keep = static_cast<size_t>(std::floor(values.size() * fraction));
    if (keep < 1) keep = 1;
    if (keep > values.size()) keep = values.size();

    // Only the `keep` smallest elements need to be in order
    std::partial_sort(values.begin(), values.begin() + keep, values.end());

    Double_t sum = 0.;
    for (size_t i = 0; i < keep; ++i) sum += values[i];

    if (nKept) *nKept = keep;
    return sum / keep;
}

std::vector<Double_t> LogSpacedBins(Int_t nBins, Double_t min, Double_t max)
{
    std::vector<Double_t> edges;
    if (nBins < 1 || min <= 0. || max <= min) return edges;

    edges.resize(nBins + 1);

    const Double_t logMin = std::log10(min);
    const Double_t logMax = std::log10(max);
    const Double_t step = (logMax - logMin) / nBins;

    for (Int_t i = 0; i <= nBins; ++i) {
        edges[i] = std::pow(10., logMin + i * step);
    }

    return edges;
}

} // namespace ana
