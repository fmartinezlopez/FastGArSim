 /***************************************************************************
 * PlotStyle.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Common plotting style and small drawing helpers shared by the analysis
 *   macros.
 *
 ***************************************************************************/

#ifndef PlotStyle_hh
#define PlotStyle_hh

#include "Rtypes.h"

class TH1;

namespace ana {

// Apply the standard FastGArSim plotting style to gStyle
void SetPlotStyle();

// Draw two histograms on the same canvas, area-normalised, and save it to
// `outputName` (any format understood by TCanvas::SaveAs).
//   labelA / labelB -- legend entries
//   logY            -- use a logarithmic vertical scale
// The histograms are scaled in place.
void DrawComparison(TH1* histA, TH1* histB,
                    const char* labelA, const char* labelB,
                    const char* outputName,
                    Bool_t logY = kTRUE);

} // namespace ana

#endif
