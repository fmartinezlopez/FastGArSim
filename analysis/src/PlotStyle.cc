 /***************************************************************************
 * PlotStyle.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the shared plotting style and drawing helpers.
 *
 ***************************************************************************/

#include "PlotStyle.hh"

#include <algorithm>

#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"

namespace ana {

void SetPlotStyle()
{
    // General plotting options
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // Set canvas margins
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.08);

    // Title and label sizes
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetTitleOffset(1.2, "Y");

    // Use better fonts
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XY");
    gStyle->SetTitleFont(42, "XY");
}

void DrawComparison(TH1* histA, TH1* histB,
                    const char* labelA, const char* labelB,
                    const char* outputName,
                    Bool_t logY)
{
    if (!histA || !histB) return;

    TCanvas canvas("canvas", outputName, 800, 600);
    if (logY) canvas.SetLogy();

    // Area-normalise both histograms so shapes can be compared
    if (histA->Integral() > 0) histA->Scale(1.0 / histA->Integral(), "width");
    if (histB->Integral() > 0) histB->Scale(1.0 / histB->Integral(), "width");

    histA->SetLineColor(kBlue);
    histA->SetLineWidth(2);
    histA->SetFillColor(kBlue - 10);
    histA->Draw("hist");

    histB->SetLineColor(kRed);
    histB->SetLineWidth(2);
    histB->SetFillColor(kRed);
    histB->SetFillStyle(3004);
    histB->Draw("hist same");

    // Make sure neither distribution is clipped
    const Double_t maximum = std::max(histA->GetMaximum(), histB->GetMaximum());
    histA->SetMaximum(logY ? maximum * 5.0 : maximum * 1.3);

    TLegend legend(0.75, 0.75, 0.90, 0.90);  // x1, y1, x2, y2 in NDC
    legend.SetBorderSize(0);
    legend.SetTextSize(0.050);
    legend.AddEntry(histA, labelA, "f");
    legend.AddEntry(histB, labelB, "f");
    legend.Draw();

    canvas.SaveAs(outputName);
}

} // namespace ana
