//
// ECalDigiAnalysis.C - Example ROOT macro to plot ECal digitization output
//
// Reads the output of GArReconstruction run with the ecal_digi.mac macro.
//
// Usage (from the build directory):
//   root -l 'ECalDigiAnalysis.C("ecal_digi.root")'
//
// Needs the DigiDataDict shared library and the common/ headers. Both come
// from the build environment:
//   source <build>/setup.sh
//

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "DigiDataTypes.hh"

void ECalDigiAnalysis(const char* filename = "ecal_digi.root") {

    std::cout << "\n==================================================" << std::endl;
    std::cout << "   ECal Digitization Analysis" << std::endl;
    std::cout << "==================================================" << std::endl;

    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("Reco");
    if (!tree) {
        std::cerr << "Error: Cannot find the Reco tree in file" << std::endl;
        return;
    }

    std::vector<digi::TileDigiHit>*  tileHits  = nullptr;
    std::vector<digi::StripDigiHit>* stripHits = nullptr;

    tree->SetBranchAddress("ECalTileDigiHits",  &tileHits);
    tree->SetBranchAddress("ECalStripDigiHits", &stripHits);

    // -----------------------------------------------------------------------
    // Histograms - tiles (HG layers)
    // -----------------------------------------------------------------------

    // Multiplicity
    TH1F* hNTiles  = new TH1F("hNTiles",
        "HG tile hits per event;N_{tile hits};Events",
        100, 0, 500);

    // ADC and energy
    TH1F* hTileADC = new TH1F("hTileADC",
        "HG tile ADC sum;ADC counts;Hits",
        100, 0, 4096);
    TH1F* hTileE   = new TH1F("hTileE",
        "HG tile true energy;True energy [MeV];Hits",
        100, 0, 20);

    // Hit time
    TH1F* hTileTime = new TH1F("hTileTime",
        "HG tile hit time;Time [ns];Hits",
        100, 0, 200);

    // Layer occupancy
    TH1F* hTileLayer = new TH1F("hTileLayer",
        "HG tile hits vs. layer;Layer;Hits",
        30, -0.5, 29.5);

    // -----------------------------------------------------------------------
    // Histograms - strips (LG layers)
    // -----------------------------------------------------------------------

    TH1F* hNStrips = new TH1F("hNStrips",
        "LG strip hits per event;N_{strip hits};Events",
        100, 0, 500);

    TH1F* hStripADCComb = new TH1F("hStripADCComb",
        "LG strip combined ADC;ADC counts (geometric mean);Hits",
        100, 0, 100);
    TH1F* hStripE = new TH1F("hStripE",
        "LG strip true energy;True energy [MeV];Hits",
        100, 0, 20);

    TH1F* hStripRecoTime = new TH1F("hStripRecoTime",
        "LG strip reconstructed time;Reco time [ns];Hits",
        100, 0, 200);

    // Reconstructed position along strip
    TH1F* hStripRecoPos = new TH1F("hStripRecoPos",
        "LG strip reconstructed position;Position along strip [cm];Hits",
        100, -150, 150);

    // ADC left vs right scatter (with attenuation: spreads off diagonal with position)
    TH2F* hStripADCLvsR = new TH2F("hStripADCLvsR",
        "LG strip ADC left vs. right;ADC left;ADC right",
        50, 0, 100, 50, 0, 100);

    // Time difference between SiPM ends (encodes position along strip)
    TH1F* hStripTimeDiff = new TH1F("hStripTimeDiff",
        "LG strip SiPM time difference;#Deltat = t_{R} - t_{L} [ns];Hits",
        100, -50, 50);

    // ADC asymmetry vs reco position (sigmoid shape with attenuation)
    TH2F* hStripAsymVsPos = new TH2F("hStripAsymVsPos",
        "LG strip ADC asymmetry vs. reco position;"
        "Reco position [cm];(ADC_{R}-ADC_{L})/(ADC_{R}+ADC_{L})",
        50, -150, 150, 50, -1.0, 1.0);

    // Layer occupancy
    TH1F* hStripLayer = new TH1F("hStripLayer",
        "LG strip hits vs. layer;Layer;Hits",
        30, -0.5, 29.5);

    // -----------------------------------------------------------------------
    // Event-level totals
    // -----------------------------------------------------------------------

    TH1F* hTotalE = new TH1F("hTotalE",
        "Total true energy per event (ECal);True energy [MeV];Events",
        100, 0, 2000);

    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------

    Long64_t nEntries = tree->GetEntries();
    std::cout << " Processing " << nEntries << " events..." << std::endl;

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        tree->GetEntry(iEntry);

        hNTiles->Fill(tileHits->size());
        hNStrips->Fill(stripHits->size());

        float evtEnergy = 0;

        for (const auto& h : *tileHits) {
            hTileADC  ->Fill(h.adcSum);
            hTileE    ->Fill(h.trueEnergy);
            hTileTime ->Fill(h.time);
            hTileLayer->Fill(h.layer);
            evtEnergy += h.trueEnergy;
        }

        for (const auto& h : *stripHits) {
            hStripADCComb ->Fill(h.adcCombined);
            hStripE       ->Fill(h.trueEnergy);
            hStripRecoTime->Fill(h.recoTime);
            hStripRecoPos ->Fill(h.recoPosition);
            hStripLayer   ->Fill(h.layer);
            evtEnergy     += h.trueEnergy;

            hStripADCLvsR ->Fill(h.adcLeft, h.adcRight);
            hStripTimeDiff->Fill(h.timeRight - h.timeLeft);
            float adcSum = h.adcLeft + h.adcRight;
            if (adcSum > 0)
                hStripAsymVsPos->Fill(h.recoPosition,
                                      (h.adcRight - h.adcLeft) / adcSum);
        }

        hTotalE->Fill(evtEnergy);
    }

    // -----------------------------------------------------------------------
    // Draw
    // -----------------------------------------------------------------------

    gStyle->SetOptStat(1110);

    // Canvas 1: hit multiplicity and event energy
    TCanvas* c1 = new TCanvas("c1", "ECal Digi Overview", 1200, 400);
    c1->Divide(3, 1);

    c1->cd(1); hNTiles ->Draw();
    c1->cd(2); hNStrips->Draw();
    c1->cd(3); hTotalE ->Draw();

    // Canvas 2: HG tile quantities
    TCanvas* c2 = new TCanvas("c2", "HG Tile Hits", 1200, 600);
    c2->Divide(3, 2);

    c2->cd(1); hTileADC  ->Draw();
    c2->cd(2); hTileE    ->Draw();
    c2->cd(3); hTileTime ->Draw();
    c2->cd(4); hTileLayer->Draw();

    // ADC-weighted energy calibration: ADC vs true energy
    TH2F* hTileCalib = new TH2F("hTileCalib",
        "HG tile: ADC vs. true energy;True energy [MeV];ADC counts",
        50, 0, 20, 50, 0, 4096);
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        tree->GetEntry(iEntry);
        for (const auto& h : *tileHits)
            hTileCalib->Fill(h.trueEnergy, h.adcSum);
    }
    c2->cd(5); hTileCalib->Draw("colz");

    // Canvas 3: LG strip quantities
    TCanvas* c3 = new TCanvas("c3", "LG Strip Hits", 1200, 800);
    c3->Divide(3, 2);

    c3->cd(1); hStripADCComb ->Draw();
    c3->cd(2); hStripE       ->Draw();
    c3->cd(3); hStripRecoTime->Draw();
    c3->cd(4); hStripRecoPos ->Draw();
    c3->cd(5); hStripLayer   ->Draw();
    c3->cd(6); hStripADCLvsR ->Draw("colz");

    // Canvas 4: strip position reconstruction diagnostics
    TCanvas* c4 = new TCanvas("c4", "LG Strip Position Reco", 1200, 400);
    c4->Divide(3, 1);
    c4->cd(1); hStripTimeDiff  ->Draw();
    c4->cd(2); hStripRecoPos   ->Draw();
    c4->cd(3); hStripAsymVsPos ->Draw("colz");
    c4->Update();

    c1->Update();
    c2->Update();
    c3->Update();

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------

    std::cout << "\n Analysis complete." << std::endl;
    std::cout << "   Events processed : " << nEntries << std::endl;
    std::cout << "   Mean tile hits/event  : "
              << hNTiles->GetMean()  << std::endl;
    std::cout << "   Mean strip hits/event : "
              << hNStrips->GetMean() << std::endl;
    std::cout << "==================================================" << std::endl;
}
