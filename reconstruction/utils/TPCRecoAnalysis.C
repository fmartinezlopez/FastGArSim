//
// TPCRecoAnalysis.C - Example ROOT macro to plot TPC reconstruction output
//
// Reads the output of GArReconstruction run with the tpc_reco.mac macro.
//
// Usage (from the build directory):
//   root -l 'TPCRecoAnalysis.C("tpc_reco.root")'
//
// Needs the DigiDataDict shared library and the common/ headers. Both come
// from the build environment:
//   source <build>/setup.sh
//

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "DigiDataTypes.hh"

void TPCRecoAnalysis(const char* filename = "tpc_reco.root") {

    std::cout << "\n==================================================" << std::endl;
    std::cout << "   TPC Reconstruction Analysis" << std::endl;
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

    std::vector<digi::TPCHit>*     hits     = nullptr;
    std::vector<digi::TPCCluster>* clusters = nullptr;

    tree->SetBranchAddress("TPCHits", &hits);
    if (tree->GetBranch("TPCClusters")) {
        tree->SetBranchAddress("TPCClusters", &clusters);
    }

    TH1F* hNHits    = new TH1F("hNHits", "Hits per event;N_{hits};Events",
                               100, 0, 5000);
    TH1F* hHitADC   = new TH1F("hHitADC", "Hit integral;#SigmaADC;Hits",
                               100, 0, 20000);
    TH1F* hHitCharge = new TH1F("hHitCharge", "Hit charge;Electrons;Hits",
                                100, 0, 2000);
    TH1F* hHitEnergy = new TH1F("hHitEnergy",
                                "Hit energy;Reconstructed energy [MeV];Hits",
                                100, 0, 0.1);
    TH1F* hResidual = new TH1F("hResidual",
                               "Hit energy residual;(E_{reco} - E_{true}) / E_{true};Hits",
                               100, -1, 1);
    TH2F* hPadPlane = new TH2F("hPadPlane", "Hit positions;x [cm];y [cm]",
                               250, -250, 250, 250, -250, 250);
    TH2F* hDriftView = new TH2F("hDriftView", "Hit positions;z [cm];x [cm]",
                                250, -250, 250, 250, -250, 250);
    TH1F* hNClusters = new TH1F("hNClusters", "Clusters per event;N_{clusters};Events",
                                50, 0, 50);
    TH1F* hClusterEnergy = new TH1F("hClusterEnergy",
                                    "Cluster energy;Energy [MeV];Clusters",
                                    100, 0, 100);
    TH1F* hClusterLength = new TH1F("hClusterLength",
                                    "Cluster length;Length [cm];Clusters",
                                    100, 0, 200);

    const Long64_t nEntries = tree->GetEntries();
    std::cout << "\n Processing " << nEntries << " events..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        hNHits->Fill(hits->size());
        for (const auto& hit : *hits) {
            hHitADC->Fill(hit.adcSum);
            hHitCharge->Fill(hit.charge);
            hHitEnergy->Fill(hit.energy);
            hPadPlane->Fill(hit.x, hit.y);
            hDriftView->Fill(hit.z, hit.x);
            if (hit.trueEnergy > 0) {
                hResidual->Fill((hit.energy - hit.trueEnergy) / hit.trueEnergy);
            }
        }

        if (!clusters) continue;
        hNClusters->Fill(clusters->size());
        for (const auto& cluster : *clusters) {
            hClusterEnergy->Fill(cluster.energy);
            hClusterLength->Fill(cluster.length);
        }
    }

    gStyle->SetOptStat(1110);

    TCanvas* c1 = new TCanvas("c1", "TPC Reconstruction", 1400, 900);
    c1->Divide(3, 3);

    c1->cd(1); hNHits->Draw();
    c1->cd(2); hHitADC->Draw();
    c1->cd(3); hHitCharge->Draw();
    c1->cd(4); hHitEnergy->Draw();
    c1->cd(5); hResidual->Draw();
    c1->cd(6); hPadPlane->Draw("COLZ");
    c1->cd(7); hDriftView->Draw("COLZ");
    c1->cd(8); hNClusters->Draw();
    c1->cd(9); hClusterLength->Draw();

    c1->Update();
    c1->SaveAs("tpc_reco.png");

    std::cout << "\n Analysis complete!" << std::endl;
    std::cout << " Total events processed: " << nEntries << std::endl;
    std::cout << "\n==================================================" << std::endl;
}
