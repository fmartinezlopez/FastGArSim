//
// RecoExample.C - Example ROOT macro to analyze reconstruction output
//

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

// Include the reconstruction data structures
#include "../include/RecoDataTypes.hh"

void RecoExample(const char* filename = "reconstruction_output.root") {

    std::cout << "\n==================================================" << std::endl;
    std::cout << "   Reconstruction Analysis Example" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Open the reconstruction file
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        return;
    }

    // Get the reconstruction tree
    TTree* tree = (TTree*)file->Get("Reco");
    if (!tree) {
        std::cerr << "Error: Cannot find the Reco tree in file" << std::endl;
        return;
    }

    // Set up branches
    RecoEvent* event = nullptr;
    std::vector<RecoTrack>* tracks = nullptr;
    std::vector<RecoCluster>* clusters = nullptr;
    std::vector<RecoVertex>* vertices = nullptr;

    tree->SetBranchAddress("RecoEvent", &event);
    tree->SetBranchAddress("Tracks", &tracks);
    tree->SetBranchAddress("Clusters", &clusters);
    tree->SetBranchAddress("Vertices", &vertices);

    // Create some example histograms
    TH1F* hNTracks = new TH1F("hNTracks", "Number of Tracks per Event;N_{tracks};Events", 20, 0, 20);
    TH1F* hNClusters = new TH1F("hNClusters", "Number of Clusters per Event;N_{clusters};Events", 50, 0, 50);
    TH1F* hTotalEnergy = new TH1F("hTotalEnergy", "Total Event Energy;Energy [GeV];Events", 100, 0, 10);
    TH1F* hTrackMomentum = new TH1F("hTrackMomentum", "Track Momentum;Momentum [GeV];Tracks", 100, 0, 5);
    TH1F* hClusterEnergy = new TH1F("hClusterEnergy", "Cluster Energy;Energy [GeV];Clusters", 100, 0, 5);

    // Loop over events
    Long64_t nEntries = tree->GetEntries();
    std::cout << "\n Processing " << nEntries << " events..." << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Fill event-level histograms
        hNTracks->Fill(event->nTracks);
        hNClusters->Fill(event->nClusters);
        hTotalEnergy->Fill(event->totalEnergy);

        // Fill track histograms
        for (const auto& track : *tracks) {
            hTrackMomentum->Fill(track.momentum);
        }

        // Fill cluster histograms
        for (const auto& cluster : *clusters) {
            hClusterEnergy->Fill(cluster.energy);
        }
    }

    // Create canvas and draw histograms
    TCanvas* c1 = new TCanvas("c1", "Reconstruction Analysis", 1200, 800);
    c1->Divide(3, 2);

    c1->cd(1);
    hNTracks->Draw();

    c1->cd(2);
    hNClusters->Draw();

    c1->cd(3);
    hTotalEnergy->Draw();

    c1->cd(4);
    hTrackMomentum->Draw();

    c1->cd(5);
    hClusterEnergy->Draw();

    c1->Update();

    std::cout << "\n Analysis complete!" << std::endl;
    std::cout << " Total events processed: " << nEntries << std::endl;
    std::cout << "\n==================================================" << std::endl;
}
