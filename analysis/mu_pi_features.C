 /***************************************************************************
 * mu_pi_features.C
 * 
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 * 
 * Created: 26/06/2025
 * 
 * Description:
 *   Extract particle-level features from calorimeter hit collections
 *   relevant for muon/pion separation
 * 
 ***************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void draw_histograms(TH1F* hMuon, TH1F* hPion, const char* title, const char* outName) {

    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", title, 800, 600);
    canvas->SetLogy(); // log-scale for y-axis

    // Draw muon histogram
    hMuon->Scale(1.0 / hMuon->Integral(), "width");
    hMuon->SetLineColor(kBlue);
    hMuon->SetLineWidth(2);
    hMuon->SetFillColor(kBlue-10);
    hMuon->Draw("hist");

    // Draw pion histogram
    hPion->Scale(1.0 / hPion->Integral(), "width");
    hPion->SetLineColor(kRed);
    hPion->SetLineWidth(2);
    hPion->SetFillColor(kRed);
    hPion->SetFillStyle(3004);
    hPion->Draw("hist same");

    // Create legend
    TLegend *legend = new TLegend(0.75, 0.75, 0.90, 0.90);  // x1, y1, x2, y2 in NDC
    legend->SetBorderSize(0);
    legend->SetTextSize(0.050);
    legend->AddEntry(hMuon, "Muon", "f");
    legend->AddEntry(hPion, "Pion", "f");
    legend->Draw();

    canvas->SaveAs(outName);
    delete canvas;

}

void mu_pi_features(const char* inputFileName, const char* outputFileName, const char* inputTreeName = "AnaTree", const char* outputTreeName = "MuonPionTree") {

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

    // Open the analysis ROOT file
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return;
    }
    
    // Get the analysis TTree
    TTree* inputTree = (TTree*)inputFile->Get(inputTreeName);
    if (!inputTree) {
        std::cerr << "Error: Could not find TTree " << inputTreeName << std::endl;
        inputFile->Close();
        return;
    }
    
    // Declare pointers to vectors
    std::vector<Short_t>* trackID = nullptr;
    std::vector<Int_t>*   pdgCode = nullptr;
    std::vector<Float_t>* momX = nullptr;
    std::vector<Float_t>* momY = nullptr;
    std::vector<Float_t>* momZ = nullptr;
    std::vector<Short_t>* ecalHitTrackID = nullptr;
    std::vector<Float_t>* ecalHitEdep = nullptr;
    
    // Set branch addresses
    inputTree->SetBranchAddress("trackID",               &trackID);
    inputTree->SetBranchAddress("pdgCode",               &pdgCode);
    inputTree->SetBranchAddress("momX",                     &momX);
    inputTree->SetBranchAddress("momY",                     &momY);
    inputTree->SetBranchAddress("momZ",                     &momZ);
    inputTree->SetBranchAddress("ecalHitTrackID", &ecalHitTrackID);
    inputTree->SetBranchAddress("ecalHitEdep",       &ecalHitEdep);

    // Create output file
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }

    // Create output TTree
    TTree* outputTree = new TTree(outputTreeName, outputTreeName);

    // Define output variables
    Int_t   outPdgCode;
    Float_t outMomentum;
    Float_t outEcalTotalE;
    Float_t outEcalMeanE;

    // Declare branches for output ntuple
    outputTree->Branch("pdgCode",       &outPdgCode);
    outputTree->Branch("momentum",     &outMomentum);
    outputTree->Branch("EcalTotalE", &outEcalTotalE);
    outputTree->Branch("EcalMeanE",   &outEcalMeanE);
    
    // Create some histograms
    TH1F* hMuonTotalE = new TH1F("hMuonTotalE", "Muon total ECal energy;Total energy [MeV];Counts", 50, 0, 500);
    TH1F* hPionTotalE = new TH1F("hPionTotalE", "Pion total ECal energy;Total energy [MeV];Counts", 50, 0, 500);

    TH1F* hMuonMeanE = new TH1F("hMuonMeanE", "Muon mean ECal energy;Mean energy [MeV];Counts", 50, 0, 1);
    TH1F* hPionMeanE = new TH1F("hPionMeanE", "Pion mean ECal energy;Mean energy [MeV];Counts", 50, 0, 1);
    
    // Get number of entries
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Total number of entries: " << nEntries << std::endl;
    
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {

        // Load the entry
        inputTree->GetEntry(iEntry);
        
        // Get number of particles and ECal deposits in event
        size_t n_particles = pdgCode->size();
        size_t n_ecal_deposits = ecalHitTrackID->size();
        
        for (size_t i = 0; i < n_particles; i++) {

            // Get current trackID and PDG code
            int this_id = trackID->at(i);
            int this_pdg = std::abs(pdgCode->at(i));

            // Skip if particle is not a muon or a charged pion
            if (! ((this_pdg == 13) || (this_pdg == 211))) continue;

            // Compute 3-momentum norm
            float this_momX = momX->at(i);
            float this_momY = momY->at(i);
            float this_momZ = momZ->at(i);
            float this_mom = std::sqrt(this_momX*this_momX + this_momY*this_momY + this_momZ*this_momZ);

            // Extract ECal edeps associated to this particle
            std::vector<float> edeps;
            for (size_t j = 0; j < n_ecal_deposits; j++) {

                if (ecalHitTrackID->at(j) == this_id) {

                    edeps.push_back(ecalHitEdep->at(j));
                
                }

            } // end loop over ECal energy deposits

            // Skip if no ECal edeps associated
            if (edeps.size() == 0) continue;

            // Compute total energy deposited
            float total_e = std::accumulate(edeps.begin(), edeps.end(), 0.0);
            float mean_e = total_e / edeps.size();

            // Add to relevant histogram
            switch (this_pdg) {
                case 13:
                    hMuonTotalE->Fill(total_e);
                    hMuonMeanE->Fill(mean_e);
                case 211:
                    hPionTotalE->Fill(total_e);
                    hPionMeanE->Fill(mean_e);
            }

            // Assign variables for output tree
            outPdgCode = this_pdg;
            outMomentum = this_mom;
            outEcalTotalE = total_e;
            outEcalMeanE = mean_e;

            outputTree->Fill();

        } // end loop over particles

    } // end loop over events
    
    // Draw histograms
    draw_histograms(hMuonTotalE, hPionTotalE, "Total ECal energy", "total_ecal_energy.png");
    draw_histograms(hMuonMeanE, hPionMeanE, "Mean ECal energy", "mean_ecal_energy.png");
    
    // Clean up
    delete hMuonTotalE;
    delete hPionTotalE;
    delete hMuonMeanE;
    delete hPionMeanE;

    // Close input file
    inputFile->Close();
    delete inputFile;

    // Write output tree and close
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
}