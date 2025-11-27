 /***************************************************************************
 * EventToNtupleConverter.C
 * 
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 * 
 * Created: June 2025
 * 
 * Description:
 *   Convert ROOT file with Event objects from FastGArSim into flat ntuple
 * 
 ***************************************************************************/

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
// #include "../include/ROOTDataTypes.hh"
#include <iostream>
#include <vector>

using namespace root;

// Converter class
class EventConverter {
private:
    TFile* inputFile;
    TFile* outputFile;
    TTree* inputTree;
    TTree* outputTree;
    TTree* geometryTree;
    Event* event;
    
    // Variables for analysis ntuple
    Int_t    eventID;

    // Particle properties
    std::vector<Int_t> trackID;
    std::vector<Int_t> pdgCode;
    std::vector<Int_t> motherID;
    std::vector<std::string> creatorProcess;
    std::vector<std::string> endProcess;

    // Start and end trajectory points
    std::vector<Float_t> startX, startY, startZ;
    std::vector<Float_t> endX, endY, endZ;

    // Initial and final momenta
    std::vector<Float_t> startPX, startPY, startPZ;
    std::vector<Float_t> endPX, endPY, endPZ;

    // TPC hits
    std::vector<Int_t> tpcHitTrackID;  // Which particle this hit belongs to
    std::vector<Bool_t> tpcHitIsSec;
    std::vector<Float_t> tpcHitX, tpcHitY, tpcHitZ;
    std::vector<Float_t> tpcHitEdep;
    std::vector<Float_t> tpcHitStepSize;

    // ECal hits (both primary and secondary)
    std::vector<Int_t> ecalHitTrackID;
    std::vector<Bool_t> ecalHitIsSec;
    std::vector<Float_t> ecalHitX, ecalHitY, ecalHitZ;
    std::vector<Float_t> ecalHitTime;
    std::vector<Float_t> ecalHitEdep;
    std::vector<Int_t> ecalHitSegment;
    std::vector<Int_t> ecalHitLayer;
    std::vector<Int_t> ecalHitDetID;

    // MuID hits (both primary and secondary)
    std::vector<Int_t> muidHitTrackID;
    std::vector<Bool_t> muidHitIsSec;
    std::vector<Float_t> muidHitX, muidHitY, muidHitZ;
    std::vector<Float_t> muidHitTime;
    std::vector<Float_t> muidHitEdep;
    std::vector<Int_t> muidHitSegment;
    std::vector<Int_t> muidHitLayer;
    std::vector<Int_t> muidHitDetID;
    
public:
    EventConverter() : inputFile(nullptr), outputFile(nullptr), 
                               inputTree(nullptr), outputTree(nullptr),
                               geometryTree(nullptr),
                               event(nullptr) {}
    
    ~EventConverter() {
        delete event;
        if (inputFile && inputFile->IsOpen()) inputFile->Close();
        if (outputFile && outputFile->IsOpen()) outputFile->Close();
        delete inputFile;
        delete outputFile;
    }
    
    Bool_t Initialize(const char* inputFileName, const char* outputFileName, 
                      const char* treeName = "Events") {

        // Open input file
        inputFile = TFile::Open(inputFileName, "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Error: Cannot open input file " << inputFileName << std::endl;
            return kFALSE;
        }
        
        // Get input tree
        inputTree = dynamic_cast<TTree*>(inputFile->Get(treeName));
        if (!inputTree) {
            std::cerr << "Error: Cannot find tree " << treeName << " in input file" << std::endl;
            return kFALSE;
        }

        // Try to get Geometry tree (optional - won't fail if it doesn't exist)
        geometryTree = dynamic_cast<TTree*>(inputFile->Get("Geometry"));
        if (geometryTree) {
            std::cout << "Found Geometry tree!" << std::endl;
        } else {
            std::cout << "Warning: No Geometry tree found in input file" << std::endl;
        }
        
        // Create Event object and set branch address
        inputTree->SetBranchAddress("Event", &event);
        
        // Create output file
        outputFile = TFile::Open(outputFileName, "RECREATE");
        if (!outputFile || outputFile->IsZombie()) {
            std::cerr << "Error: Cannot create output file " << outputFileName << std::endl;
            return kFALSE;
        }

        // Set compression level
        outputFile->SetCompressionLevel(5);
        
        // Create output tree with analysis structure
        outputTree = new TTree("AnaTree", "Analysis Tree");
        
        // Define branches for ana ntuple
        outputTree->Branch("eventID", &eventID);
        outputTree->Branch("trackID", &trackID);
        outputTree->Branch("pdgCode", &pdgCode);
        outputTree->Branch("motherID", &motherID);
        outputTree->Branch("creatorProcess", &creatorProcess);
        outputTree->Branch("endProcess", &endProcess);
        
        // Trajectory endpoints
        outputTree->Branch("startX", &startX);
        outputTree->Branch("startY", &startY);
        outputTree->Branch("startZ", &startZ);
        outputTree->Branch("endX",   &endX);
        outputTree->Branch("endY",   &endY);
        outputTree->Branch("endZ",   &endZ);
        
        // Initial and final momenta
        outputTree->Branch("startPX", &startPX);
        outputTree->Branch("startPY", &startPY);
        outputTree->Branch("startPZ", &startPZ);
        outputTree->Branch("endPX", &endPX);
        outputTree->Branch("endPY", &endPY);
        outputTree->Branch("endPZ", &endPZ);

        // TPC hits
        outputTree->Branch("tpcHitTrackID", &tpcHitTrackID);
        outputTree->Branch("tpcHitIsSec", &tpcHitIsSec);
        outputTree->Branch("tpcHitX", &tpcHitX);
        outputTree->Branch("tpcHitY", &tpcHitY);
        outputTree->Branch("tpcHitZ", &tpcHitZ);
        outputTree->Branch("tpcHitEdep", &tpcHitEdep);
        outputTree->Branch("tpcHitStepSize", &tpcHitStepSize);
        
        // ECal hits
        outputTree->Branch("ecalHitTrackID", &ecalHitTrackID);
        outputTree->Branch("ecalHitIsSec", &ecalHitIsSec);
        outputTree->Branch("ecalHitX", &ecalHitX);
        outputTree->Branch("ecalHitY", &ecalHitY);
        outputTree->Branch("ecalHitZ", &ecalHitZ);
        outputTree->Branch("ecalHitTime", &ecalHitTime);
        outputTree->Branch("ecalHitEdep", &ecalHitEdep);
        outputTree->Branch("ecalHitSegment", &ecalHitSegment);
        outputTree->Branch("ecalHitLayer", &ecalHitLayer);
        outputTree->Branch("ecalHitDetID", &ecalHitDetID);

        // MuID hits
        outputTree->Branch("muidHitTrackID", &muidHitTrackID);
        outputTree->Branch("muidHitIsSec", &muidHitIsSec);
        outputTree->Branch("muidHitX", &muidHitX);
        outputTree->Branch("muidHitY", &muidHitY);
        outputTree->Branch("muidHitZ", &muidHitZ);
        outputTree->Branch("muidHitTime", &muidHitTime);
        outputTree->Branch("muidHitEdep", &muidHitEdep);
        outputTree->Branch("muidHitSegment", &muidHitSegment);
        outputTree->Branch("muidHitLayer", &muidHitLayer);
        outputTree->Branch("muidHitDetID", &muidHitDetID);
        
        return kTRUE;
    }

    void ClearVectors() {

        trackID.clear();
        pdgCode.clear();
        motherID.clear();
        creatorProcess.clear();
        endProcess.clear();
        
        startX.clear(); startY.clear(); startZ.clear();
        endX.clear(); endY.clear(); endZ.clear();

        startPX.clear(); startPY.clear(); startPZ.clear();
        endPX.clear(); endPY.clear(); endPZ.clear();
        
        tpcHitTrackID.clear();
        tpcHitIsSec.clear();
        tpcHitX.clear(); tpcHitY.clear(); tpcHitZ.clear();
        tpcHitEdep.clear(); tpcHitStepSize.clear();
        
        ecalHitTrackID.clear();
        ecalHitIsSec.clear();
        ecalHitX.clear(); ecalHitY.clear(); ecalHitZ.clear();
        ecalHitTime.clear();
        ecalHitEdep.clear();
        ecalHitSegment.clear();
        ecalHitLayer.clear();
        ecalHitDetID.clear();

        muidHitTrackID.clear();
        muidHitIsSec.clear();
        muidHitX.clear(); muidHitY.clear(); muidHitZ.clear();
        muidHitTime.clear();
        muidHitEdep.clear();
        muidHitSegment.clear();
        muidHitLayer.clear();
        muidHitDetID.clear();

    }
    
    void ProcessEvents() {
        if (!inputTree || !outputTree) {
            std::cerr << "Error: Trees not properly initialized" << std::endl;
            return;
        }
        
        Long64_t nEntries = inputTree->GetEntries();
        std::cout << "Processing " << nEntries << " events..." << std::endl;
        
        // Progress counter
        Int_t reportEvery = nEntries / 10;
        if (reportEvery == 0) reportEvery = 1;
 
        // Loop over all events
        for (Long64_t i = 0; i < nEntries; i++) {
            if (i % reportEvery == 0) {
                std::cout << "Processing event " << i << " (" 
                          << (100.0 * i / nEntries) << "%)" << std::endl;
            }
            
            // Get entry
            inputTree->GetEntry(i);

            // Clear vectors
            ClearVectors();
            
            // Extract information from Event object
            eventID = event->eventID;

            // Loop over particles in the event
            const auto& particles = event->particles;
            for (size_t particleIdx = 0; particleIdx < particles.size(); ++particleIdx) {
                const auto& particle = particles[particleIdx];
                
                // Store particle properties
                trackID.push_back(particle.trackID);
                pdgCode.push_back(particle.pdgCode);
                motherID.push_back(particle.motherID);
                creatorProcess.push_back(particle.creatorProcess.Data());
                endProcess.push_back(particle.endProcess.Data());
                
                // Get trajectory points
                const auto& trajPoints = particle.trajectory.points;
                if (!trajPoints.empty()) {
                    // First trajectory point
                    const auto& firstPoint = trajPoints.front();
                    startX.push_back(firstPoint.x);
                    startY.push_back(firstPoint.y);
                    startZ.push_back(firstPoint.z);
                    // Last trajectory point
                    const auto& lastPoint = trajPoints.back();
                    endX.push_back(lastPoint.x);
                    endY.push_back(lastPoint.y);
                    endZ.push_back(lastPoint.z);
                } else {
                    // No trajectory points - fill with bogus
                    startX.push_back(-9999.0);
                    startY.push_back(-9999.0);
                    startZ.push_back(-9999.0);
                    endX.push_back(-9999.0);
                    endY.push_back(-9999.0);
                    endZ.push_back(-9999.0);
                }

                const auto& momPoints = particle.trajectory.mom_points;
                if (!momPoints.empty()) {
                    // First trajectory point
                    const auto& firstPoint = momPoints.front();
                    startPX.push_back(firstPoint.x);
                    startPY.push_back(firstPoint.y);
                    startPZ.push_back(firstPoint.z);
                    // Last trajectory point
                    const auto& lastPoint = momPoints.back();
                    endPX.push_back(lastPoint.x);
                    endPY.push_back(lastPoint.y);
                    endPZ.push_back(lastPoint.z);
                } else {
                    // No momentum points - fill with bogus
                    startPX.push_back(-9999.0);
                    startPY.push_back(-9999.0);
                    startPZ.push_back(-9999.0);
                    endPX.push_back(-9999.0);
                    endPY.push_back(-9999.0);
                    endPZ.push_back(-9999.0);
                }
                
                // Process TPC hits
                const auto& tpcHits = particle.tpcHits;
                for (const auto& hit : tpcHits) {
                    tpcHitTrackID.push_back(particle.trackID);
                    tpcHitIsSec.push_back(false);
                    tpcHitX.push_back(hit.x);
                    tpcHitY.push_back(hit.y);
                    tpcHitZ.push_back(hit.z);
                    tpcHitEdep.push_back(hit.energyDeposit);
                    tpcHitStepSize.push_back(hit.stepSize);
                }

                // Process secondary TPC hits
                const auto& sec_tpcHits = particle.sec_tpcHits;
                for (const auto& hit : sec_tpcHits) {
                    tpcHitTrackID.push_back(particle.trackID);
                    tpcHitIsSec.push_back(true);
                    tpcHitX.push_back(hit.x);
                    tpcHitY.push_back(hit.y);
                    tpcHitZ.push_back(hit.z);
                    tpcHitEdep.push_back(hit.energyDeposit);
                    tpcHitStepSize.push_back(hit.stepSize);
                }
                
                // Process ECal hits
                const auto& ecalHits = particle.ecalHits;
                for (const auto& hit : ecalHits) {
                    ecalHitTrackID.push_back(particle.trackID);
                    ecalHitIsSec.push_back(false);
                    ecalHitX.push_back(hit.x);
                    ecalHitY.push_back(hit.y);
                    ecalHitZ.push_back(hit.z);
                    ecalHitTime.push_back(hit.time);
                    ecalHitEdep.push_back(hit.energyDeposit);
                    ecalHitSegment.push_back(hit.segment);
                    ecalHitLayer.push_back(hit.layer);
                    ecalHitDetID.push_back(hit.detID);
                }

                // Process secondary ECal hits
                const auto& sec_ecalHits = particle.sec_ecalHits;
                for (const auto& hit : sec_ecalHits) {
                    ecalHitTrackID.push_back(particle.trackID);
                    ecalHitIsSec.push_back(true);
                    ecalHitX.push_back(hit.x);
                    ecalHitY.push_back(hit.y);
                    ecalHitZ.push_back(hit.z);
                    ecalHitTime.push_back(hit.time);
                    ecalHitEdep.push_back(hit.energyDeposit);
                    ecalHitSegment.push_back(hit.segment);
                    ecalHitLayer.push_back(hit.layer);
                    ecalHitDetID.push_back(hit.detID);
                }

                // Process MuID hits
                const auto& muidHits = particle.muidHits;
                for (const auto& hit : muidHits) {
                    muidHitTrackID.push_back(particle.trackID);
                    muidHitIsSec.push_back(false);
                    muidHitX.push_back(hit.x);
                    muidHitY.push_back(hit.y);
                    muidHitZ.push_back(hit.z);
                    muidHitTime.push_back(hit.time);
                    muidHitEdep.push_back(hit.energyDeposit);
                    muidHitSegment.push_back(hit.segment);
                    muidHitLayer.push_back(hit.layer);
                    muidHitDetID.push_back(hit.detID);
                }

                // Process secondary MuID hits
                const auto& sec_muidHits = particle.sec_muidHits;
                for (const auto& hit : sec_muidHits) {
                    muidHitTrackID.push_back(particle.trackID);
                    muidHitIsSec.push_back(true);
                    muidHitX.push_back(hit.x);
                    muidHitY.push_back(hit.y);
                    muidHitZ.push_back(hit.z);
                    muidHitTime.push_back(hit.time);
                    muidHitEdep.push_back(hit.energyDeposit);
                    muidHitSegment.push_back(hit.segment);
                    muidHitLayer.push_back(hit.layer);
                    muidHitDetID.push_back(hit.detID);
                }

            }
            
            // Fill output tree
            outputTree->Fill();
        }
        
        std::cout << "Processing complete!" << std::endl;
    }
    
    void Finalize() {
        if (outputFile && outputFile->IsOpen()) {
            outputFile->cd();

            // Write the analysis tree
            outputTree->OptimizeBaskets();
            outputTree->Write();

            // Clone the geometry tree if it exists
            if (geometryTree) {
                std::cout << "Copying Geometry tree to output file..." << std::endl;
                TTree* geometryClone = geometryTree->CloneTree();
                geometryClone->SetName("GeoTree");
                geometryClone->Write();
                std::cout << "Geometry tree copied successfully" << std::endl;
            }
            
            // Print summary
            std::cout << "\n=== Summary ===" << std::endl;
            std::cout << "Input file:      " << inputFile->GetName() << std::endl;
            std::cout << "Output file:     " << outputFile->GetName() << std::endl;
            std::cout << "File size:       " << outputFile->GetSize()/1024.0/1024.0 << " MB" << std::endl;
            std::cout << "Output events:   " << outputTree->GetEntries() << std::endl;
            
            outputFile->Close();
        }
    }
    
    void Run(const char* inputFileName, const char* outputFileName, 
             const char* treeName = "Events") {
        if (Initialize(inputFileName, outputFileName, treeName)) {
            ProcessEvents();
            Finalize();
        }
    }
};

// Main macro function
void EventToNtupleConverter(const char* inputFile, const char* outputFile, 
                           const char* treeName = "Events") {

    // Load the dictionary if not already loaded
    gSystem->Load("libROOTDataDict");

    EventConverter converter;
    converter.Run(inputFile, outputFile, treeName);

}
