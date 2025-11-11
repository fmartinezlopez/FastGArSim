 /***************************************************************************
 * neutron_analysis.C
 * 
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 * 
 * Created: 11/11/2025
 * 
 * Description:
 *   Neutron identification analysis
 *   Uses DBSCAN algorithm to cluster ECal hits
 * 
 ***************************************************************************/

#include <iostream>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <cmath>
#include <tuple>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

/* -------------------------------------------------------------------------- */
/*                              Useful constants                              */
/* -------------------------------------------------------------------------- */

// Fixed max dimension to read arrays from gst TTree
const Int_t kNPmax = 1000;

// Parameters for DBSCAN
const Double_t eps = 0.5;  // cluster radius (in cm)
const Int_t minPts = 5;    // minimum neighbors

/* -------------------------------------------------------------------------- */
/*                            Clustering algorithm                            */
/* -------------------------------------------------------------------------- */

struct Cluster {
    std::vector<Int_t> indices;
};

struct GridKey {
    Int_t x, y, z;
    Bool_t operator==(const GridKey &other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct KeyHash {
    size_t operator()(const GridKey &k) const {
        return ((51 + std::hash<Int_t>()(k.x)) * 51 + std::hash<Int_t>()(k.y)) * 51 + std::hash<Int_t>()(k.z);
    }
};

// Assign points to a 3D grid for fast neighbor lookup
std::unordered_map<GridKey, std::vector<Int_t>, KeyHash> buildGrid(const std::vector<TVector3> &points, Float_t eps) {
    std::unordered_map<GridKey, vector<Int_t>, KeyHash> grid;
    for (size_t i = 0; i < points.size(); ++i) {
        GridKey key{
            Int_t(floor(points[i].X() / eps)),
            Int_t(floor(points[i].Y() / eps)),
            Int_t(floor(points[i].Z() / eps))
        };
        grid[key].push_back(i);
    }
    return grid;
}

// Retrieve neighbors using the grid
std::vector<Int_t> regionQuery(const std::vector<TVector3> &points, const std::unordered_map<GridKey, vector<Int_t>, KeyHash> &grid,
                               Int_t idx, Double_t eps) {
    std::vector<Int_t> neighbors;
    const auto &p = points[idx];
    Int_t gx = Int_t(floor(p.X() / eps));
    Int_t gy = Int_t(floor(p.Y() / eps));
    Int_t gz = Int_t(floor(p.Z() / eps));

    for (Int_t dx = -1; dx <= 1; ++dx)
        for (Int_t dy = -1; dy <= 1; ++dy)
            for (Int_t dz = -1; dz <= 1; ++dz) {
                GridKey neighborKey{gx + dx, gy + dy, gz + dz};
                auto it = grid.find(neighborKey);
                if (it == grid.end()) continue;
                for (Int_t j : it->second) {
                    if ((points[j] - p).Mag() <= eps)
                        neighbors.push_back(j);
                }
            }
    return neighbors;
}

// DBSCAN-like clustering
std::vector<Cluster> dbscan3D(const std::vector<TVector3> &points, Double_t eps, Int_t minPts) {
    Int_t n = points.size();
    std::vector<Int_t> labels(n, -1); // -1 = unvisited, -2 = noise
    Int_t clusterID = 0;

    auto grid = buildGrid(points, eps);
    std::vector<Cluster> clusters;

    for (Int_t i = 0; i < n; ++i) {
        if (labels[i] != -1) continue;

        auto neighbors = regionQuery(points, grid, i, eps);
        if ((Int_t)neighbors.size() < minPts) {
            labels[i] = -2; // noise
            continue;
        }

        Cluster cluster;
        cluster.indices.push_back(i);
        labels[i] = clusterID;

        std::vector<Int_t> seeds = neighbors;
        for (size_t j = 0; j < seeds.size(); ++j) {
            Int_t idx = seeds[j];
            if (labels[idx] == -2)
                labels[idx] = clusterID; // previously noise, now part of cluster
            if (labels[idx] != -1)
                continue;

            labels[idx] = clusterID;
            cluster.indices.push_back(idx);

            auto subNeighbors = regionQuery(points, grid, idx, eps);
            if ((Int_t)subNeighbors.size() >= minPts)
                seeds.insert(seeds.end(), subNeighbors.begin(), subNeighbors.end());
        }

        clusters.push_back(cluster);
        clusterID++;
    }

    return clusters;
}

/* -------------------------------------------------------------------------- */
/*                               Analyzer class                               */
/* -------------------------------------------------------------------------- */

class Analyzer {
private:
    // Declare TFiles and TTrees
    TFile* inputFileGENIE;
    TFile* inputFileG4;
    TFile* outputFile;
    TTree* inputTreeGENIE;
    TTree* inputTreeG4;
    TTree* outputTree;

    // GST variables
    Int_t iev;             // event number
    Int_t neu;             // neutrino PDG code
    Bool_t qel;            // is QE event?
    Bool_t mec;            // is MEC event?
    Bool_t res;            // is RES event?
    Bool_t dis;            // is DIS event?
    Bool_t coh;            // is COH event?
    Int_t resid;           // baryon resonance ID
    Bool_t cc;             // is CC?
    Bool_t nc;             // is NC?
    Double_t Ev;           // neutrino energy
    Double_t pxv;          // neutrino px
    Double_t pyv;          // neutrino py
    Double_t pzv;          // neutrino pz
    Double_t Q2;           // 4-momentum transfer
    Double_t W;            // hadronic invariant mass
    Double_t x;            // Bjorken x
    Double_t y;            // inelasticity
    Double_t El;           // primary lepton energy
    Double_t pxl;          // primary lepton px
    Double_t pyl;          // primary lepton py
    Double_t pzl;          // primary lepton pz
    Int_t nf;              // number final state hadrons
    Int_t pdgf[kNPmax];    // i-th hadron PDG code
    Double_t Ef[kNPmax];   // i-th hadron energy
    Double_t pxf[kNPmax];  // i-th hadron px
    Double_t pyf[kNPmax];  // i-th hadron py
    Double_t pzf[kNPmax];  // i-th hadron pz

    // FastGArSim variables
    Int_t eventID;
    std::vector<Int_t> *trackID = nullptr;
    std::vector<Int_t> *pdgCode = nullptr;
    std::vector<Int_t> *motherID = nullptr;
    std::vector<std::string> *creatorProcess = nullptr;
    std::vector<std::string> *endProcess = nullptr;
    std::vector<Float_t> *startX = nullptr, *startY = nullptr, *startZ = nullptr;
    std::vector<Float_t> *endX = nullptr, *endY = nullptr, *endZ = nullptr;
    std::vector<Float_t> *startPX = nullptr, *startPY = nullptr, *startPZ = nullptr;
    std::vector<Float_t> *endPX = nullptr, *endPY = nullptr, *endPZ = nullptr;
    std::vector<Int_t> *tpcHitTrackID = nullptr;
    std::vector<Bool_t> *tpcHitIsSec = nullptr;
    std::vector<Float_t> *tpcHitX = nullptr, *tpcHitY = nullptr, *tpcHitZ = nullptr;
    std::vector<Float_t> *tpcHitEdep = nullptr;
    std::vector<Float_t> *tpcHitStepSize = nullptr;
    std::vector<Int_t> *ecalHitTrackID = nullptr;
    std::vector<Bool_t> *ecalHitIsSec = nullptr;
    std::vector<Float_t> *ecalHitX = nullptr, *ecalHitY = nullptr, *ecalHitZ = nullptr;
    std::vector<Float_t> *ecalHitTime = nullptr;
    std::vector<Float_t> *ecalHitEdep = nullptr;
    std::vector<Int_t> *ecalHitLayer = nullptr;
    std::vector<Int_t> *ecalHitDetID = nullptr;
    std::vector<Int_t> *muidHitTrackID = nullptr;
    std::vector<Bool_t> *muidHitIsSec = nullptr;
    std::vector<Float_t> *muidHitX = nullptr, *muidHitY = nullptr, *muidHitZ = nullptr;
    std::vector<Float_t> *muidHitTime = nullptr;
    std::vector<Float_t> *muidHitEdep = nullptr;
    std::vector<Int_t> *muidHitLayer = nullptr;
    std::vector<Int_t> *muidHitDetID = nullptr;

public:
    Analyzer() : inputFileGENIE(nullptr), inputFileG4(nullptr), outputFile(nullptr),
                 inputTreeGENIE(nullptr), inputTreeG4(nullptr), outputTree(nullptr)
    {}

    ~Analyzer()
    {
        if (inputFileGENIE && inputFileGENIE->IsOpen()) inputFileGENIE->Close();
        if (inputFileG4 && inputFileG4->IsOpen()) inputFileG4->Close();
        if (outputFile && outputFile->IsOpen()) outputFile->Close();
        delete inputFileGENIE;
        delete inputFileG4;
        delete outputFile;
    }

    Bool_t Initialize(const char* inputFileNameGENIE, const char* inputFileNameG4, const char* outputFileName,
                      const char* inputTreeNameGENIE, const char* inputTreeNameG4, const char* outputTreeName)
    {

        // Open the GENIE file
        inputFileGENIE = TFile::Open(inputFileNameGENIE, "READ");
        if (!inputFileGENIE || inputFileGENIE->IsZombie()) {
            std::cerr << "Error: Could not open file " << inputFileNameGENIE << std::endl;
            return kFALSE;
        }

        // Get the GENIE GST TTree
        inputTreeGENIE = (TTree*)inputFileGENIE->Get(inputTreeNameGENIE);
        if (!inputTreeGENIE) {
            std::cerr << "Error: Could not find TTree " << inputTreeNameGENIE << std::endl;
            inputFileGENIE->Close();
            return kFALSE;
        }

        // Open the FastGArSim ntuple file
        inputFileG4 = TFile::Open(inputFileNameG4, "READ");
        if (!inputFileG4 || inputFileG4->IsZombie()) {
            std::cerr << "Error: Could not open file " << inputFileNameG4 << std::endl;
            return kFALSE;
        }

        // Get the FastGArSim analysis TTree
        inputTreeG4 = (TTree*)inputFileG4->Get(inputTreeNameG4);
        if (!inputTreeG4) {
            std::cerr << "Error: Could not find TTree " << inputTreeNameG4 << std::endl;
            inputFileG4->Close();
            return kFALSE;
        }
        
        // Set branch addresses for GST TTree
        inputTreeGENIE->SetBranchAddress("iev",     &iev);
        inputTreeGENIE->SetBranchAddress("neu",     &neu);
        inputTreeGENIE->SetBranchAddress("qel",     &qel);
        inputTreeGENIE->SetBranchAddress("mec",     &mec);
        inputTreeGENIE->SetBranchAddress("res",     &res);
        inputTreeGENIE->SetBranchAddress("dis",     &dis);
        inputTreeGENIE->SetBranchAddress("coh",     &coh);
        inputTreeGENIE->SetBranchAddress("resid", &resid);
        inputTreeGENIE->SetBranchAddress("cc",       &cc);
        inputTreeGENIE->SetBranchAddress("nc",       &nc);
        inputTreeGENIE->SetBranchAddress("Ev",       &Ev);
        inputTreeGENIE->SetBranchAddress("pxv",     &pxv);
        inputTreeGENIE->SetBranchAddress("pyv",     &pyv);
        inputTreeGENIE->SetBranchAddress("pzv",     &pzv);
        inputTreeGENIE->SetBranchAddress("Q2",       &Q2);
        inputTreeGENIE->SetBranchAddress("W",         &W);
        inputTreeGENIE->SetBranchAddress("x",         &x);
        inputTreeGENIE->SetBranchAddress("y",         &y);
        inputTreeGENIE->SetBranchAddress("El",       &El);
        inputTreeGENIE->SetBranchAddress("pxl",     &pxl);
        inputTreeGENIE->SetBranchAddress("pyl",     &pyl);
        inputTreeGENIE->SetBranchAddress("pzl",     &pzl);
        inputTreeGENIE->SetBranchAddress("nf",       &nf);
        inputTreeGENIE->SetBranchAddress("pdgf",   &pdgf);
        inputTreeGENIE->SetBranchAddress("Ef",       &Ef);
        inputTreeGENIE->SetBranchAddress("pxf",     &pxf);
        inputTreeGENIE->SetBranchAddress("pyf",     &pyf);
        inputTreeGENIE->SetBranchAddress("pzf",     &pzf);

        // Set branch addresses for FastGArSim TTree
        inputTreeG4->SetBranchAddress("eventID", &eventID);
        inputTreeG4->SetBranchAddress("trackID", &trackID);
        inputTreeG4->SetBranchAddress("pdgCode", &pdgCode);
        inputTreeG4->SetBranchAddress("motherID", &motherID);
        inputTreeG4->SetBranchAddress("creatorProcess", &creatorProcess);
        inputTreeG4->SetBranchAddress("endProcess", &endProcess);
        inputTreeG4->SetBranchAddress("startX", &startX);
        inputTreeG4->SetBranchAddress("startY", &startY);
        inputTreeG4->SetBranchAddress("startZ", &startZ);
        inputTreeG4->SetBranchAddress("endX",   &endX);
        inputTreeG4->SetBranchAddress("endY",   &endY);
        inputTreeG4->SetBranchAddress("endZ",   &endZ);
        inputTreeG4->SetBranchAddress("startPX", &startPX);
        inputTreeG4->SetBranchAddress("startPY", &startPY);
        inputTreeG4->SetBranchAddress("startPZ", &startPZ);
        inputTreeG4->SetBranchAddress("endPX", &endPX);
        inputTreeG4->SetBranchAddress("endPY", &endPY);
        inputTreeG4->SetBranchAddress("endPZ", &endPZ);
        inputTreeG4->SetBranchAddress("tpcHitTrackID", &tpcHitTrackID);
        inputTreeG4->SetBranchAddress("tpcHitIsSec", &tpcHitIsSec);
        inputTreeG4->SetBranchAddress("tpcHitX", &tpcHitX);
        inputTreeG4->SetBranchAddress("tpcHitY", &tpcHitY);
        inputTreeG4->SetBranchAddress("tpcHitZ", &tpcHitZ);
        inputTreeG4->SetBranchAddress("tpcHitEdep", &tpcHitEdep);
        inputTreeG4->SetBranchAddress("tpcHitStepSize", &tpcHitStepSize);
        inputTreeG4->SetBranchAddress("ecalHitTrackID", &ecalHitTrackID);
        inputTreeG4->SetBranchAddress("ecalHitIsSec", &ecalHitIsSec);
        inputTreeG4->SetBranchAddress("ecalHitX", &ecalHitX);
        inputTreeG4->SetBranchAddress("ecalHitY", &ecalHitY);
        inputTreeG4->SetBranchAddress("ecalHitZ", &ecalHitZ);
        inputTreeG4->SetBranchAddress("ecalHitTime", &ecalHitTime);
        inputTreeG4->SetBranchAddress("ecalHitEdep", &ecalHitEdep);
        inputTreeG4->SetBranchAddress("ecalHitLayer", &ecalHitLayer);
        inputTreeG4->SetBranchAddress("ecalHitDetID", &ecalHitDetID);
        inputTreeG4->SetBranchAddress("muidHitTrackID", &muidHitTrackID);
        inputTreeG4->SetBranchAddress("muidHitIsSec", &muidHitIsSec);
        inputTreeG4->SetBranchAddress("muidHitX", &muidHitX);
        inputTreeG4->SetBranchAddress("muidHitY", &muidHitY);
        inputTreeG4->SetBranchAddress("muidHitZ", &muidHitZ);
        inputTreeG4->SetBranchAddress("muidHitTime", &muidHitTime);
        inputTreeG4->SetBranchAddress("muidHitEdep", &muidHitEdep);
        inputTreeG4->SetBranchAddress("muidHitLayer", &muidHitLayer);
        inputTreeG4->SetBranchAddress("muidHitDetID", &muidHitDetID);

        return kTRUE;

    }

    void ProcessEvents() {

        // Get number of events
        Long64_t nEvents = inputTreeG4->GetEntries();
        std::cout << "Number of events: " << nEvents << std::endl;

        // Progress counter
        Int_t reportEvery = nEvents / 10;
        if (reportEvery == 0) reportEvery = 1;

        // Main event loop
        for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {

            if (iEvent > 0) break;

            // Print progress
            if (iEvent % reportEvery == 0) {
                std::cout << "Processing event " << iEvent << " (" 
                        << (100.0 * iEvent / nEvents) << "%)" << std::endl;
            }

            // Load current event
            inputTreeG4->GetEntry(iEvent);

            // Load corresponding entry in GENIE tree
            inputTreeGENIE->GetEntry(eventID);

            // Group all ECal hits in the different layers
            std::unordered_map<Int_t, std::vector<TVector3>> ecalHitLayerMap;
            Long64_t nECalHits = ecalHitTrackID->size();
            for (size_t k = 0; k < nECalHits; ++k) {
                if (ecalHitDetID->at(k) != 2) continue; // only for hits in the barrel
                Int_t layer = ecalHitLayer->at(k);
                ecalHitLayerMap[layer].emplace_back(ecalHitX->at(k), ecalHitY->at(k), ecalHitZ->at(k));
            } // end loop over ECal hits

            // Form clusters in every layer
            for (const auto& [layer, points] : ecalHitLayerMap) {
                auto clusters = dbscan3D(points, eps, minPts);
                std::cout << "Layer:              " << layer << "\n"
                          << "Number of hits:     " << points.size() << "\n"
                          << "Number of clusters: " << clusters.size() << "\n"
                          << std::endl;
            } // end loop over ECal hit groups

        } // end loop over events
    }

    void Finalize() {
        if (outputFile && outputFile->IsOpen()) {
            
            outputFile->Close();
        }
    }

    void Run(const char* inputFileNameGENIE, const char* inputFileNameG4, const char* outputFileName,
             const char* inputTreeNameGENIE, const char* inputTreeNameG4, const char* outputTreeName) {
        if (Initialize(inputFileNameGENIE, inputFileNameG4, outputFileName, inputTreeNameGENIE, inputTreeNameG4, outputTreeName)) {
            ProcessEvents();
            Finalize();
        }
    }

};

/* -------------------------------------------------------------------------- */
/*                                Main function                               */
/* -------------------------------------------------------------------------- */

void neutron_analysis(const char* inputFileNameGENIE, const char* inputFileNameG4, const char* outputFileName,
                      const char* inputTreeNameGENIE = "gst",                  
                      const char* inputTreeNameG4 = "AnaTree",
                      const char* outputTreeName = "NeutronAnalysis") {

    std::cout << "Running neutron_analysis.C\n"
              << "    inputFileNameGENIE: " << inputFileNameGENIE << "\n"
              << "    inputFileNameG4:    " << inputFileNameG4 << "\n"
              << "    outputFileName:     " << outputFileName << "\n"
              << std::endl;

    Analyzer ana;
    ana.Run(inputFileNameGENIE, inputFileNameG4, outputFileName,
            inputTreeNameGENIE, inputTreeNameG4, outputTreeName);

}