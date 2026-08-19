 /***************************************************************************
 * ExampleAnalysis.C
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Example FastGArSim analysis macro, and the starting point for writing
 *   new ones.
 *
 *   Everything that used to be copy-pasted into each macro -- opening files,
 *   finding trees, wiring branches, the event loop, progress reporting and
 *   writing the output -- lives in the compiled library (ana::AnalysisBase).
 *   All this macro does is derive from it and implement Run(), which is
 *   called once per event with `sim`, `genie` and `geo` already filled.
 *
 *   The analysis itself clusters the ECal barrel hits of each layer with
 *   DBSCAN and writes one output entry per cluster.
 *
 * Usage:
 *   From the build directory (rootlogon.C sets everything up):
 *
 *     root -l 'macros/ExampleAnalysis.C("ntuple.root", "example_out.root")'
 *
 *   With the GENIE truth record attached:
 *
 *     root -l 'macros/ExampleAnalysis.C("ntuple.root", "example_out.root", "genie.gst.root")'
 *
 *   From anywhere else, source build/setup.sh first.
 *
 ***************************************************************************/

R__LOAD_LIBRARY(libGArAnalysis)

#include <iostream>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector3.h"

#include "AnalysisBase.hh"
#include "Clustering.hh"
#include "PlotStyle.hh"

/* -------------------------------------------------------------------------- */
/*                              Analysis parameters                           */
/* -------------------------------------------------------------------------- */

// Parameters for DBSCAN
const Double_t kEps = 0.5;   // cluster radius (in cm)
const Int_t kMinPts = 5;     // minimum neighbours

/* -------------------------------------------------------------------------- */
/*                                  Analysis                                  */
/* -------------------------------------------------------------------------- */

class ECalClusterAnalysis : public ana::AnalysisBase {
protected:

    /* ------------------------------------------------------------------ */
    /* Called once, before the event loop. Book output branches and        */
    /* histograms here -- histograms are attached to the output file       */
    /* automatically and written out at the end of the job.                */
    /* ------------------------------------------------------------------ */
    void BeginJob() override
    {
        // One entry per reconstructed cluster. Output() is null when
        // AnalysisConfig::outputFile is left empty.
        if (TTree* out = Output()) {
            out->Branch("eventID",       &fEventID);
            out->Branch("layer",         &fLayer);
            out->Branch("nHitsInLayer",  &fNHitsInLayer);
            out->Branch("nClusters",     &fNClustersInLayer);
            out->Branch("clusterNHits",  &fClusterNHits);
            out->Branch("clusterEdep",   &fClusterEdep);
            out->Branch("clusterRadius", &fClusterRadius);
            out->Branch("clusterX",      &fClusterX);
            out->Branch("clusterY",      &fClusterY);
            out->Branch("clusterZ",      &fClusterZ);
            out->Branch("nuEnergy",      &fNuEnergy);
        }

        fHClustersPerLayer = new TH1F("hClustersPerLayer",
                                      "Clusters per ECal barrel layer;Clusters;Layers",
                                      20, 0, 20);
        fHClusterEdep = new TH1F("hClusterEdep",
                                 "Cluster energy;Energy deposit [MeV];Clusters",
                                 50, 0, 50);
        fHPrimaryEcalE = new TH1F("hPrimaryEcalE",
                                  "Primary particle ECal energy;Energy deposit [MeV];Particles",
                                  50, 0, 500);
    }

    /* ------------------------------------------------------------------ */
    /* Called once per event. This is the only method an analysis must     */
    /* implement.                                                          */
    /* ------------------------------------------------------------------ */
    void Run() override
    {
        // The hit branches are optional, so check before dereferencing
        if (!sim.ecalHitDetID || !sim.ecalHitLayer || !sim.ecalHitEdep) return;

        fEventID = sim.eventID;

        // GENIE truth is only available when a gst file was supplied
        fNuEnergy = HasGenie() ? genie.Ev : -1.;

        /* --------------- Cluster the ECal barrel hits by layer --------- */

        std::unordered_map<Int_t, std::vector<TVector3>> pointsByLayer;
        std::unordered_map<Int_t, std::vector<Float_t>> edepByLayer;

        for (size_t k = 0; k < sim.NECalHits(); ++k) {
            if (sim.ecalHitDetID->at(k) != ana::kBarrel) continue;  // barrel only
            const Int_t layer = sim.ecalHitLayer->at(k);
            pointsByLayer[layer].push_back(sim.ECalHitPosition(k));
            edepByLayer[layer].push_back(sim.ecalHitEdep->at(k));
        } // end loop over ECal hits

        for (const auto& [layer, points] : pointsByLayer) {

            const std::vector<Float_t>& edeps = edepByLayer[layer];
            const std::vector<ana::Cluster> clusters = ana::DBSCAN3D(points, kEps, kMinPts);

            fLayer = layer;
            fNHitsInLayer = static_cast<Int_t>(points.size());
            fNClustersInLayer = static_cast<Int_t>(clusters.size());
            fHClustersPerLayer->Fill(fNClustersInLayer);

            for (const ana::Cluster& cluster : clusters) {
                const TVector3 centroid = cluster.Centroid(points, edeps);

                fClusterNHits  = static_cast<Int_t>(cluster.Size());
                fClusterEdep   = cluster.TotalWeight(edeps);
                fClusterRadius = cluster.Radius(points);
                fClusterX = centroid.X();
                fClusterY = centroid.Y();
                fClusterZ = centroid.Z();

                fHClusterEdep->Fill(fClusterEdep);
                fTotalClusters++;

                Fill();  // one output entry per cluster
            } // end loop over clusters
        } // end loop over layers

        /* ------------- Energy deposited by the primary particles ------- */

        for (size_t i = 0; i < sim.NParticles(); ++i) {
            if (sim.motherID && sim.motherID->at(i) != 0) continue;  // primaries only
            const Double_t ecalEdep = sim.ECalEdepOfTrack(sim.trackID->at(i));
            if (ecalEdep > 0.) fHPrimaryEcalE->Fill(ecalEdep);
        } // end loop over particles
    }

    /* ------------------------------------------------------------------ */
    /* Called once, after the event loop and before the output is written. */
    /* ------------------------------------------------------------------ */
    void EndJob() override
    {
        std::cout << "\n=== Summary ===\n"
                  << "Events processed:   " << NEvents() << "\n"
                  << "Clusters found:     " << fTotalClusters << "\n"
                  << "ECal barrel layers: " << geo.ECalBarrelLayers() << "\n"
                  << std::endl;

        // Shared plotting style, then save a quick look at the results
        ana::SetPlotStyle();

        TCanvas canvas("cClusterEdep", "Cluster energy", 800, 600);
        canvas.SetLogy();
        fHClusterEdep->SetLineWidth(2);
        fHClusterEdep->Draw("hist");
        canvas.SaveAs("example_cluster_energy.png");
    }

private:
    // Output branch variables
    Int_t fEventID = 0;
    Int_t fLayer = 0;
    Int_t fNHitsInLayer = 0;
    Int_t fNClustersInLayer = 0;
    Int_t fClusterNHits = 0;
    Double_t fClusterEdep = 0.;
    Double_t fClusterRadius = 0.;
    Double_t fClusterX = 0., fClusterY = 0., fClusterZ = 0.;
    Double_t fNuEnergy = -1.;

    // Histograms (owned by the output file)
    TH1F* fHClustersPerLayer = nullptr;
    TH1F* fHClusterEdep = nullptr;
    TH1F* fHPrimaryEcalE = nullptr;

    // Counters
    Long64_t fTotalClusters = 0;
};

/* -------------------------------------------------------------------------- */
/*                                Main function                               */
/* -------------------------------------------------------------------------- */

void ExampleAnalysis(const char* inputFilesG4,
                     const char* outputFileName,
                     const char* inputFilesGENIE = "",
                     Long64_t maxEvents = -1)
{
    ana::AnalysisConfig config;
    config.simFiles = inputFilesG4;          // path, wildcard or comma-separated list
    config.outputFile = outputFileName;
    config.genieFiles = inputFilesGENIE;     // leave empty for gun samples
    config.outputTreeName = "ECalClusters";
    config.maxEvents = maxEvents;            // -1 processes the whole file

    ECalClusterAnalysis analysis;
    analysis.Execute(config);
}
