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
 *   DBSCAN and writes one output entry per cluster. It also shows how to
 *   reach a reconstruction product: the input may or may not hold one, so it
 *   is asked for by name in BeginJob() and the handle is tested before use.
 *
 *   The input is a simulation file, or a reconstruction file, which by
 *   default carries the simulation trees along with the reconstruction.
 *
 * Usage:
 *   The macro is run by GArAnalysis, which compiles it and hands it the job:
 *
 *     GArAnalysis -a ExampleAnalysis.C -i sim.root -o example_out.root
 *
 *   With the GENIE truth record attached, and with the parameters below set
 *   from a job macro rather than left at their defaults:
 *
 *     GArAnalysis -a ExampleAnalysis.C -i sim.root -o example_out.root \
 *                 -g genie.gst.root -m ExampleAnalysis.mac
 *
 *   From anywhere other than the build directory, source build/setup.sh
 *   first, or give the macros by their full path.
 *
 * Parameters (see ExampleAnalysis.mac):
 *   /ana/eps      cluster radius handed to DBSCAN [cm]
 *   /ana/minPts   neighbours a hit needs to be a core hit
 *
 ***************************************************************************/

// GArAnalysis has libGArAnalysis loaded before it compiles this macro, so
// there is no R__LOAD_LIBRARY here. To compile it by hand in ROOT instead
// (.L ExampleAnalysis.C+), start ROOT from the build directory or with its
// rootlogon.C, which is what loads the library there.

#include <iostream>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector3.h"

#include "AnalysisBase.hh"
#include "Clustering.hh"
#include "DigiDataTypes.hh"
#include "PlotStyle.hh"
#include "SimDataTypes.hh"

/* -------------------------------------------------------------------------- */
/*                                  Analysis                                  */
/* -------------------------------------------------------------------------- */

class ECalClusterAnalysis : public ana::AnalysisBase {
protected:

    /* ------------------------------------------------------------------ */
    /* Called once, before anything is opened, with whatever the job macro */
    /* set. A parameter the macro did not mention leaves its member alone, */
    /* so the value it is declared with below is the default -- there is   */
    /* no second place where the defaults have to be kept in step.         */
    /* ------------------------------------------------------------------ */
    void Configure(const ana::ParameterSet& params) override
    {
        params.Get("eps",    fEps);
        params.Get("minPts", fMinPts);
    }

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
        fHRecoClusterE = new TH1F("hRecoClusterE",
                                  "Reconstructed TPC cluster energy;Energy [MeV];Clusters",
                                  50, 0, 10);

        // Reconstruction products are asked for by name, because which ones a
        // file holds depends on how the reconstruction was configured.
        // Optional() gives back a handle that is simply never valid when the
        // product is absent, so this analysis also runs on plain simulation
        // files; Require() would end the job instead, before the first event.
        fClusters = Optional<std::vector<digi::TPCCluster>>("TPCClusters");
    }

    /* ------------------------------------------------------------------ */
    /* Called once per event. This is the only method an analysis must     */
    /* implement.                                                          */
    /* ------------------------------------------------------------------ */
    void Run() override
    {
        if (!sim.IsValid()) return;

        fEventID = sim.eventID;

        // GENIE truth is only available when a gst file was supplied
        fNuEnergy = HasGenie() ? genie.Ev : -1.;

        /* --------------- Cluster the ECal barrel hits by layer --------- */

        std::unordered_map<Int_t, std::vector<TVector3>> pointsByLayer;
        std::unordered_map<Int_t, std::vector<Float_t>> edepByLayer;

        for (size_t k = 0; k < sim.NECalHits(); ++k) {
            const root::ECalHit& hit = sim.ECalHit(k);
            if (hit.detID != ana::kBarrel) continue;  // barrel only
            pointsByLayer[hit.layer].push_back(sim.ECalHitPosition(k));
            edepByLayer[hit.layer].push_back(hit.energyDeposit);
        } // end loop over ECal hits

        for (const auto& [layer, points] : pointsByLayer) {

            const std::vector<Float_t>& edeps = edepByLayer[layer];
            const std::vector<ana::Cluster> clusters = ana::DBSCAN3D(points, fEps, fMinPts);

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
            if (sim.MotherID(i) != 0) continue;  // primaries only
            const Double_t ecalEdep = sim.ECalEdepOfTrack(sim.TrackID(i));
            if (ecalEdep > 0.) fHPrimaryEcalE->Fill(ecalEdep);
        } // end loop over particles

        /* --------- Reconstructed TPC clusters, when the file has them --- */

        // fClusters was asked for with Optional(), so it is only valid when
        // the input has been through a reconstruction that produced them
        if (fClusters) {
            for (const digi::TPCCluster& cluster : *fClusters) {
                fHRecoClusterE->Fill(cluster.energy);
            }
        }
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
                  << "Reco TPC clusters:  "
                  << (fClusters.IsValid() ? "read from the input"
                                          : "not in the input") << "\n"
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
    // Parameters for DBSCAN, with the values a job macro that says nothing
    // about them gets
    Double_t fEps = 0.5;   // cluster radius [cm]
    Int_t fMinPts = 5;     // minimum neighbours

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

    // Reconstruction products this analysis can make use of
    ana::Handle<std::vector<digi::TPCCluster>> fClusters;

    // Histograms (owned by the output file)
    TH1F* fHClustersPerLayer = nullptr;
    TH1F* fHClusterEdep = nullptr;
    TH1F* fHPrimaryEcalE = nullptr;
    TH1F* fHRecoClusterE = nullptr;

    // Counters
    Long64_t fTotalClusters = 0;
};

/* -------------------------------------------------------------------------- */
/*                       What GArAnalysis runs from here                      */
/* -------------------------------------------------------------------------- */

// The one line every analysis macro ends with. Which files to read, where to
// write and how many events to do are the job's business, not the analysis's,
// and come from the GArAnalysis command line and the job macro.
ANA_ANALYSIS(ECalClusterAnalysis)
