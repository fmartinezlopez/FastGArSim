 /***************************************************************************
 * TrackingEfficiency.C
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Tracking efficiency of primary particles in the gas TPC, for
 *   particle-gun samples.
 *
 *   The question this answers is not "did a fitter succeed" -- there is no
 *   pattern recognition in the chain yet -- but the one before it: does the
 *   reconstruction leave behind enough of a primary for a track to be found
 *   at all? Every primary that ionised the gas is followed through to the
 *   reconstruction output and asked whether what survives is track-like.
 *
 *   For each primary (motherID == 0):
 *
 *     1. TRUE QUANTITIES. Its momentum, direction and the path it took
 *        through the gas, Sum(stepSize) over its own TPC deposits. The
 *        magnetic field runs along z, so the helix radius is set by the
 *        transverse momentum and the efficiency is plotted against both |p|
 *        and pT, as well as against cos(theta) and the true path length.
 *
 *     2. DENOMINATOR. Primaries whose true path in the gas reaches
 *        `minTrueLength` are the ones a tracker could reasonably be asked to
 *        find. Neutrals and particles that never really traversed the gas
 *        drop out here on their own, so what is left measures the
 *        reconstruction rather than the acceptance.
 *
 *     3. HIT MATCHING. A reconstructed hit belongs to the primary when the
 *        primary appears among its trackIDs with a trackFraction of at least
 *        `minHitFraction` -- i.e. when the primary's ionisation dominates the
 *        pulse that made it.
 *
 *     4. CANDIDATE. Clusters are ranked by how many matched hits they hold,
 *        and every cluster the primary dominates -- at least
 *        `minClusterHits` matched hits and a hit purity of at least
 *        `minClusterPurity` -- joins the candidate. A track broken in two by
 *        a gap is therefore still one candidate, and the number of clusters
 *        it took is recorded, so fragmentation stays visible rather than
 *        being hidden by the merge.
 *
 *     5. NUMERATOR. The candidate counts as a findable track when it has at
 *        least `minHits` hits, spans at least `minLength` in 3D, and is at
 *        least `minPurity` pure.
 *
 *   EXTENT
 *   The 3D span is the distance between the two most widely separated hits
 *   of the candidate, found with the usual two-pass farthest-point walk. For
 *   a straight track that is its length; for one that curls up in the field
 *   it is the chord across the spiral, which is the right thing to cut on --
 *   a tight curler covers little ground however many hits it leaves. The
 *   extent along the candidate's principal axis is written out alongside it.
 *
 *   SECONDARIES
 *   The true path length counts only the primary's own deposits, so that it
 *   is the length of the particle's track and not of the delta rays around
 *   it. The matching does the opposite and keeps them, because the
 *   digitisation books a delta ray's charge under the parent that made it
 *   and no detector could separate the two.
 *
 *   The output tree holds one entry per primary in the denominator, with the
 *   raw quantities behind every cut, so the working point can be moved
 *   afterwards without running the job again.
 *
 * Usage:
 *   The input has to have been through the TPC reconstruction, since it needs
 *   TPCHits and TPCClusters:
 *
 *     GArReconstruction -i gun.root -m macros/tpc_reco.mac -o tpc_reco.root
 *     GArAnalysis -a TrackingEfficiency.C -i tpc_reco.root -o tracking_eff.root
 *
 *   The working point is set by a job macro, so asking for longer tracks or a
 *   purer candidate is a matter of editing TrackingEfficiency.mac:
 *
 *     GArAnalysis -a TrackingEfficiency.C -i tpc_reco.root -o eff.root \
 *                 -m TrackingEfficiency.mac
 *
 * Parameters (see TrackingEfficiency.mac):
 *   /ana/minHits             hits the candidate needs to count as found
 *   /ana/minLength           3D span the candidate has to cover [cm]
 *   /ana/minPurity           fraction of the candidate's hits owned by the primary
 *   /ana/minTrueLength       true path in the gas needed to enter the denominator [cm]
 *   /ana/minHitFraction      share of a hit's true energy the primary needs to own it
 *   /ana/minClusterPurity    hit purity a cluster needs to join the candidate
 *   /ana/minClusterHits      matched hits a cluster needs to join the candidate
 *   /ana/excludeSecondaries  leave delta-ray deposits out of the true path length
 *
 ***************************************************************************/

// GArAnalysis has libGArAnalysis loaded before it compiles this macro, so
// there is no R__LOAD_LIBRARY here. To compile it by hand in ROOT instead
// (.L TrackingEfficiency.C+), start ROOT from the build directory or with its
// rootlogon.C, which is what loads the library there.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TPad.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "AnalysisBase.hh"
#include "AnalysisMath.hh"
#include "DigiDataTypes.hh"
#include "PlotStyle.hh"
#include "SimDataTypes.hh"

/* -------------------------------------------------------------------------- */
/*                              Unit conversions                              */
/* -------------------------------------------------------------------------- */

// The simulation stores momenta in MeV/c; the plots use GeV/c
const Double_t kMeVToGeV = 1.e-3;

/* -------------------------------------------------------------------------- */
/*                               Small helpers                                */
/* -------------------------------------------------------------------------- */

namespace {

// Distance between the two most widely separated points, by the two-pass
// farthest-point walk: the point furthest from an arbitrary start is an end of
// the diameter, and the point furthest from that is the other. Exact for the
// straight and gently curved sets of points a track leaves, and linear rather
// than quadratic in the number of hits.
Double_t Extent(const std::vector<TVector3>& points)
{
    if (points.size() < 2) return 0.;

    auto farthestFrom = [&points](const TVector3& from) {
        size_t best = 0;
        Double_t bestDistance = -1.;
        for (size_t i = 0; i < points.size(); ++i) {
            const Double_t distance = (points[i] - from).Mag();
            if (distance > bestDistance) { bestDistance = distance; best = i; }
        }
        return best;
    };

    const size_t a = farthestFrom(points[0]);
    const size_t b = farthestFrom(points[a]);
    return (points[b] - points[a]).Mag();
}

// Extent along the principal axis of the points, i.e. the length of the
// straight line they are strung out along. Shorter than Extent() for a track
// that bends, equal to it for one that does not.
Double_t PrincipalExtent(const std::vector<TVector3>& points)
{
    if (points.size() < 2) return 0.;

    TVector3 centroid;
    for (const TVector3& point : points) centroid += point;
    centroid *= 1. / points.size();

    // Covariance of the points about their centroid
    Double_t covariance[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    for (const TVector3& point : points) {
        const TVector3 d = point - centroid;
        for (Int_t a = 0; a < 3; ++a) {
            for (Int_t b = 0; b < 3; ++b) covariance[a][b] += d[a] * d[b];
        }
    }

    // Leading eigenvector, by power iteration
    TVector3 axis(1., 1., 1.);
    for (Int_t iteration = 0; iteration < 32; ++iteration) {
        TVector3 next;
        for (Int_t a = 0; a < 3; ++a) {
            for (Int_t b = 0; b < 3; ++b) next[a] += covariance[a][b] * axis[b];
        }
        if (next.Mag() <= 0.) break;
        axis = next.Unit();
    }

    Double_t low = 0., high = 0.;
    Bool_t first = kTRUE;
    for (const TVector3& point : points) {
        const Double_t along = (point - centroid).Dot(axis);
        if (first || along < low)  low = along;
        if (first || along > high) high = along;
        first = kFALSE;
    }
    return high - low;
}

} // namespace

/* -------------------------------------------------------------------------- */
/*                                  Analysis                                  */
/* -------------------------------------------------------------------------- */

// A denominator/numerator pair for one variable, binned over the range the
// sample turned out to cover.
//
// It lives outside the analysis class on purpose: GArAnalysis compiles this
// macro with ACLiC, which generates a ROOT dictionary for it, and a
// dictionary cannot be made for a std::vector of a type nested inside a
// class -- the generated code has no access to it. A type a member vector is
// built on therefore belongs at file scope.
struct EffPlot {
    TH1D* all = nullptr;
    TH1D* found = nullptr;
    std::string name;
    std::string axis;
    Bool_t logX = kFALSE;
};

class TrackingEfficiencyAnalysis : public ana::AnalysisBase {
protected:

    /* ------------------------------------------------------------------ */
    /* The working point, read from the job macro before anything is       */
    /* opened. A cut the macro does not mention keeps the value its member */
    /* is declared with under "Settings" at the bottom of the class.       */
    /* ------------------------------------------------------------------ */
    void Configure(const ana::ParameterSet& params) override
    {
        params.Get("minHits",            fMinHits);
        params.Get("minLength",          fMinLength);
        params.Get("minPurity",          fMinPurity);
        params.Get("minTrueLength",      fMinTrueLength);
        params.Get("minHitFraction",     fMinHitFraction);
        params.Get("minClusterPurity",   fMinClusterPurity);
        params.Get("minClusterHits",     fMinClusterHits);
        params.Get("excludeSecondaries", fExcludeSecondaries);
    }

    void BeginJob() override
    {
        // Both are needed: the clusters carry indices into the hit
        // collection, so one without the other says nothing
        fHits     = Require<std::vector<digi::TPCHit>>("TPCHits");
        fClusters = Require<std::vector<digi::TPCCluster>>("TPCClusters");

        // One entry per primary that ionised the gas
        if (TTree* out = Output()) {
            out->Branch("eventID",        &fEventID);
            out->Branch("trackID",        &fTrackID);
            out->Branch("pdgCode",        &fPdgCode);

            // Truth
            out->Branch("momentum",       &fMomentum);       // GeV/c
            out->Branch("momentumT",      &fMomentumT);      // GeV/c
            out->Branch("cosTheta",       &fCosTheta);       // w.r.t. z
            out->Branch("trueLength",     &fTrueLength);     // cm
            out->Branch("trueEdep",       &fTrueEdep);       // MeV
            out->Branch("nTrueHits",      &fNTrueHits);
            out->Branch("startZ",         &fStartZ);         // cm
            out->Branch("startR",         &fStartR);         // cm

            // What the reconstruction left behind
            out->Branch("nMatchedHits",   &fNMatchedHits);   // over the whole event
            out->Branch("nCandClusters",  &fNCandClusters);
            out->Branch("nCandHits",      &fNCandHits);
            out->Branch("nCandMatched",   &fNCandMatched);
            out->Branch("candLength",     &fCandLength);     // cm, 3D span
            out->Branch("candPcaLength",  &fCandPcaLength);  // cm, along the axis
            out->Branch("candPurity",     &fCandPurity);
            out->Branch("candCompleteness", &fCandCompleteness);
            out->Branch("candCharge",     &fCandCharge);     // electrons
            out->Branch("candEnergy",     &fCandEnergy);     // MeV

            // The leading cluster on its own, for comparison
            out->Branch("leadHits",       &fLeadHits);
            out->Branch("leadMatched",    &fLeadMatched);
            out->Branch("leadLength",     &fLeadLength);     // cm
            out->Branch("leadPurity",     &fLeadPurity);

            out->Branch("isFound",        &fIsFound);
        }

        // The efficiency histograms are booked in EndJob() instead, once the
        // momentum and length the sample actually covers are known: a gun
        // sample spans whatever range it was thrown over, and a fixed axis
        // would be mostly empty for one and clipped for the next
    }

    /* ------------------------------------------------------------------ */
    /* One event = one gun particle (or a handful of them)                 */
    /* ------------------------------------------------------------------ */
    void Run() override
    {
        if (!sim.IsValid()) return;

        fEventID = sim.eventID;

        for (size_t i = 0; i < sim.NParticles(); ++i) {

            if (sim.MotherID(i) != 0) continue;  // primaries only
            fNPrimaries++;

            fTrackID = sim.TrackID(i);
            fPdgCode = sim.PdgCode(i);

            /* -------------------- True path through the gas ------------ */

            fTrueLength = 0.;
            fTrueEdep = 0.;
            fNTrueHits = 0;

            for (size_t k : sim.TPCHitsOfTrack(fTrackID)) {
                // Delta rays folded back onto the parent wander away from it,
                // so they are no part of the parent's own path
                if (fExcludeSecondaries && sim.TPCHitIsSecondary(k)) continue;

                const root::TPCHit& hit = sim.TPCHit(k);
                fTrueLength += hit.stepSize;
                fTrueEdep += hit.energyDeposit;
                fNTrueHits++;
            } // end loop over the primary's TPC deposits

            // The denominator: primaries that actually crossed the gas
            if (fTrueLength < fMinTrueLength) {
                fNTooShort++;
                continue;
            }
            fNInDenominator++;

            /* ------------------------- Kinematics ---------------------- */

            const TVector3 momentum = sim.StartMomentum(i);
            fMomentum  = momentum.Mag() * kMeVToGeV;
            fMomentumT = momentum.Perp() * kMeVToGeV;   // field is along z
            fCosTheta  = (momentum.Mag() > 0.) ? momentum.CosTheta() : 0.;

            const TVector3 start = sim.StartPosition(i);
            fStartZ = start.Z();
            fStartR = start.Perp();

            /* --------------------- Match hits to the primary ----------- */

            // matched[h] is true when hit h is dominated by this primary
            std::vector<Bool_t> matched(fHits->size(), kFALSE);
            fNMatchedHits = 0;

            for (size_t h = 0; h < fHits->size(); ++h) {
                if (!IsMatched((*fHits)[h])) continue;
                matched[h] = kTRUE;
                fNMatchedHits++;
            } // end loop over reconstructed hits

            /* ------------- Collect the clusters this primary owns ------ */

            std::vector<TVector3> candidatePoints;
            fNCandClusters = 0;
            fNCandHits = 0;
            fNCandMatched = 0;
            fCandCharge = 0.;
            fCandEnergy = 0.;

            fLeadHits = 0;
            fLeadMatched = 0;
            fLeadLength = 0.;
            fLeadPurity = 0.;

            std::vector<TVector3> leadPoints;

            for (const digi::TPCCluster& cluster : *fClusters) {

                Int_t nMatched = 0;
                for (Int_t index : cluster.hitIndices) {
                    if (index >= 0 && index < static_cast<Int_t>(matched.size())
                        && matched[index]) nMatched++;
                }
                if (nMatched == 0) continue;

                const Int_t nTotal = static_cast<Int_t>(cluster.hitIndices.size());
                const Double_t purity = (nTotal > 0)
                                      ? static_cast<Double_t>(nMatched) / nTotal : 0.;

                // The best cluster is worth reporting on its own, whether or
                // not it is pure enough to be believed
                if (nMatched > fLeadMatched) {
                    fLeadMatched = nMatched;
                    fLeadHits = nTotal;
                    fLeadPurity = purity;
                    leadPoints = HitPositions(cluster);
                }

                // A cluster joins the candidate when this primary dominates
                // it; several may, so a track split by a gap survives whole
                if (nMatched < fMinClusterHits || purity < fMinClusterPurity) continue;

                fNCandClusters++;
                fNCandHits += nTotal;
                fNCandMatched += nMatched;
                fCandCharge += cluster.charge;
                fCandEnergy += cluster.energy;

                const std::vector<TVector3> points = HitPositions(cluster);
                candidatePoints.insert(candidatePoints.end(),
                                       points.begin(), points.end());
            } // end loop over clusters

            fLeadLength = Extent(leadPoints);

            fCandLength    = Extent(candidatePoints);
            fCandPcaLength = PrincipalExtent(candidatePoints);
            fCandPurity    = (fNCandHits > 0)
                           ? static_cast<Double_t>(fNCandMatched) / fNCandHits : 0.;
            fCandCompleteness = (fNMatchedHits > 0)
                              ? static_cast<Double_t>(fNCandMatched) / fNMatchedHits : 0.;

            /* ------------------ Is that enough for a track? ------------ */

            fIsFound = (fNCandHits >= fMinHits)
                    && (fCandLength >= fMinLength)
                    && (fCandPurity >= fMinPurity);

            if (fIsFound) {
                fNFound++;
                fNCandClustersTotal += fNCandClusters;
                fCompletenessTotal += fCandCompleteness;
            }

            // A track shorter than the span the candidate is asked to cover
            // cannot pass however well it was reconstructed, so the efficiency
            // over these is reported separately from the one over everything
            if (fTrueLength >= fMinLength) {
                fNTrackable++;
                if (fIsFound) fNFoundTrackable++;
            }

            // Kept so that the plots can be binned over the range the
            // sample covers; all four are parallel to fFound
            fMomenta.push_back(fMomentum);
            fMomentaT.push_back(fMomentumT);
            fCosThetas.push_back(fCosTheta);
            fTrueLengths.push_back(fTrueLength);
            fFound.push_back(fIsFound);

            Fill();  // one output entry per primary
        } // end loop over particles
    }

    /* ------------------------------------------------------------------ */
    /* Called once, after the event loop and before the output is written. */
    /* ------------------------------------------------------------------ */
    void EndJob() override
    {
        // EndJob() runs even when the job never started -- a required
        // product missing from the input aborts during BeginJob() -- and an
        // efficiency of zero over nothing would read as a result
        if (fNPrimaries == 0) {
            std::cout << "\nNo primaries were processed, so there is no "
                         "efficiency to report.\n" << std::endl;
            return;
        }

        const Double_t efficiency = (fNInDenominator > 0)
            ? static_cast<Double_t>(fNFound) / fNInDenominator : 0.;
        const Double_t efficiencyTrackable = (fNTrackable > 0)
            ? static_cast<Double_t>(fNFoundTrackable) / fNTrackable : 0.;

        std::cout << "\n=== Tracking efficiency ===\n"
                  << "Events processed:            " << NEvents() << "\n"
                  << "Primaries:                   " << fNPrimaries << "\n"
                  << "  under " << fMinTrueLength << " cm in the gas:     "
                  << fNTooShort << "\n"
                  << "Denominator (ionised gas):   " << fNInDenominator << "\n"
                  << "Found as tracks:             " << fNFound << "\n"
                  << "Efficiency:                  " << efficiency << "\n"
                  << std::endl;

        // A primary whose true path is shorter than the span demanded of the
        // candidate is out of reach whatever the reconstruction does, so the
        // first number above folds that in and the second one does not
        std::cout << "Of those whose true path reaches the " << fMinLength
                  << " cm the candidate\nis asked to span:\n"
                  << "  in range:                  " << fNTrackable << "\n"
                  << "  found:                     " << fNFoundTrackable << "\n"
                  << "  efficiency:                " << efficiencyTrackable << "\n"
                  << std::endl;

        if (fNFound > 0) {
            // Fragmentation: how many clusters it took to hold one track, and
            // how much of the track's charge those clusters actually gathered
            std::cout << "Candidates that were found:\n"
                      << "  clusters per track:        "
                      << static_cast<Double_t>(fNCandClustersTotal) / fNFound << "\n"
                      << "  completeness:              "
                      << fCompletenessTotal / fNFound << "\n"
                      << std::endl;
        }

        std::cout << "A candidate is found when it has >= " << fMinHits
                  << " hits, spans >= " << fMinLength << " cm and is >= "
                  << fMinPurity << " pure.\n"
                  << std::endl;

        fPlots.push_back(MakePlot(fMomenta, "P", "Momentum [GeV/c]", kTRUE));
        fPlots.push_back(MakePlot(fMomentaT, "Pt",
                                  "Transverse momentum [GeV/c]", kTRUE));
        fPlots.push_back(MakePlot(fCosThetas, "CosTheta", "cos#theta", kFALSE,
                                  -1., 1.));
        fPlots.push_back(MakePlot(fTrueLengths, "Length",
                                  "True track length in gas [cm]", kFALSE, 0.));

        for (const EffPlot& plot : fPlots) {
            MakeEfficiency(plot, Form("eff%s", plot.name.c_str()));
        }

        Draw();
    }

private:

    /* --------------------------- Truth matching ------------------------ */

    // Does this primary dominate the pulse the hit was made from?
    Bool_t IsMatched(const digi::TPCHit& hit) const
    {
        for (size_t t = 0; t < hit.trackIDs.size(); ++t) {
            if (hit.trackIDs[t] != fTrackID) continue;
            return t < hit.trackFractions.size()
                && hit.trackFractions[t] >= fMinHitFraction;
        }
        return kFALSE;
    }

    std::vector<TVector3> HitPositions(const digi::TPCCluster& cluster) const
    {
        std::vector<TVector3> points;
        points.reserve(cluster.hitIndices.size());
        for (Int_t index : cluster.hitIndices) {
            if (index < 0 || index >= static_cast<Int_t>(fHits->size())) continue;
            const digi::TPCHit& hit = (*fHits)[index];
            points.emplace_back(hit.x, hit.y, hit.z);
        }
        return points;
    }

    /* ------------------------- Histogram plumbing ---------------------- */

    // `low` and `high` pin an edge that should not follow the data, e.g. the
    // -1 to 1 of a cosine or the zero a length starts at
    EffPlot MakePlot(const std::vector<Double_t>& values,
                     const char* name, const char* axis, Bool_t logX,
                     Double_t low = std::numeric_limits<Double_t>::quiet_NaN(),
                     Double_t high = std::numeric_limits<Double_t>::quiet_NaN())
    {
        EffPlot plot;
        plot.name = name;
        plot.axis = axis;
        plot.logX = logX;

        if (values.empty()) return plot;

        const auto range = std::minmax_element(values.begin(), values.end());
        if (std::isnan(low))  low  = *range.first;
        if (std::isnan(high)) high = *range.second;

        // A logarithmic axis cannot start at zero, and a sample that covers
        // no range at all still has to produce a drawable histogram
        if (logX) low = std::max(low, 1.e-4);
        if (!(high > low)) high = low + (logX ? low : 1.);

        // Leave a little room, so the outermost entries are not on an edge
        if (logX) { low *= 0.9; high *= 1.1; }
        else      { const Double_t margin = 0.02 * (high - low);
                    low -= margin; high += margin; }

        const Int_t nBins = 20;
        const TString title = Form(";%s;Primaries", axis);

        if (logX) {
            const std::vector<Double_t> bins = ana::LogSpacedBins(nBins, low, high);
            plot.all   = new TH1D(Form("hAll%s", name),   title,
                                  nBins, bins.data());
            plot.found = new TH1D(Form("hFound%s", name), title,
                                  nBins, bins.data());
        } else {
            plot.all   = new TH1D(Form("hAll%s", name),   title, nBins, low, high);
            plot.found = new TH1D(Form("hFound%s", name), title, nBins, low, high);
        }

        for (size_t i = 0; i < values.size(); ++i) {
            plot.all->Fill(values[i]);
            if (fFound[i]) plot.found->Fill(values[i]);
        }
        return plot;
    }

    // TEfficiency carries the binomial uncertainties the ratio of two
    // histograms would lose
    void MakeEfficiency(const EffPlot& plot, const char* name)
    {
        if (!plot.all || !plot.found) return;
        if (!TEfficiency::CheckConsistency(*plot.found, *plot.all)) {
            std::cerr << "Warning: cannot build " << name
                      << ", the numerator is not a subset of the denominator"
                      << std::endl;
            return;
        }
        TEfficiency* efficiency = new TEfficiency(*plot.found, *plot.all);
        efficiency->SetNameTitle(name, Form("Tracking efficiency;%s;Efficiency",
                                            plot.axis.c_str()));
        efficiency->SetStatisticOption(TEfficiency::kFCP);  // Clopper-Pearson
        if (OutputFile()) efficiency->Write();
    }

    void Draw()
    {
        ana::SetPlotStyle();

        TCanvas canvas("cTrackingEfficiency", "Tracking efficiency", 1200, 900);
        canvas.Divide(2, 2);

        for (size_t i = 0; i < fPlots.size() && i < 4; ++i) {
            const EffPlot& plot = fPlots[i];
            if (!plot.all || !plot.found) continue;
            if (!TEfficiency::CheckConsistency(*plot.found, *plot.all)) continue;

            canvas.cd(static_cast<Int_t>(i) + 1);
            if (plot.logX) gPad->SetLogx();

            TEfficiency* efficiency = new TEfficiency(*plot.found, *plot.all);
            efficiency->SetTitle(Form(";%s;Efficiency", plot.axis.c_str()));
            efficiency->SetMarkerStyle(20);
            efficiency->SetLineWidth(2);
            efficiency->Draw("AP");

            gPad->Update();
            if (auto* graph = efficiency->GetPaintedGraph()) {
                graph->SetMinimum(0.);
                graph->SetMaximum(1.05);
                graph->GetXaxis()->SetLimits(plot.all->GetXaxis()->GetXmin(),
                                             plot.all->GetXaxis()->GetXmax());
                // A gun sample often spans less than a decade, which leaves a
                // logarithmic axis with no labelled tick at all
                if (plot.logX) {
                    graph->GetXaxis()->SetMoreLogLabels();
                    graph->GetXaxis()->SetNoExponent();
                }
            }
            gPad->Update();
        }

        canvas.SaveAs("tracking_efficiency.png");
    }

    /* ------------------------------ Settings --------------------------- */

    Int_t fMinHits = 20;                 // in the candidate
    Double_t fMinLength = 10.;           // cm, 3D span of the candidate
    Double_t fMinPurity = 0.5;           // of the candidate
    Double_t fMinTrueLength = 5.;        // cm, true path in the gas
    Double_t fMinHitFraction = 0.5;      // to own a hit
    Double_t fMinClusterPurity = 0.5;    // to own a cluster
    Int_t fMinClusterHits = 3;           // matched, to own a cluster
    Bool_t fExcludeSecondaries = kTRUE;  // in the true path length

    /* --------------------------- Input products ------------------------ */

    ana::Handle<std::vector<digi::TPCHit>> fHits;
    ana::Handle<std::vector<digi::TPCCluster>> fClusters;

    /* -------------------------- Output branches ------------------------ */

    Int_t fEventID = 0;
    Int_t fTrackID = 0;
    Int_t fPdgCode = 0;

    Double_t fMomentum = 0.;
    Double_t fMomentumT = 0.;
    Double_t fCosTheta = 0.;
    Double_t fTrueLength = 0.;
    Double_t fTrueEdep = 0.;
    Int_t fNTrueHits = 0;
    Double_t fStartZ = 0.;
    Double_t fStartR = 0.;

    Int_t fNMatchedHits = 0;
    Int_t fNCandClusters = 0;
    Int_t fNCandHits = 0;
    Int_t fNCandMatched = 0;
    Double_t fCandLength = 0.;
    Double_t fCandPcaLength = 0.;
    Double_t fCandPurity = 0.;
    Double_t fCandCompleteness = 0.;
    Double_t fCandCharge = 0.;
    Double_t fCandEnergy = 0.;

    Int_t fLeadHits = 0;
    Int_t fLeadMatched = 0;
    Double_t fLeadLength = 0.;
    Double_t fLeadPurity = 0.;

    Bool_t fIsFound = kFALSE;

    /* ------------------------------ Histograms ------------------------- */

    // One entry per primary in the denominator, all parallel to fFound, kept
    // so that EndJob() can bin them over the range the sample covers
    std::vector<Double_t> fMomenta;
    std::vector<Double_t> fMomentaT;
    std::vector<Double_t> fCosThetas;
    std::vector<Double_t> fTrueLengths;
    std::vector<Bool_t> fFound;

    std::vector<EffPlot> fPlots;

    /* ------------------------------- Counters -------------------------- */

    Long64_t fNPrimaries = 0;
    Long64_t fNTooShort = 0;
    Long64_t fNInDenominator = 0;
    Long64_t fNFound = 0;
    Long64_t fNTrackable = 0;
    Long64_t fNFoundTrackable = 0;
    Long64_t fNCandClustersTotal = 0;
    Double_t fCompletenessTotal = 0.;
};

/* -------------------------------------------------------------------------- */
/*                       What GArAnalysis runs from here                      */
/* -------------------------------------------------------------------------- */

// The output tree is named in TrackingEfficiency.mac, with /ana/global/outputTree
ANA_ANALYSIS(TrackingEfficiencyAnalysis)
