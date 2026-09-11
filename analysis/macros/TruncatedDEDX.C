 /***************************************************************************
 * TruncatedDEDX.C
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Created: 29/07/2026
 *
 * Description:
 *   Truncated mean dE/dx of primary particles in the gas TPC, plotted
 *   against several quantities. Intended for particle-gun samples.
 *
 *   For every primary particle (motherID == 0) the macro:
 *     0. optionally keeps only the tracks that leave the TPC gas volume,
 *        i.e. punch-through tracks, judged by where the track stops;
 *     1. collects its TPC hits, optionally dropping the deposits flagged as
 *        coming from secondaries (delta rays folded back onto the parent);
 *     2. forms one dE/dx sample per hit, edep / stepSize, in keV/cm;
 *     3. averages the lowest 60% of those samples (truncated mean);
 *     4. records it against four quantities:
 *          - the initial momentum        |p|          [GeV/c], log axis
 *          - the track length in the gas Sum(step)    [cm]
 *          - the kinetic energy on leaving the TPC    [GeV]
 *          - the kinetic energy lost inside the TPC   [GeV]
 *
 *   Note that the truncation is applied to the dE/dx samples, not to the raw
 *   energy deposits: Geant4 step lengths are not uniform, so the two orderings
 *   are not the same, and it is the dE/dx distribution whose Landau tail the
 *   truncation is meant to remove. To truncate on the raw deposits instead,
 *   sort on the hits' energyDeposit and average edep/step over the survivors.
 *
 *   HOW THE EXIT KINETIC ENERGY IS OBTAINED
 *   The last stored trajectory point is where the track finally stops -- in
 *   the calorimeter or beyond for a punch-through particle, not where it
 *   crosses the TPC boundary. The kinetic energy at the boundary is therefore
 *   not stored and is reconstructed from energy conservation:
 *
 *       KE_exit = KE_initial - (energy deposited in the gas)
 *
 *   summed over all of the track's TPC hits, delta rays included, since that
 *   energy also came out of the parent. KE_initial is built from the initial
 *   momentum and the PDG mass. The energy lost is then KE_initial - KE_exit.
 *   This ignores energy carried out of the gas by escaping secondaries and by
 *   radiated photons, so it is a lower bound on the loss; it is exact only to
 *   the extent that everything the particle lost was deposited locally. For
 *   the true boundary-crossing value the simulation would have to store the
 *   momentum at the last trajectory point inside the TPC.
 *
 *   The simulation stores energies in MeV and lengths in cm; this macro
 *   converts to the units conventionally used for dE/dx plots, GeV and
 *   keV/cm, on the way in. The output tree is in those units too.
 *
 * Usage:
 *   The macro is run by GArAnalysis, which compiles it and hands it the job:
 *
 *     GArAnalysis -a TruncatedDEDX.C -i sim.root -o dedx.root
 *
 *   The input accepts wildcards and comma-separated lists, and files under
 *   /pnfs are streamed over XRootD automatically:
 *
 *     GArAnalysis -a TruncatedDEDX.C -i '/pnfs/dune/scratch/users/me/gun_*.root' \
 *                 -o dedx.root
 *
 *   The parameters below come from a job macro, and only the ones that are to
 *   differ from the defaults have to appear in it:
 *
 *     GArAnalysis -a TruncatedDEDX.C -i sim.root -o dedx.root -m TruncatedDEDX.mac
 *
 * Parameters (see TruncatedDEDX.mac):
 *   /ana/truncation           fraction of the dE/dx samples kept, lowest first
 *   /ana/excludeSecondaries   drop the deposits booked as delta rays
 *   /ana/requireExit          keep only the tracks that leave the TPC
 *   /ana/minHits              samples a primary needs to be used at all
 *   /ana/pMin, /ana/pMax      momentum axis [GeV/c]; negative takes it from the data
 *   /ana/dedxMin, /ana/dedxMax  dE/dx axis [keV/cm]; negative takes it from the data
 *
 *   The dE/dx and momentum axes are the only ones that can be restricted; the
 *   track-length and kinetic-energy axes always follow the data.
 *
 ***************************************************************************/

// GArAnalysis has libGArAnalysis loaded before it compiles this macro, so
// there is no R__LOAD_LIBRARY here. To compile it by hand in ROOT instead
// (.L TruncatedDEDX.C+), start ROOT from the build directory or with its
// rootlogon.C, which is what loads the library there.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TParticlePDG.h"
#include "TProfile.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "AnalysisBase.hh"
#include "AnalysisMath.hh"
#include "PlotStyle.hh"
#include "SimDataTypes.hh"

/* -------------------------------------------------------------------------- */
/*                              Unit conversions                              */
/* -------------------------------------------------------------------------- */

// The simulation stores momenta in MeV/c, energy deposits in MeV and step
// sizes in cm, so dE/dx comes out in MeV/cm
const Double_t kMeVToGeV = 1.e-3;
const Double_t kMeVPerCmToKeVPerCm = 1.e3;

/* -------------------------------------------------------------------------- */
/*                                  Analysis                                  */
/* -------------------------------------------------------------------------- */

class TruncatedDEDXAnalysis : public ana::AnalysisBase {
protected:

    /* ------------------------------------------------------------------ */
    /* The parameters this analysis takes, read from the job macro before  */
    /* anything is opened. One the macro does not set keeps the value its  */
    /* member is declared with at the bottom of the class.                 */
    /* ------------------------------------------------------------------ */
    void Configure(const ana::ParameterSet& params) override
    {
        params.Get("truncation",         fTruncation);
        params.Get("excludeSecondaries", fExcludeSecondaries);
        params.Get("requireExit",        fRequireExit);
        params.Get("minHits",            fMinHits);
        params.Get("pMin",               fPMin);
        params.Get("pMax",               fPMax);
        params.Get("dedxMin",            fDEDXMin);
        params.Get("dedxMax",            fDEDXMax);

        // A truncation outside (0, 1] would silently become "keep one sample"
        // or "keep everything" inside ana::TruncatedMean, so say so instead
        if (fTruncation <= 0. || fTruncation > 1.) {
            std::cerr << "TruncatedDEDX: /ana/truncation has to be in (0, 1], not "
                      << fTruncation << std::endl;
            Abort();
        }
    }

    void BeginJob() override
    {
        // The punch-through selection needs the TPC dimensions, which are
        // known by now -- the geometry record is read before BeginJob() runs
        // -- so a job that cannot apply the requested selection is stopped
        // here rather than quietly producing an unselected plot.
        if (fRequireExit && !SetUpExitSelection()) {
            fSelectionUsable = kFALSE;
            Abort();
            return;
        }

        // One entry per primary particle
        if (TTree* out = Output()) {
            out->Branch("eventID",     &fEventID);
            out->Branch("pdgCode",     &fPdgCode);
            out->Branch("momentum",    &fMomentum);     // GeV/c
            out->Branch("dedx",        &fDEDX);         // keV/cm
            out->Branch("nHits",       &fNHits);        // samples available
            out->Branch("nHitsKept",   &fNHitsKept);    // samples averaged
            out->Branch("trackLength", &fTrackLength);  // cm
            out->Branch("keInitial",   &fKEInitial);    // GeV
            out->Branch("keExit",      &fKEExit);       // GeV, -1 if unknown
            out->Branch("keLoss",      &fKELoss);       // GeV, -1 if unknown
        }
    }

    /* ------------------------------------------------------------------ */
    /* One event = one gun particle (or a handful of them)                 */
    /* ------------------------------------------------------------------ */
    void Run() override
    {
        if (!sim.IsValid()) return;

        for (size_t i = 0; i < sim.NParticles(); ++i) {

            /* ------------------- Identify the primary ------------------ */

            if (sim.MotherID(i) != 0) continue;  // primaries only

            const Int_t id = sim.TrackID(i);

            /* ------------ Punch-through selection (optional) ----------- */

            if (fRequireExit) {
                // A particle with no stored trajectory has no known fate
                if (!sim.HasTrajectory(i)) {
                    fNNoEndPoint++;
                    continue;
                }
                if (!ExitsTPC(sim.EndPosition(i))) {
                    fNStoppedInTPC++;
                    continue;
                }
            }

            /* ------------- Collect this particle's dE/dx samples ------- */

            std::vector<Double_t> samples;
            Double_t trackLength = 0.;

            for (size_t k : sim.TPCHitsOfTrack(id)) {

                // Deposits attributed to unstored secondaries of this track,
                // i.e. delta rays. Excluded by default: they are what the
                // truncated mean is meant to be insensitive to.
                if (fExcludeSecondaries && sim.TPCHitIsSecondary(k)) continue;

                const root::TPCHit& hit = sim.TPCHit(k);

                const Double_t step = hit.stepSize;
                if (step <= 0.) continue;  // guard against zero-length steps

                // MeV / cm -> keV / cm
                samples.push_back(hit.energyDeposit / step * kMeVPerCmToKeVPerCm);
                trackLength += step;
            } // end loop over TPC hits of this track

            // Too few samples for the truncated mean to mean anything
            const Int_t nHits = static_cast<Int_t>(samples.size());
            if (nHits < fMinHits) {
                fNTooFewHits++;
                continue;
            }

            /* ------------------- Initial momentum ---------------------- */

            const TVector3 momentum = sim.StartMomentum(i);

            // The null vector, when the particle has no stored trajectory
            if (momentum.Mag() <= 0.) {
                fNNoMomentum++;
                continue;
            }

            /* --------------------- Truncated mean ---------------------- */

            // Averages the lowest `fTruncation` of the samples; sorts
            // `samples` in place
            size_t kept = 0;
            fDEDX = ana::TruncatedMean(samples, fTruncation, &kept);

            fEventID = sim.eventID;
            fPdgCode = sim.PdgCode(i);
            fMomentum = momentum.Mag() * kMeVToGeV;   // MeV/c -> GeV/c
            fNHits = nHits;
            fNHitsKept = static_cast<Int_t>(kept);
            fTrackLength = trackLength;

            /* -------------------- Kinetic energies --------------------- */

            // Energy the track left in the gas. Taken over all of its TPC
            // hits, delta rays included, because that energy came out of the
            // parent whether or not the samples above used it.
            const Double_t edepTPC = sim.TPCEdepOfTrack(id) * kMeVToGeV;

            fKEInitial = KineticEnergy(fPdgCode, fMomentum);
            if (fKEInitial >= 0.) {
                fKEExit = std::max(0., fKEInitial - edepTPC);
                fKELoss = fKEInitial - fKEExit;
            } else {
                // Unknown PDG code, e.g. a nucleus: no mass, so no kinetic
                // energy. The momentum and length plots are still filled.
                fKEExit = -1.;
                fKELoss = -1.;
                fNNoMass++;
            }

            fDEDXValues.push_back(fDEDX);
            fMomenta.push_back(fMomentum);
            fTrackLengths.push_back(fTrackLength);
            fKEExitValues.push_back(fKEExit);
            fKELossValues.push_back(fKELoss);
            fNPrimaries++;

            Fill();
        } // end loop over particles
    }

    /* ------------------------------------------------------------------ */
    /* Histograms are built here, once the range of the data is known      */
    /* ------------------------------------------------------------------ */
    void EndJob() override
    {
        // BeginJob() already explained why the job could not run
        if (!fSelectionUsable) return;

        std::cout << "\n=== Summary ===\n"
                  << "Primaries used:                " << fNPrimaries << "\n"
                  << "Rejected, fewer than " << fMinHits << " hits:   " << fNTooFewHits << "\n"
                  << "Rejected, no initial momentum: " << fNNoMomentum << "\n";

        if (fRequireExit) {
            std::cout << "Rejected, stopped in the TPC:  " << fNStoppedInTPC << "\n"
                      << "Rejected, no end point:        " << fNNoEndPoint << "\n";
        }
        if (fNNoMass > 0) {
            std::cout << "No PDG mass, no kinetic energy: " << fNNoMass << "\n";
        }

        std::cout << "Truncation fraction:           " << fTruncation << "\n"
                  << "Secondary deposits:            "
                  << (fExcludeSecondaries ? "excluded" : "included") << "\n"
                  << "Track selection:               "
                  << (fRequireExit ? "exiting the TPC only" : "all primaries") << "\n"
                  << std::endl;

        if (fDEDXValues.empty()) {
            std::cerr << "No primary particles passed the selection; nothing to plot."
                      << std::endl;
            return;
        }

        ana::SetPlotStyle();

        // The four views of the same truncated dE/dx measurement
        MakePlot({"P", "Momentum",
                  "Initial momentum [GeV/c]",
                  "truncated_dedx_vs_momentum.png",
                  kTRUE, fPMin, fPMax, &fMomenta});

        MakePlot({"Length", "TrackLength",
                  "Track length in the TPC [cm]",
                  "truncated_dedx_vs_length.png",
                  kFALSE, -1., -1., &fTrackLengths});

        MakePlot({"KEExit", "KineticEnergyExit",
                  "Kinetic energy leaving the TPC [GeV]",
                  "truncated_dedx_vs_ke_exit.png",
                  kFALSE, -1., -1., &fKEExitValues});

        MakePlot({"KELoss", "KineticEnergyLoss",
                  "Kinetic energy lost in the TPC [GeV]",
                  "truncated_dedx_vs_ke_loss.png",
                  kFALSE, -1., -1., &fKELossValues});
    }

private:

    /* ------------------------------------------------------------------ */
    /* One 2D view of the measurement: what to plot the dE/dx against      */
    /* ------------------------------------------------------------------ */
    struct PlotSpec {
        const char* tag;        // short name, used for the histograms and folder
        const char* title;      // plot title
        const char* axisTitle;  // x axis title, with units
        const char* imageName;  // file the canvas is saved to
        Bool_t logX;            // logarithmic x axis, and logarithmic binning
        Double_t userMin;       // negative -> take the range from the data
        Double_t userMax;
        const std::vector<Double_t>* values;  // parallel to fDEDXValues; < 0 skips
    };

    /* ------------------------------------------------------------------ */
    /* Build, write and draw the 2D histogram, its profile and its slices  */
    /* for one view                                                        */
    /* ------------------------------------------------------------------ */
    void MakePlot(const PlotSpec& spec)
    {
        // Entries for which this quantity is defined
        std::vector<Double_t> x, y;
        x.reserve(spec.values->size());
        y.reserve(spec.values->size());
        for (size_t i = 0; i < spec.values->size(); ++i) {
            if (spec.values->at(i) < 0.) continue;  // undefined for this particle
            x.push_back(spec.values->at(i));
            y.push_back(fDEDXValues[i]);
        }

        if (x.empty()) {
            std::cerr << "No entries with a defined " << spec.axisTitle
                      << "; skipping that plot." << std::endl;
            return;
        }

        /* ------------------------- Axis ranges ------------------------- */

        Double_t xLo = 0., xHi = 0., yLo = 0., yHi = 0.;
        if (!DetermineRanges(x, y, spec.userMin, spec.userMax, spec.logX,
                             xLo, xHi, yLo, yHi, spec.axisTitle)) {
            return;
        }

        std::cout << spec.title << ": x from " << xLo << " to " << xHi
                  << (spec.logX ? "  (logarithmic)" : "")
                  << ", dE/dx from " << yLo << " to " << yHi << " keV/cm"
                  << std::endl;

        const std::vector<Double_t> xBins = BinEdges(xLo, xHi, spec.logX);
        if (xBins.empty()) {
            std::cerr << "Could not build bins for " << spec.axisTitle << std::endl;
            return;
        }

        // Attach the histograms to the output file
        if (OutputFile()) OutputFile()->cd();

        /* -------------------------- Histograms ------------------------- */

        const TString plotTitle = TString::Format(
            "%s;%s;Truncated mean dE/dx [keV/cm]", spec.title, spec.axisTitle);

        TH2F* hDEDX = new TH2F(TString::Format("hDEDXvs%s", spec.tag).Data(),
                               plotTitle.Data(),
                               kNXBins, xBins.data(),
                               kNDEDXBins, yLo, yHi);

        // Passing the dE/dx range to the profile as well, so that entries
        // outside the visible region are left out of the mean instead of
        // silently pulling it
        TProfile* pDEDX = new TProfile(TString::Format("pDEDXvs%s", spec.tag).Data(),
                                       plotTitle.Data(),
                                       kNXBins, xBins.data(), yLo, yHi);

        for (size_t i = 0; i < x.size(); ++i) {
            hDEDX->Fill(x[i], y[i]);
            pDEDX->Fill(x[i], y[i]);
        }

        WriteSlices(hDEDX, spec);

        /* ---------------------------- Drawing -------------------------- */

        TCanvas canvas(TString::Format("cDEDXvs%s", spec.tag).Data(),
                       spec.title, 800, 600);
        if (spec.logX) canvas.SetLogx();
        canvas.SetRightMargin(0.15);  // room for the colour scale

        hDEDX->Draw("colz");

        pDEDX->SetLineColor(kRed);
        pDEDX->SetLineWidth(2);
        pDEDX->SetMarkerColor(kRed);
        pDEDX->SetMarkerStyle(20);
        pDEDX->SetMarkerSize(0.8);
        pDEDX->Draw("same");

        canvas.SaveAs(spec.imageName);
    }

    /* ------------------------------------------------------------------ */
    /* The 1D dE/dx distribution of each x bin, as Y projections of the 2D */
    /* histogram, in a folder of their own                                 */
    /* ------------------------------------------------------------------ */
    void WriteSlices(TH2F* hDEDX, const PlotSpec& spec)
    {
        if (!OutputFile()) return;

        const TString dirName = TString::Format("%s_slices", spec.tag);
        TDirectory* sliceDir = OutputFile()->mkdir(dirName.Data());
        if (!sliceDir) {
            std::cerr << "Could not create the '" << dirName
                      << "' directory in the output file" << std::endl;
            return;
        }

        // Histograms are attached to the current directory, so the
        // projections end up inside the folder and are written out with the
        // rest of the file
        sliceDir->cd();

        Int_t nNonEmpty = 0;
        for (Int_t bin = 1; bin <= kNXBins; ++bin) {
            const Double_t lo = hDEDX->GetXaxis()->GetBinLowEdge(bin);
            const Double_t hi = hDEDX->GetXaxis()->GetBinUpEdge(bin);

            const TString name = TString::Format("hDEDX_%sbin%02d", spec.tag, bin);
            const TString title = TString::Format(
                "Truncated mean dE/dx, %.3g < %s < %.3g;"
                "Truncated mean dE/dx [keV/cm];Particles",
                lo, spec.axisTitle, hi);

            // Single x bin -> the dE/dx distribution of that bin
            TH1D* slice = hDEDX->ProjectionY(name.Data(), bin, bin);
            slice->SetTitle(title.Data());

            if (slice->GetEntries() > 0) nNonEmpty++;
        }

        std::cout << "    " << kNXBins << " slices in '" << dirName << "/' ("
                  << nNonEmpty << " non-empty)" << std::endl;

        OutputFile()->cd();
    }

    /* ------------------------------------------------------------------ */
    /* Bin edges spanning [lo, hi], logarithmic or uniform                 */
    /* ------------------------------------------------------------------ */
    static std::vector<Double_t> BinEdges(Double_t lo, Double_t hi, Bool_t logarithmic)
    {
        if (logarithmic) return ana::LogSpacedBins(kNXBins, lo, hi);

        std::vector<Double_t> edges(kNXBins + 1);
        const Double_t step = (hi - lo) / kNXBins;
        for (Int_t i = 0; i <= kNXBins; ++i) edges[i] = lo + i * step;
        return edges;
    }

    /* ------------------------------------------------------------------ */
    /* Resolve the plot ranges from the settings and, where a setting was  */
    /* left automatic, from the data. Returns kFALSE if no usable range    */
    /* could be established.                                              */
    /* ------------------------------------------------------------------ */
    Bool_t DetermineRanges(const std::vector<Double_t>& xValues,
                           const std::vector<Double_t>& yValues,
                           Double_t xUserMin, Double_t xUserMax,
                           Bool_t logX,
                           Double_t& xLo, Double_t& xHi,
                           Double_t& yLo, Double_t& yHi,
                           const char* xLabel) const
    {
        // Extent of `values` over the entries whose companion quantity lies
        // within [otherLo, otherHi]. A negative bound means "unbounded", so
        // restricting one axis narrows the automatic range of the other
        // rather than leaving it stretched by entries that are not drawn.
        auto extent = [](const std::vector<Double_t>& values,
                         const std::vector<Double_t>& other,
                         Double_t otherLo, Double_t otherHi,
                         Double_t& lo, Double_t& hi) -> Bool_t {
            Bool_t found = kFALSE;
            for (size_t i = 0; i < values.size(); ++i) {
                if (otherLo >= 0. && other[i] < otherLo) continue;
                if (otherHi >= 0. && other[i] > otherHi) continue;
                if (!found || values[i] < lo) lo = values[i];
                if (!found || values[i] > hi) hi = values[i];
                found = kTRUE;
            }
            return found;
        };

        // The x quantity, narrowed to the requested dE/dx window
        xLo = xUserMin;
        xHi = xUserMax;
        if (xUserMin < 0. || xUserMax < 0.) {
            Double_t lo = 0., hi = 0.;
            if (!extent(xValues, yValues, fDEDXMin, fDEDXMax, lo, hi)) {
                std::cerr << "No entries inside the requested dE/dx range ["
                          << fDEDXMin << ", " << fDEDXMax << "] keV/cm; "
                          << "cannot set the " << xLabel << " range automatically."
                          << std::endl;
                return kFALSE;
            }
            if (logX) {
                // Multiplicative padding, so the margin is even on a log axis
                if (xUserMin < 0.) xLo = lo * 0.9;
                if (xUserMax < 0.) xHi = hi * 1.1;
            } else {
                // Additive padding, not letting a non-negative quantity go
                // below zero
                const Double_t margin = 0.05 * (hi - lo);
                if (xUserMin < 0.) xLo = (lo >= 0.) ? std::max(0., lo - margin) : lo - margin;
                if (xUserMax < 0.) xHi = hi + margin;
                // A single-valued quantity would give an empty range
                if (xHi <= xLo) xHi = xLo + 1.;
            }
        }

        // dE/dx, narrowed to the x window just settled on.
        // An energy-loss axis conventionally starts at zero.
        yLo = (fDEDXMin < 0.) ? 0. : fDEDXMin;
        yHi = fDEDXMax;
        if (fDEDXMax < 0.) {
            Double_t lo = 0., hi = 0.;
            if (!extent(yValues, xValues, xLo, xHi, lo, hi)) {
                std::cerr << "No entries inside the " << xLabel << " range ["
                          << xLo << ", " << xHi << "]; "
                          << "cannot set the dE/dx range automatically." << std::endl;
                return kFALSE;
            }
            yHi = hi * 1.1;
        }

        // A logarithmic axis cannot reach zero
        if (logX && xLo <= 0.) {
            std::cerr << "The " << xLabel << " axis is logarithmic and must start above "
                      << "zero; got " << xLo << std::endl;
            return kFALSE;
        }
        if (xHi <= xLo) {
            std::cerr << "Empty " << xLabel << " range: " << xLo << " to " << xHi
                      << std::endl;
            return kFALSE;
        }
        if (yHi <= yLo) {
            std::cerr << "Empty dE/dx range: " << yLo << " to " << yHi << " keV/cm"
                      << std::endl;
            return kFALSE;
        }

        return kTRUE;
    }

    /* ------------------------------------------------------------------ */
    /* Kinetic energy in GeV from the momentum in GeV/c and the PDG mass.  */
    /* Returns -1 for a code the PDG database does not know, such as a     */
    /* nucleus.                                                            */
    /* ------------------------------------------------------------------ */
    static Double_t KineticEnergy(Int_t pdgCode, Double_t momentum)
    {
        TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
        if (!particle) return -1.;

        const Double_t mass = particle->Mass();  // GeV
        return std::sqrt(momentum * momentum + mass * mass) - mass;
    }

    /* ------------------------------------------------------------------ */
    /* Prepare the punch-through selection. Returns kFALSE, having said     */
    /* why, if the input cannot support it.                                */
    /* ------------------------------------------------------------------ */
    Bool_t SetUpExitSelection()
    {
        if (!HasGeometry()) {
            std::cerr << "Punch-through selection requested, but the input has no "
                      << "geometry tree, so the TPC dimensions are unknown."
                      << std::endl;
            return kFALSE;
        }

        if (geo.geometry_type != ana::kGArLike) {
            std::cerr << "Punch-through selection requested, but the geometry is not "
                      << "GArLike; there is no gas TPC cylinder to leave." << std::endl;
            return kFALSE;
        }

        fTPCRadius = geo.gar_tpc_radius;
        // The geometry tree stores the full length of the gas volume, which is
        // centred on the origin
        fTPCHalfLength = 0.5 * geo.gar_tpc_length;

        if (fTPCRadius <= 0. || fTPCHalfLength <= 0.) {
            std::cerr << "Punch-through selection requested, but the TPC dimensions "
                      << "read from the geometry tree are not usable: radius "
                      << fTPCRadius << " cm, length " << geo.gar_tpc_length << " cm."
                      << std::endl;
            return kFALSE;
        }

        std::cout << "Keeping only tracks that leave the TPC gas volume, i.e. that stop "
                  << "at r > " << fTPCRadius << " cm or |z| > " << fTPCHalfLength
                  << " cm\n" << std::endl;

        return kTRUE;
    }

    // Whether a track ended outside the TPC gas cylinder. Positions are global
    // and in cm, and the gas volume is centred on the origin.
    Bool_t ExitsTPC(const TVector3& end) const
    {
        const Double_t r = std::hypot(end.X(), end.Y());
        return (r > fTPCRadius) || (std::fabs(end.Z()) > fTPCHalfLength);
    }

    // Binning of the summary plots
    static constexpr Int_t kNXBins = 60;
    static constexpr Int_t kNDEDXBins = 100;

    // Settings
    Double_t fTruncation = 0.6;             // keep the lowest 60% of the samples
    Bool_t fExcludeSecondaries = kTRUE;     // drop delta-ray deposits
    Bool_t fRequireExit = kFALSE;           // keep only tracks that leave the TPC
    Int_t fMinHits = 10;                    // minimum samples per particle
    Double_t fPMin = -1., fPMax = -1.;      // GeV/c,  negative -> from the data
    Double_t fDEDXMin = -1., fDEDXMax = -1.;  // keV/cm, negative -> from the data

    // TPC gas volume, filled in when the punch-through selection is set up
    Double_t fTPCRadius = 0.;      // cm
    Double_t fTPCHalfLength = 0.;  // cm
    Bool_t fSelectionUsable = kTRUE;

    // Output branch variables
    Int_t fEventID = 0;
    Int_t fPdgCode = 0;
    Double_t fMomentum = 0.;     // GeV/c
    Double_t fDEDX = 0.;         // keV/cm
    Int_t fNHits = 0;
    Int_t fNHitsKept = 0;
    Double_t fTrackLength = 0.;  // cm
    Double_t fKEInitial = -1.;   // GeV
    Double_t fKEExit = -1.;      // GeV
    Double_t fKELoss = -1.;      // GeV

    // Kept so the histogram ranges can be taken from the data. All parallel
    // to fDEDXValues; a negative entry means the quantity is undefined for
    // that particle and is left out of its plot.
    std::vector<Double_t> fDEDXValues;
    std::vector<Double_t> fMomenta;
    std::vector<Double_t> fTrackLengths;
    std::vector<Double_t> fKEExitValues;
    std::vector<Double_t> fKELossValues;

    // Counters
    Long64_t fNPrimaries = 0;
    Long64_t fNTooFewHits = 0;
    Long64_t fNNoMomentum = 0;
    Long64_t fNStoppedInTPC = 0;
    Long64_t fNNoEndPoint = 0;
    Long64_t fNNoMass = 0;
};

/* -------------------------------------------------------------------------- */
/*                       What GArAnalysis runs from here                      */
/* -------------------------------------------------------------------------- */

// The output tree is named in TruncatedDEDX.mac, with /ana/global/outputTree
ANA_ANALYSIS(TruncatedDEDXAnalysis)
