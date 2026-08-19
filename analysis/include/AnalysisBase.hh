 /***************************************************************************
 * AnalysisBase.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Base class for FastGArSim analyses. It owns all of the boilerplate that
 *   used to be copy-pasted into every analysis macro: opening the input
 *   files, locating the trees, wiring up every branch, driving the event
 *   loop, reporting progress and writing the output file.
 *
 *   A concrete analysis only has to derive from AnalysisBase and implement
 *   Run(), which is called once per event with all input variables already
 *   loaded. Two optional hooks, BeginJob() and EndJob(), are available for
 *   booking output branches/histograms and for producing summaries.
 *
 *   Usage from a ROOT macro:
 *
 *     #include "AnalysisBase.hh"
 *
 *     class MyAnalysis : public ana::AnalysisBase {
 *     protected:
 *         void BeginJob() override { Output()->Branch("nhits", &fNHits); }
 *         void Run() override      { fNHits = sim.NECalHits(); Fill(); }
 *     private:
 *         Int_t fNHits = 0;
 *     };
 *
 *     void my_analysis(const char* simFiles, const char* outFile) {
 *         ana::AnalysisConfig cfg;
 *         cfg.simFiles = simFiles;   // path, wildcard or comma-separated list
 *         cfg.outputFile = outFile;
 *         MyAnalysis().Execute(cfg);
 *     }
 *
 ***************************************************************************/

#ifndef AnalysisBase_hh
#define AnalysisBase_hh

#include <string>
#include <vector>

#include "Rtypes.h"

#include "AnalysisEvent.hh"

class TChain;
class TFile;
class TTree;

namespace ana {

/* -------------------------------------------------------------------------- */
/*                             Job configuration                              */
/* -------------------------------------------------------------------------- */

struct AnalysisConfig {
    // FastGArSim flat ntuples, as produced by EventToNtupleConverter.C.
    // Required; must contain the analysis tree, and normally the geometry tree.
    //
    // This is an input specification rather than a single path: a
    // comma-separated list whose entries may be plain paths or shell
    // wildcard patterns. The files are chained, and the geometry is taken
    // from the first of them. See ana::ExpandInputFiles.
    std::string simFiles;

    // Output ROOT file. Leave empty to run without producing one (for
    // analyses that only print or draw).
    std::string outputFile;

    // GENIE gst files holding the truth record of the simulated events.
    // Optional: leave empty for particle-gun samples. When set, the entry
    // matching SimEvent::eventID is loaded before every call to Run().
    // Accepts the same wildcards and lists as simFiles, but note that the
    // eventID has to be unique across the whole chain for the lookup to be
    // meaningful.
    std::string genieFiles;

    // Rewrite resolved /pnfs paths as XRootD URLs, so that dCache files are
    // streamed instead of being read through the NFS mount. Wildcards are
    // still expanded against the mount, which therefore has to be visible.
    Bool_t xrootdForPnfs = kTRUE;

    // Tree names
    std::string simTreeName    = "AnaTree";
    std::string geoTreeName    = "GeoTree";
    std::string genieTreeName  = "gst";
    std::string outputTreeName = "AnaOutput";

    // Event range: process `maxEvents` entries starting at `firstEvent`.
    // maxEvents < 0 means "all remaining entries".
    Long64_t firstEvent = 0;
    Long64_t maxEvents  = -1;

    // Number of progress messages printed over the whole job
    Int_t nReports = 10;

    // Print the detector configuration at start-up
    Bool_t printGeometry = kTRUE;
};

/* -------------------------------------------------------------------------- */
/*                               Analyzer base                                */
/* -------------------------------------------------------------------------- */

class AnalysisBase {
public:
    AnalysisBase();
    virtual ~AnalysisBase();

    AnalysisBase(const AnalysisBase&) = delete;
    AnalysisBase& operator=(const AnalysisBase&) = delete;

    // Run the whole job: open inputs, BeginJob(), loop calling Run() once per
    // event, EndJob(), write and close the output. Returns kFALSE if the job
    // could not be set up.
    Bool_t Execute(const AnalysisConfig& config);

    // Shorthand for the common case
    Bool_t Execute(const char* simFileSpec,
                   const char* outputFileName = "",
                   const char* genieFileSpec = "");

protected:
    /* ---------------------------- Hooks to implement ---------------------- */

    // Called once, after the inputs are open and the output file exists.
    // Book output branches and histograms here.
    virtual void BeginJob() {}

    // Called once per selected event, with sim (and genie, when available)
    // already loaded. This is the only method an analysis must implement.
    virtual void Run() = 0;

    // Called once, after the event loop and before the output is written.
    virtual void EndJob() {}

    /* ------------------------------- Input data --------------------------- */

    SimEvent     sim;    // FastGArSim AnaTree, reloaded every event
    GenieEvent   genie;  // GENIE gst record, reloaded every event (if present)
    GeometryInfo geo;    // Detector configuration, loaded once

    /* --------------------------- Output and context ----------------------- */

    // Output tree; nullptr when AnalysisConfig::outputFile is empty
    TTree* Output() const { return fOutputTree; }

    // Output file; histograms created after BeginJob() starts are attached to
    // it automatically and written out at the end of the job
    TFile* OutputFile() const { return fOutputFile; }

    // Fill one entry of the output tree (no-op without an output tree)
    void Fill();

    Bool_t HasGenie() const { return fGenieTree != nullptr; }
    Bool_t HasGeometry() const { return fGeoTree != nullptr; }

    // Number of input files the analysis tree was chained from
    size_t NInputFiles() const { return fSimFileNames.size(); }
    const std::vector<std::string>& InputFiles() const { return fSimFileNames; }

    // Entry number of the event currently being processed, and the number of
    // events this job will process in total
    Long64_t CurrentEntry() const { return fCurrentEntry; }
    Long64_t NEvents() const { return fNEvents; }

    const AnalysisConfig& Config() const { return fConfig; }

    // Stop the event loop early. Called from Run() it ends the loop after the
    // current event; called from BeginJob() the loop is skipped entirely, so
    // it doubles as "this job cannot run". EndJob() is invoked either way.
    void Abort() { fAbort = kTRUE; }

private:
    Bool_t Initialize();
    void ProcessEvents();
    void Finalize();
    void CloseInputs();

    AnalysisConfig fConfig;

    // Kept open for the geometry tree, which is read from the first input file
    TFile* fGeoFile = nullptr;
    TFile* fOutputFile = nullptr;

    // The input trees are chains, so that a job can span many files
    TChain* fSimTree   = nullptr;
    TChain* fGenieTree = nullptr;
    TTree* fGeoTree    = nullptr;
    TTree* fOutputTree = nullptr;

    std::vector<std::string> fSimFileNames;

    Long64_t fCurrentEntry = 0;
    Long64_t fNEvents = 0;
    Bool_t fAbort = kFALSE;
};

} // namespace ana

#endif
