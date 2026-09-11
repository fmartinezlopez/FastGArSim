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
 *   loaded. Three optional hooks are available: Configure() for reading the
 *   analysis's own parameters out of the job macro, BeginJob() for booking
 *   output branches and histograms, and EndJob() for summaries and plots.
 *
 *   The input is read as it was written: the simulation's root::Event objects
 *   from the Events tree, and the reconstruction's products from the Reco
 *   tree, with no flat ntuple stage in between. Since the reconstruction is
 *   modular, which products are there depends on how it was run, so they are
 *   asked for by name in BeginJob() -- see Require() and Optional() below.
 *
 *   Usage from a ROOT macro:
 *
 *     #include "AnalysisBase.hh"
 *
 *     class MyAnalysis : public ana::AnalysisBase {
 *     protected:
 *         void Configure(const ana::ParameterSet& p) override {
 *             p.Get("minHits", fMinHits);   // /ana/minHits in the job macro
 *         }
 *         void BeginJob() override {
 *             fClusters = Require<std::vector<digi::TPCCluster>>("TPCClusters");
 *             Output()->Branch("nhits", &fNHits);
 *         }
 *         void Run() override {
 *             fNHits = sim.NECalHits();
 *             for (const digi::TPCCluster& c : *fClusters) { ... }
 *             Fill();
 *         }
 *     private:
 *         ana::Handle<std::vector<digi::TPCCluster>> fClusters;
 *         Int_t fMinHits = 10;
 *         Int_t fNHits = 0;
 *     };
 *
 *     ANA_ANALYSIS(MyAnalysis)
 *
 *   which GArAnalysis then compiles and runs:
 *
 *     GArAnalysis -a MyAnalysis.C -i sim.root -o out.root -m my.mac
 *
 ***************************************************************************/

#ifndef AnalysisBase_hh
#define AnalysisBase_hh

#include <string>
#include <vector>

#include "Rtypes.h"

#include "ProductSchema.hh"

#include "AnalysisEvent.hh"
#include "ParameterSet.hh"
#include "ProductStore.hh"

class TChain;
class TFile;
class TTree;

namespace ana {

/* -------------------------------------------------------------------------- */
/*                             Job configuration                              */
/* -------------------------------------------------------------------------- */

struct AnalysisConfig {
    // FastGArSim ROOT files: simulation output, or reconstruction output,
    // which by default also carries the simulation trees. Required.
    //
    // This is an input specification rather than a single path: a
    // comma-separated list whose entries may be plain paths or shell
    // wildcard patterns. The files are chained, and the geometry is taken
    // from the first of them. See ana::ExpandInputFiles.
    std::string inputFiles;

    // Output ROOT file. Leave empty to run without producing one (for
    // analyses that only print or draw).
    std::string outputFile;

    // GENIE gst files holding the truth record of the simulated events.
    // Optional: leave empty for particle-gun samples. When set, the entry
    // matching SimEvent::eventID is loaded before every call to Run().
    // Accepts the same wildcards and lists as inputFiles, but note that the
    // eventID has to be unique across the whole chain for the lookup to be
    // meaningful.
    std::string genieFiles;

    // Rewrite resolved /pnfs paths as XRootD URLs, so that dCache files are
    // streamed instead of being read through the NFS mount. Wildcards are
    // still expanded against the mount, which therefore has to be visible.
    Bool_t xrootdForPnfs = kTRUE;

    // Tree names
    std::string simTreeName    = "Events";
    std::string recoTreeName   = "Reco";
    std::string geoTreeName    = "Geometry";
    std::string genieTreeName  = "gst";
    std::string outputTreeName = "AnaOutput";

    // Branch of simTreeName holding the simulated event
    std::string simBranchName  = "Event";

    // Stop with an error when the input has no simulation tree. Analyses that
    // only read reconstruction products can set this false and run over files
    // written with /reco/global/copyInput false.
    Bool_t requireSim = kTRUE;

    // List the reconstruction products the input holds at start-up
    Bool_t printProducts = kTRUE;

    // Event range: process `maxEvents` entries starting at `firstEvent`.
    // maxEvents < 0 means "all remaining entries".
    Long64_t firstEvent = 0;
    Long64_t maxEvents  = -1;

    // Number of progress messages printed over the whole job
    Int_t nReports = 10;

    // Print the detector configuration at start-up
    Bool_t printGeometry = kTRUE;

    // The analysis's own parameters -- everything the job macro set that is
    // not one of the fields above. They are handed to Configure() before the
    // job starts; see ParameterSet.hh.
    ParameterSet params;
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
    Bool_t Execute(const char* inputFileSpec,
                   const char* outputFileName = "",
                   const char* genieFileSpec = "");

protected:
    /* ---------------------------- Hooks to implement ---------------------- */

    // Called once, before anything is opened, with the parameters the job
    // macro set. Read them into the analysis's own members here:
    //
    //     void Configure(const ana::ParameterSet& p) override {
    //         p.Get("minHits", fMinHits);
    //     }
    //
    // A parameter the macro did not set leaves the member alone, so its
    // declared value is the default. Nothing else has happened yet, so a
    // parameter that makes no sense can be rejected with Abort() before a
    // single file is opened.
    virtual void Configure(const ParameterSet& params) { (void)params; }

    // Called once, after the inputs are open and the output file exists.
    // Book output branches and histograms here.
    virtual void BeginJob() {}

    // Called once per selected event, with sim (and genie, when available)
    // already loaded. This is the only method an analysis must implement.
    virtual void Run() = 0;

    // Called once, after the event loop and before the output is written.
    virtual void EndJob() {}

    /* ------------------------------- Input data --------------------------- */

    SimEvent     sim;    // FastGArSim Events tree, reloaded every event
    GenieEvent   genie;  // GENIE gst record, reloaded every event (if present)
    GeometryInfo geo;    // Detector configuration, loaded once

    // Reconstruction products, whatever this file happens to hold. Use the
    // Require()/Optional() helpers below rather than binding through this
    // directly; it is here for Print(), Has() and Keys().
    ProductStore reco;

    /* --------------------------- Reconstruction products ------------------ */

    // Ask for a reconstruction product by name. Call these from BeginJob():
    // the handles they return stay valid for the whole job and are refilled
    // every event.
    //
    // Require() ends the job before the first event if the product is not
    // there, or is there with a different type, saying what the files do
    // hold. Optional() returns a handle that is simply never valid, so the
    // analysis can carry on without that product:
    //
    //     if (fWaveforms) { ... }
    //
    template <class T>
    Handle<T> Require(const std::string& key) { return reco.Bind<T>(key, kTRUE); }

    template <class T>
    Handle<T> Optional(const std::string& key) { return reco.Bind<T>(key, kFALSE); }

    // Whether the input holds a reconstruction tree at all
    Bool_t HasReco() const { return fRecoTree != nullptr; }

    // What the first input file says it contains, and how it was produced.
    // Files written before the schema record existed describe themselves by
    // their trees instead, so this is never empty for a readable file.
    const fastgarsim::ProductSchema& Schema() const { return fSchema; }

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

    // Number of input files the trees were chained from
    size_t NInputFiles() const { return fInputFileNames.size(); }
    const std::vector<std::string>& InputFiles() const { return fInputFileNames; }

    // Entry number of the event currently being processed, and the number of
    // events this job will process in total
    Long64_t CurrentEntry() const { return fCurrentEntry; }
    Long64_t NEvents() const { return fNEvents; }

    const AnalysisConfig& Config() const { return fConfig; }

    // The parameters this job was configured with, for an analysis that wants
    // to read one outside Configure()
    const ParameterSet& Params() const { return fConfig.params; }

    // Stop the event loop early. Called from Run() it ends the loop after the
    // current event; called from BeginJob() the loop is skipped entirely, so
    // it doubles as "this job cannot run". EndJob() is invoked either way.
    void Abort() { fAbort = kTRUE; }

private:
    // Report unreadable and unasked-for parameters after Configure(); false
    // when the job cannot go ahead
    Bool_t CheckParameters();

    Bool_t Initialize();
    Bool_t OpenInputs();
    // Report the products BeginJob() asked for and did not get, and stop the
    // job. Returns false when something was missing.
    Bool_t CheckRequestedProducts();
    void ProcessEvents();
    void Finalize();
    void CloseInputs();

    AnalysisConfig fConfig;

    // Kept open for the geometry tree and the schema, both read from the
    // first input file
    TFile* fFirstFile = nullptr;
    TFile* fOutputFile = nullptr;

    // The input trees are chains, so that a job can span many files
    TChain* fSimTree   = nullptr;
    TChain* fRecoTree  = nullptr;
    TChain* fGenieTree = nullptr;
    TTree* fGeoTree    = nullptr;
    TTree* fOutputTree = nullptr;

    fastgarsim::ProductSchema fSchema;

    std::vector<std::string> fInputFileNames;

    Long64_t fCurrentEntry = 0;
    Long64_t fNEvents = 0;
    Int_t fRecoTreeNumber = -1;
    Bool_t fAbort = kFALSE;
};

} // namespace ana

/* -------------------------------------------------------------------------- */
/*                          Making the analysis runnable                      */
/* -------------------------------------------------------------------------- */

// Every analysis macro ends with this line, naming its analysis class:
//
//     ANA_ANALYSIS(TruncatedDEDXAnalysis)
//
// GArAnalysis compiles the macro and calls the function it defines to get the
// analysis object; everything else about the job -- which files to read, where
// to write and what the analysis's own parameters are -- comes from the
// command line and the job macro, so this is all a macro has to declare.
//
// The class therefore has to be default-constructible: parameters arrive
// through Configure(), not through a constructor.
#define ANA_ANALYSIS(CLASS)                                                    \
    ana::AnalysisBase* ana_MakeAnalysis() { return new CLASS(); }

#endif
