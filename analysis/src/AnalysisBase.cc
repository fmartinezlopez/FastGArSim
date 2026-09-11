 /***************************************************************************
 * AnalysisBase.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the analysis driver: input/output handling, branch
 *   wiring and the event loop.
 *
 ***************************************************************************/

#include "AnalysisBase.hh"

#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "AnalysisInput.hh"

namespace ana {

AnalysisBase::AnalysisBase() = default;

AnalysisBase::~AnalysisBase()
{
    CloseInputs();

    if (fOutputFile) {
        if (fOutputFile->IsOpen()) fOutputFile->Close();
        delete fOutputFile;
        fOutputFile = nullptr;
    }
}

void AnalysisBase::CloseInputs()
{
    // Deleting a chain closes whichever of its files is currently open
    delete fSimTree;
    fSimTree = nullptr;
    delete fRecoTree;
    fRecoTree = nullptr;
    delete fGenieTree;
    fGenieTree = nullptr;

    if (fFirstFile) {
        if (fFirstFile->IsOpen()) fFirstFile->Close();
        delete fFirstFile;
        fFirstFile = nullptr;
    }
    fGeoTree = nullptr;
}

/* -------------------------------------------------------------------------- */
/*                                   Driver                                   */
/* -------------------------------------------------------------------------- */

Bool_t AnalysisBase::Execute(const AnalysisConfig& config)
{
    fConfig = config;
    fAbort = kFALSE;

    // The analysis reads its parameters first, before anything is opened, so
    // that a job configured with a value it cannot use costs nothing
    Configure(fConfig.params);
    if (!CheckParameters()) return kFALSE;

    if (!Initialize()) {
        CloseInputs();
        return kFALSE;
    }

    // BeginJob() is where an analysis asks for the reconstruction products it
    // needs, so the check that it got them comes straight after
    BeginJob();
    const Bool_t haveProducts = CheckRequestedProducts();
    if (!haveProducts) Abort();

    ProcessEvents();
    EndJob();
    Finalize();

    // EndJob() still ran, so any summary the analysis wanted was printed, but
    // a job that could not read what it asked for did not succeed
    return haveProducts;
}

Bool_t AnalysisBase::Execute(const char* inputFileSpec,
                             const char* outputFileName,
                             const char* genieFileSpec)
{
    AnalysisConfig config;
    config.inputFiles = inputFileSpec ? inputFileSpec : "";
    config.outputFile = outputFileName ? outputFileName : "";
    config.genieFiles = genieFileSpec ? genieFileSpec : "";
    return Execute(config);
}

/* -------------------------------------------------------------------------- */
/*                                Initialization                              */
/* -------------------------------------------------------------------------- */

//---------------------------------------------------------------------------
// Report what Configure() made of the parameters it was given. Nothing checks
// parameter names at compile time, so this is where both kinds of mistake in
// a job macro surface: a value of the wrong type, which stops the job, and a
// name no analysis ever asked for, which is almost always a misspelling and
// is reported but left to the user to judge -- the same macro may be shared
// between several analyses.
//---------------------------------------------------------------------------
Bool_t AnalysisBase::CheckParameters()
{
    const std::vector<std::string> unused = fConfig.params.UnusedKeys();
    if (!unused.empty()) {
        std::cout << "\nWarning: " << unused.size() << " parameter"
                  << (unused.size() == 1 ? "" : "s") << " the analysis never asked for";
        for (const std::string& name : unused) std::cout << "\n   /ana/" << name;
        std::cout << "\n   (ignored -- check the spelling against what the analysis reads)\n"
                  << std::endl;
    }

    if (fConfig.params.NErrors() > 0) {
        std::cerr << "\nError: " << fConfig.params.NErrors() << " parameter"
                  << (fConfig.params.NErrors() == 1 ? " was" : "s were")
                  << " set to a value that could not be read; nothing was run.\n"
                  << std::endl;
        return kFALSE;
    }

    // Configure() can also turn a job down itself, for a combination of
    // parameters that is readable but makes no sense
    if (fAbort) {
        std::cerr << "\nError: the analysis rejected its configuration; nothing was run.\n"
                  << std::endl;
        return kFALSE;
    }
    return kTRUE;
}

Bool_t AnalysisBase::Initialize()
{
    std::cout << "FastGArSim analysis\n"
              << "    input:       " << fConfig.inputFiles << "\n"
              << "    GENIE input: "
              << (fConfig.genieFiles.empty() ? "(none)" : fConfig.genieFiles) << "\n"
              << "    output:      "
              << (fConfig.outputFile.empty() ? "(none)" : fConfig.outputFile) << "\n"
              << std::endl;

    if (fConfig.inputFiles.empty()) {
        std::cerr << "Error: no input given" << std::endl;
        return kFALSE;
    }

    if (!OpenInputs()) return kFALSE;

    // Chain the GENIE files, if any were requested
    if (!fConfig.genieFiles.empty()) {
        const std::vector<std::string> genieFileNames =
            ExpandInputFiles(fConfig.genieFiles, fConfig.xrootdForPnfs);
        if (genieFileNames.empty()) return kFALSE;

        fGenieTree = new TChain(fConfig.genieTreeName.c_str());
        for (const std::string& file : genieFileNames) fGenieTree->Add(file.c_str());

        if (fGenieTree->LoadTree(0) < 0) {
            std::cerr << "Error: Could not read TTree " << fConfig.genieTreeName
                      << " from the GENIE files." << std::endl;
            ExplainOpenFailure(genieFileNames.front());
            return kFALSE;
        }
        genie.Connect(fGenieTree);

        // The GENIE entry is looked up by eventID, which is only meaningful
        // across a chain if the simulation numbered its events globally
        if (fInputFileNames.size() > 1 || genieFileNames.size() > 1) {
            std::cout << "Warning: the GENIE record is looked up by eventID, so with "
                      << "several chained files the eventID must be unique across the "
                      << "whole chain. Check this before trusting the truth-level "
                      << "quantities." << std::endl;
        }
    }

    // Create the output file and tree. Doing this last leaves the output file
    // as the current directory, so histograms booked in BeginJob() are
    // attached to it and written out automatically.
    if (!fConfig.outputFile.empty()) {
        fOutputFile = new TFile(fConfig.outputFile.c_str(), "RECREATE");
        if (fOutputFile->IsZombie()) {
            std::cerr << "Error: Could not create output file "
                      << fConfig.outputFile << std::endl;
            return kFALSE;
        }
        fOutputFile->cd();
        fOutputTree = new TTree(fConfig.outputTreeName.c_str(),
                                fConfig.outputTreeName.c_str());
    }

    return kTRUE;
}

Bool_t AnalysisBase::OpenInputs()
{
    // Resolve wildcards and rewrite /pnfs paths as XRootD URLs
    fInputFileNames = ExpandInputFiles(fConfig.inputFiles, fConfig.xrootdForPnfs);
    if (fInputFileNames.empty()) return kFALSE;

    std::cout << "Chaining " << fInputFileNames.size() << " input file(s)" << std::endl;

    const std::string& first = fInputFileNames.front();

    // The first file is opened directly as well as through the chains, for
    // the geometry and the schema record
    fFirstFile = TFile::Open(first.c_str(), "READ");
    if (!fFirstFile || fFirstFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << first << std::endl;
        ExplainOpenFailure(first);
        return kFALSE;
    }

    // What this file says it holds. Files written before the schema record
    // existed are described by walking their trees instead.
    fSchema = fastgarsim::ProductSchema::Load(fFirstFile);

    /* ----------------------------- Simulation ----------------------------- */

    const Bool_t hasSim = fFirstFile->Get(fConfig.simTreeName.c_str()) != nullptr;

    if (!hasSim) {
        // The one wrong guess worth naming: a flat ntuple from the old
        // converter, which this framework no longer reads
        if (fFirstFile->Get("AnaTree")) {
            std::cerr << "Error: " << first << " holds an 'AnaTree', which is the flat "
                         "ntuple format\n       the analysis framework used to read. It "
                         "now reads the simulation and\n       reconstruction objects "
                         "directly, so give it the simulation or\n       reconstruction "
                         "file instead. MakeNtuple still produces flat ntuples,\n"
                         "       for reading outside this framework." << std::endl;
            return kFALSE;
        }
        if (fConfig.requireSim) {
            std::cerr << "Error: no '" << fConfig.simTreeName << "' tree in " << first
                      << ".\n       Set AnalysisConfig::requireSim false to run on "
                         "reconstruction products alone." << std::endl;
            return kFALSE;
        }
        std::cout << "Note: no '" << fConfig.simTreeName << "' tree in the input; "
                  << "the simulation record will be empty." << std::endl;
    } else {
        fSimTree = new TChain(fConfig.simTreeName.c_str());
        for (const std::string& file : fInputFileNames) fSimTree->Add(file.c_str());

        // Force the first file open, so that a missing or empty tree is
        // reported here rather than mid-loop
        if (fSimTree->LoadTree(0) < 0) {
            std::cerr << "Error: Could not read TTree " << fConfig.simTreeName
                      << " from the input files." << std::endl;
            ExplainOpenFailure(first);
            return kFALSE;
        }
        sim.Connect(fSimTree, fConfig.simBranchName.c_str());
    }

    /* --------------------------- Reconstruction --------------------------- */

    if (fFirstFile->Get(fConfig.recoTreeName.c_str())) {
        fRecoTree = new TChain(fConfig.recoTreeName.c_str());
        for (const std::string& file : fInputFileNames) fRecoTree->Add(file.c_str());

        if (fRecoTree->LoadTree(0) < 0) {
            std::cerr << "Warning: could not read TTree " << fConfig.recoTreeName
                      << " from the input files; carrying on without "
                         "reconstruction products." << std::endl;
            delete fRecoTree;
            fRecoTree = nullptr;
        }
    }

    if (fRecoTree) {
        // The two chains are stepped together by entry number, which only
        // means anything if they have the same number of entries. They do
        // when the reconstruction wrote its tree into a copy of its input,
        // which is what it does by default.
        if (fSimTree && fSimTree->GetEntries() != fRecoTree->GetEntries()) {
            std::cerr << "Error: '" << fConfig.simTreeName << "' has "
                      << fSimTree->GetEntries() << " entries but '"
                      << fConfig.recoTreeName << "' has " << fRecoTree->GetEntries()
                      << ".\n       They are read entry by entry together, so they "
                         "have to match. This\n       usually means simulation and "
                         "reconstruction files have been mixed in\n       one job."
                      << std::endl;
            return kFALSE;
        }

        reco.Connect(fRecoTree, &fSchema);
        fRecoTreeNumber = fRecoTree->GetTreeNumber();

        if (fConfig.printProducts) {
            reco.Print(std::cout);
            std::cout << std::endl;
        }
    } else if (fConfig.printProducts) {
        std::cout << "No '" << fConfig.recoTreeName << "' tree in the input: "
                  << "reconstruction products are not available." << std::endl;
    }

    /* ------------------------------ Geometry ------------------------------ */

    // Missing geometry is not fatal: the analysis simply sees the
    // default-constructed GeometryInfo.
    fGeoTree = dynamic_cast<TTree*>(fFirstFile->Get(fConfig.geoTreeName.c_str()));
    if (fGeoTree) {
        geo.Connect(fGeoTree);
        if (fGeoTree->GetEntries() > 0) fGeoTree->GetEntry(0);
        if (fConfig.printGeometry) {
            if (fInputFileNames.size() > 1) {
                std::cout << "Geometry taken from " << first << std::endl;
            }
            geo.Print();
        }
    } else {
        std::cout << "Warning: Could not find TTree " << fConfig.geoTreeName
                  << " in " << first
                  << "; geometry parameters keep their default values." << std::endl;
    }

    return kTRUE;
}

Bool_t AnalysisBase::CheckRequestedProducts()
{
    if (reco.Errors().empty()) return kTRUE;

    std::cerr << "\nThis analysis cannot run on these files:" << std::endl;
    for (const std::string& error : reco.Errors()) {
        std::cerr << "    " << error << std::endl;
    }

    // The list of products has already been printed at start-up unless it was
    // switched off, in which case this is the first chance to see it
    if (!fConfig.printProducts) {
        std::cerr << "\n";
        reco.Print(std::cerr);
    }

    std::cerr << "\nRun DumpSchema on an input file to see how it was produced."
              << std::endl;

    return kFALSE;
}

/* -------------------------------------------------------------------------- */
/*                                 Event loop                                 */
/* -------------------------------------------------------------------------- */

void AnalysisBase::ProcessEvents()
{
    // BeginJob() may have aborted the job, e.g. because the input lacks
    // something the analysis needs
    if (fAbort) {
        std::cout << "Event loop skipped: the analysis aborted during BeginJob()"
                  << std::endl;
        return;
    }

    // Either chain can be the one driving the loop: an analysis may run on
    // reconstruction products with no simulation record alongside them
    TChain* driver = fSimTree ? fSimTree : fRecoTree;
    if (!driver) {
        std::cerr << "Error: nothing to loop over -- the input holds neither a '"
                  << fConfig.simTreeName << "' nor a '" << fConfig.recoTreeName
                  << "' tree." << std::endl;
        return;
    }

    const Long64_t nEntries = driver->GetEntries();

    Long64_t first = fConfig.firstEvent > 0 ? fConfig.firstEvent : 0;
    if (first > nEntries) first = nEntries;

    Long64_t last = nEntries;
    if (fConfig.maxEvents >= 0 && first + fConfig.maxEvents < last) {
        last = first + fConfig.maxEvents;
    }

    fNEvents = last - first;

    std::cout << "Number of events in file: " << nEntries << "\n"
              << "Number of events to process: " << fNEvents << std::endl;

    // Progress counter
    Long64_t reportEvery = (fConfig.nReports > 0) ? fNEvents / fConfig.nReports : 0;
    if (reportEvery <= 0) reportEvery = 1;

    // Main event loop
    for (Long64_t iEvent = first; iEvent < last; iEvent++) {

        fCurrentEntry = iEvent;

        // Print progress
        if ((iEvent - first) % reportEvery == 0) {
            std::cout << "Processing event " << iEvent << " ("
                      << (100.0 * (iEvent - first) / (fNEvents > 0 ? fNEvents : 1))
                      << "%)" << std::endl;
        }

        // Load current event
        if (fSimTree) {
            fSimTree->GetEntry(iEvent);
            sim.Update();
        }

        if (fRecoTree) {
            fRecoTree->GetEntry(iEvent);

            // A chain moving on to another file may be moving on to one with
            // a different set of products; the store makes the handles for
            // anything missing report themselves invalid rather than hand
            // back the previous file's contents
            if (fRecoTree->GetTreeNumber() != fRecoTreeNumber) {
                fRecoTreeNumber = fRecoTree->GetTreeNumber();
                reco.OnNewTree(fRecoTree->GetTree());
            }
        }

        // Load the corresponding entry in the GENIE tree
        if (fGenieTree) {
            if (sim.eventID >= 0 && sim.eventID < fGenieTree->GetEntries()) {
                fGenieTree->GetEntry(sim.eventID);
            } else {
                std::cerr << "Warning: eventID " << sim.eventID
                          << " is out of range for tree " << fConfig.genieTreeName
                          << " (" << fGenieTree->GetEntries() << " entries); "
                          << "skipping event " << iEvent << std::endl;
                continue;
            }
        }

        // Hand over to the concrete analysis
        Run();

        if (fAbort) {
            std::cout << "Event loop aborted at event " << iEvent << std::endl;
            break;
        }

    } // end loop over events

    std::cout << "Event loop finished" << std::endl;
}

/* -------------------------------------------------------------------------- */
/*                                  Output                                    */
/* -------------------------------------------------------------------------- */

void AnalysisBase::Fill()
{
    if (fOutputTree) fOutputTree->Fill();
}

void AnalysisBase::Finalize()
{
    if (fOutputFile && fOutputFile->IsOpen()) {
        fOutputFile->cd();
        // Writes the output tree together with every histogram attached to
        // the file since it was opened
        fOutputFile->Write();
        std::cout << "Wrote " << fConfig.outputFile;
        if (fOutputTree) std::cout << " (" << fOutputTree->GetEntries() << " entries)";
        std::cout << std::endl;
        fOutputFile->Close();
    }

    // Closing the file already deleted the tree and any attached histograms
    delete fOutputFile;
    fOutputFile = nullptr;
    fOutputTree = nullptr;

    CloseInputs();
}

} // namespace ana
