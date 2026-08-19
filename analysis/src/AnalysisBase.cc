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
    delete fGenieTree;
    fGenieTree = nullptr;

    if (fGeoFile) {
        if (fGeoFile->IsOpen()) fGeoFile->Close();
        delete fGeoFile;
        fGeoFile = nullptr;
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

    if (!Initialize()) {
        CloseInputs();
        return kFALSE;
    }

    BeginJob();
    ProcessEvents();
    EndJob();
    Finalize();

    return kTRUE;
}

Bool_t AnalysisBase::Execute(const char* simFileSpec,
                             const char* outputFileName,
                             const char* genieFileSpec)
{
    AnalysisConfig config;
    config.simFiles = simFileSpec ? simFileSpec : "";
    config.outputFile = outputFileName ? outputFileName : "";
    config.genieFiles = genieFileSpec ? genieFileSpec : "";
    return Execute(config);
}

/* -------------------------------------------------------------------------- */
/*                                Initialization                              */
/* -------------------------------------------------------------------------- */

Bool_t AnalysisBase::Initialize()
{
    std::cout << "FastGArSim analysis\n"
              << "    simulation input: " << fConfig.simFiles << "\n"
              << "    GENIE input:      "
              << (fConfig.genieFiles.empty() ? "(none)" : fConfig.genieFiles) << "\n"
              << "    output:           "
              << (fConfig.outputFile.empty() ? "(none)" : fConfig.outputFile) << "\n"
              << std::endl;

    if (fConfig.simFiles.empty()) {
        std::cerr << "Error: no simulation input given" << std::endl;
        return kFALSE;
    }

    // Resolve wildcards and rewrite /pnfs paths as XRootD URLs
    fSimFileNames = ExpandInputFiles(fConfig.simFiles, fConfig.xrootdForPnfs);
    if (fSimFileNames.empty()) return kFALSE;

    std::cout << "Chaining " << fSimFileNames.size() << " simulation file(s)"
              << std::endl;

    // Chain the FastGArSim analysis trees
    fSimTree = new TChain(fConfig.simTreeName.c_str());
    for (const std::string& file : fSimFileNames) fSimTree->Add(file.c_str());

    // Force the first file open, so that the branches can be inspected and
    // a missing or empty tree is reported here rather than mid-loop
    if (fSimTree->LoadTree(0) < 0) {
        std::cerr << "Error: Could not read TTree " << fConfig.simTreeName
                  << " from the input files." << std::endl;
        ExplainOpenFailure(fSimFileNames.front());
        return kFALSE;
    }
    sim.Connect(fSimTree);

    // Get the FastGArSim geometry TTree from the first input file. Missing
    // geometry is not fatal: the analysis simply sees the default-constructed
    // GeometryInfo.
    fGeoFile = TFile::Open(fSimFileNames.front().c_str(), "READ");
    if (!fGeoFile || fGeoFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << fSimFileNames.front() << std::endl;
        ExplainOpenFailure(fSimFileNames.front());
        return kFALSE;
    }

    fGeoTree = dynamic_cast<TTree*>(fGeoFile->Get(fConfig.geoTreeName.c_str()));
    if (fGeoTree) {
        geo.Connect(fGeoTree);
        if (fGeoTree->GetEntries() > 0) fGeoTree->GetEntry(0);
        if (fConfig.printGeometry) {
            if (fSimFileNames.size() > 1) {
                std::cout << "Geometry taken from " << fSimFileNames.front() << std::endl;
            }
            geo.Print();
        }
    } else {
        std::cout << "Warning: Could not find TTree " << fConfig.geoTreeName
                  << " in " << fSimFileNames.front()
                  << "; geometry parameters keep their default values." << std::endl;
    }

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
        if (fSimFileNames.size() > 1 || genieFileNames.size() > 1) {
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

    const Long64_t nEntries = fSimTree->GetEntries();

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
        fSimTree->GetEntry(iEvent);
        sim.Update();

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
