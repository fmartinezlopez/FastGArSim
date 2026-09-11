 /***************************************************************************
 * MakeNtuple.C
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Flatten any FastGArSim ROOT file into vector ntuples, one flat tree per
 *   input tree, keeping the tree names. Which columns come out is not written
 *   down anywhere: it is derived from the ROOT dictionaries of whatever the
 *   file holds, so this works on simulation output, on reconstruction output
 *   with any set of modules, and on files whose products did not exist when
 *   this macro was written.
 *
 *   The output carries a Schema tree saying which column came from which
 *   branch, and a Provenance entry recording the conversion.
 *
 * Usage:
 *   MakeNtuple <input> <output> [trees] [maxDepth] [maxEntries]
 *
 *     trees       comma-separated list of trees to flatten; empty means all
 *     maxDepth    how far to follow nested collections (default 4)
 *     maxEntries  stop after this many entries of each tree (default: all)
 *
 *   For example
 *
 *     MakeNtuple reco.root ntuple.root
 *     MakeNtuple reco.root clusters.root Reco
 *     root -l 'MakeNtuple.C("sim.root", "ntuple.root")'
 *
 ***************************************************************************/

R__LOAD_LIBRARY(libGArCommon)

#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "TClass.h"
#include "TFile.h"
#include "TKey.h"
#include "TList.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "ProductSchema.hh"
#include "TreeFlattener.hh"

namespace {

std::vector<std::string> SplitList(const std::string& text)
{
    std::vector<std::string> parts;
    size_t start = 0;
    while (start <= text.size()) {
        const size_t comma = text.find(',', start);
        const size_t end = (comma == std::string::npos) ? text.size() : comma;
        std::string part = text.substr(start, end - start);

        const size_t first = part.find_first_not_of(" \t");
        const size_t last = part.find_last_not_of(" \t");
        if (first != std::string::npos) parts.push_back(part.substr(first, last - first + 1));

        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return parts;
}

// Every tree in the file, in the order they were written, skipping the
// metadata trees this tool writes itself
std::vector<std::string> TreesIn(TFile* file)
{
    std::vector<std::string> names;
    std::set<std::string> seen;

    TIter next(file->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next())) {
        const std::string name = key->GetName();
        if (!seen.insert(name).second) continue;                 // older cycle
        if (fastgarsim::ProductSchema::IsMetadataTree(name)) continue;

        const TClass* cls = TClass::GetClass(key->GetClassName());
        if (!cls || !cls->InheritsFrom(TTree::Class())) continue;

        names.push_back(name);
    }
    return names;
}

} // anonymous namespace

void MakeNtuple(const char* inputFileName,
                const char* outputFileName,
                const char* treeNames = "",
                Int_t maxDepth = 4,
                Long64_t maxEntries = -1)
{
#ifdef __CLING__
    gSystem->Load("libSimDataDict");
    gSystem->Load("libDigiDataDict");
#endif

    TFile* input = TFile::Open(inputFileName, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "Error: cannot open input file " << inputFileName << std::endl;
        return;
    }

    std::vector<std::string> trees = SplitList(treeNames ? treeNames : "");
    if (trees.empty()) trees = TreesIn(input);

    if (trees.empty()) {
        std::cerr << "Error: no trees to flatten in " << inputFileName << std::endl;
        input->Close();
        return;
    }

    TFile* output = TFile::Open(outputFileName, "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error: cannot create output file " << outputFileName << std::endl;
        input->Close();
        return;
    }
    output->SetCompressionLevel(5);

    // The input file already describes itself, when it was written by a stage
    // that knew how to; that is what says which module produced what
    const fastgarsim::ProductSchema inputSchema = fastgarsim::ProductSchema::Load(input);
    fastgarsim::ProductSchema outputSchema;

    for (const std::string& treeName : trees) {

        TTree* inputTree = dynamic_cast<TTree*>(input->Get(treeName.c_str()));
        if (!inputTree) {
            std::cerr << "Warning: no tree '" << treeName << "' in " << inputFileName
                      << "; skipping it." << std::endl;
            continue;
        }

        std::cout << "\n---- " << treeName << " ("
                  << inputTree->GetEntries() << " entries) ----" << std::endl;

        fastgarsim::FlattenOptions options;
        options.maxDepth = maxDepth;

        fastgarsim::TreeFlattener flattener(options);
        if (!flattener.Connect(inputTree)) {
            std::cerr << "Warning: nothing in '" << treeName
                      << "' could be flattened; skipping it." << std::endl;
            continue;
        }
        flattener.PrintLayout(std::cout);

        output->cd();
        TTree* outputTree = new TTree(treeName.c_str(),
                                      TString::Format("%s, flattened",
                                                      inputTree->GetTitle()).Data());
        flattener.Book(outputTree);

        Long64_t nEntries = inputTree->GetEntries();
        if (maxEntries >= 0 && maxEntries < nEntries) nEntries = maxEntries;

        Long64_t reportEvery = nEntries / 10;
        if (reportEvery <= 0) reportEvery = 1;

        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            if (entry % reportEvery == 0) {
                std::cout << "   entry " << entry << " / " << nEntries << std::endl;
            }
            if (flattener.FillEntry(entry) < 0) {
                std::cerr << "Error: could not read entry " << entry
                          << " of '" << treeName << "'; stopping here." << std::endl;
                break;
            }
            outputTree->Fill();
        }

        outputTree->OptimizeBaskets();
        output->cd();
        outputTree->Write();

        // Say where each column came from, and keep whatever the input file
        // recorded about the branch it was flattened from
        for (fastgarsim::ProductInfo product : flattener.Schema(treeName)) {
            if (const fastgarsim::ProductInfo* source =
                    inputSchema.Find(product.producer, treeName)) {
                product.parameters = source->parameters;
                if (!source->producer.empty()) {
                    // Name the module rather than the branch, when it is known
                    product.producer = source->producer + " / " + product.producer;
                    product.producerType = source->producerType;
                }
            }
            outputSchema.AddProduct(product);
        }

        inputTree->ResetBranchAddresses();
    }

    // Carry the history over, so a flat ntuple still says how it was made
    for (const fastgarsim::JobInfo& job : inputSchema.Jobs()) outputSchema.AddJob(job);

    fastgarsim::JobInfo job;
    job.stage = "ntuple";
    job.timestamp = fastgarsim::ProductSchema::Now();
    job.input = inputFileName;
    job.modules = "MakeNtuple, maxDepth " + std::to_string(maxDepth);
    outputSchema.AddJob(job);
    outputSchema.Write(output);

    std::cout << "\n=== Summary ===\n"
              << "Input:   " << inputFileName << "\n"
              << "Output:  " << outputFileName << "\n"
              << "Size:    " << output->GetSize() / 1024.0 / 1024.0 << " MB\n"
              << "Columns: " << outputSchema.Products().size() << std::endl;

    output->Close();
    input->Close();
    delete output;
    delete input;
}
