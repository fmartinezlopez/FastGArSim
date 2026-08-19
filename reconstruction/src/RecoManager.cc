//
// RecoManager.cc - Implementation of reconstruction manager
//

#include "RecoManager.hh"
#include "RecoDataTypes.hh"
#include "RecoModule.hh"
#include "RecoStore.hh"
#include "ModuleFactory.hh"
#include "MacroParser.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <iostream>
#include <set>

RecoManager::RecoManager()
    : fInputFile(nullptr), fOutputFile(nullptr),
      fInputTree(nullptr), fOutputTree(nullptr),
      fEvent(nullptr), fStore(nullptr),
      fMacroParser(nullptr), fVerbose(true)
{
    fEvent       = new RecoEvent();
    fStore       = new RecoStore();
    fMacroParser = new MacroParser();
}

RecoManager::~RecoManager()
{
    for (auto module : fModules) {
        delete module;
    }
    fModules.clear();

    if (fInputFile) {
        fInputFile->Close();
        delete fInputFile;
    }

    if (fOutputFile) {
        fOutputFile->Close();
        delete fOutputFile;
    }

    delete fEvent;
    delete fStore;
    delete fMacroParser;
}

bool RecoManager::LoadMacro(const std::string& macroFile)
{
    if (!fMacroParser->ParseFile(macroFile)) {
        return false;
    }

    const auto& configs = fMacroParser->GetModuleConfigs();
    for (const auto& config : configs) {
        if (config.enabled) {
            CreateModule(config.name, config.type);

            RecoModule* module = fModules.back();
            for (const auto& param : config.parameters) {
                module->SetParameter(param.first, param.second);
            }
        }
    }

    return true;
}

void RecoManager::AddModule(RecoModule* module)
{
    if (module) {
        fModules.push_back(module);
    }
}

void RecoManager::CreateModule(const std::string& name, const std::string& type)
{
    RecoModule* module = ModuleFactory::Instance().Create(type);

    if (!module) {
        std::cerr << "Error: Unknown module type '" << type << "'." << std::endl;
        std::cerr << "       Registered types:";
        for (const auto& t : ModuleFactory::Instance().GetRegisteredTypes()) {
            std::cerr << " " << t;
        }
        std::cerr << std::endl;
        return;
    }

    // Override the module name with the one given in the macro
    // (the factory-created module carries the class default name; we want
    //  the instance name so parameters are routed correctly)
    module->SetParameter("__instanceName__", name);

    AddModule(module);
    std::cout << "   Created module: " << name << " (type: " << type << ")" << std::endl;
}

bool RecoManager::CheckConsistency() const
{
    std::cout << "\n Checking module dependency consistency..." << std::endl;

    // Seed the available-object set with every branch in the input TTree
    std::set<std::string> available;
    if (fInputTree) {
        TObjArray* branches = fInputTree->GetListOfBranches();
        for (int i = 0; i < branches->GetEntries(); ++i) {
            available.insert(branches->At(i)->GetName());
        }
        std::cout << "   Simulation branches available: " << available.size() << std::endl;
    } else {
        std::cout << "   Warning: No input tree open yet; "
                     "simulation branch names are not checked." << std::endl;
    }

    bool ok = true;
    for (const auto* module : fModules) {
        if (!module->IsEnabled()) continue;

        const std::string& mName = module->GetName();

        // Check that all required inputs are already available
        for (const auto& spec : module->GetInputSpec()) {
            if (available.count(spec.key) == 0) {
                if (spec.optional) {
                    std::cout << "   Warning: Module '" << mName
                              << "' has optional input '" << spec.key
                              << "' (" << spec.typeName
                              << ") that is not (yet) available." << std::endl;
                } else {
                    std::cerr << "   Error:   Module '" << mName
                              << "' requires '" << spec.key
                              << "' (" << spec.typeName
                              << ") but it is not available at this point." << std::endl;
                    ok = false;
                }
            }
        }

        // After this module runs, its outputs become available
        for (const auto& spec : module->GetOutputSpec()) {
            if (available.count(spec.key)) {
                std::cerr << "   Error:   Module '" << mName
                          << "' declares output '" << spec.key
                          << "' which is already produced by an earlier module." << std::endl;
                ok = false;
            }
            available.insert(spec.key);
        }

        if (ok) {
            std::cout << "   [OK] " << mName;
            const auto& ins  = module->GetInputSpec();
            const auto& outs = module->GetOutputSpec();
            if (!ins.empty() || !outs.empty()) {
                std::cout << "  (reads:";
                for (const auto& s : ins)  std::cout << " " << s.key;
                std::cout << "  writes:";
                for (const auto& s : outs) std::cout << " " << s.key;
                std::cout << ")";
            }
            std::cout << std::endl;
        }
    }

    if (!ok) {
        std::cerr << "\n Consistency check FAILED. "
                     "Fix the module order or missing dependencies in your macro." << std::endl;
    } else {
        std::cout << " Consistency check passed." << std::endl;
    }

    return ok;
}

bool RecoManager::OpenInputFile(const std::string& inputFile)
{
    fInputFile = TFile::Open(inputFile.c_str(), "READ");
    if (!fInputFile || fInputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file: " << inputFile << std::endl;
        return false;
    }

    // Tree name is configurable via /reco/global/inputTreeName in the macro
    std::string treeName = fMacroParser->GetGlobalParameter("inputTreeName", "Events");
    fInputTree = (TTree*)fInputFile->Get(treeName.c_str());
    if (!fInputTree) {
        std::cerr << "Error: Cannot find tree '" << treeName << "' in input file" << std::endl;
        return false;
    }

    return true;
}

void RecoManager::CreateOutputFile(const std::string& outputFile)
{
    fOutputFile = new TFile(outputFile.c_str(), "RECREATE");
}

void RecoManager::InitializeOutput()
{
    fOutputTree = new TTree("RecoTree", "Reconstruction Output");
    fOutputTree->Branch("Event", &fEvent);
}

void RecoManager::ResetEvent()
{
    *fEvent = RecoEvent();
}

void RecoManager::FillEvent()
{
    fOutputTree->Fill();
}

void RecoManager::InitializeModules()
{
    std::cout << "\n Initializing reconstruction modules..." << std::endl;

    for (auto module : fModules) {
        if (module->IsEnabled()) {
            module->SetInputFile(fInputFile);
            module->SetInputTree(fInputTree);
            module->SetOutputTree(fOutputTree);
            module->SetEvent(fEvent);
            module->SetStore(fStore);
            module->Initialize();
        }
    }
}

void RecoManager::ExecuteModules()
{
    for (auto module : fModules) {
        if (module->IsEnabled()) {
            module->Execute();
        }
    }
}

void RecoManager::FinalizeModules()
{
    std::cout << "\n Finalizing reconstruction modules..." << std::endl;

    for (auto module : fModules) {
        if (module->IsEnabled()) {
            module->Finalize();
        }
    }
}

bool RecoManager::RunReconstruction(const std::string& inputFile, const std::string& outputFile)
{
    if (fModules.empty()) {
        std::cerr << "Warning: No reconstruction modules configured!" << std::endl;
        std::cerr << "         Use LoadMacro() or AddModule() to configure reconstruction" << std::endl;
        return false;
    }

    if (!OpenInputFile(inputFile)) {
        return false;
    }

    Long64_t nEntries = fInputTree->GetEntries();
    std::cout << " Found " << nEntries << " events to process" << std::endl;

    // Verify dependency graph before touching events
    if (!CheckConsistency()) {
        return false;
    }

    CreateOutputFile(outputFile);
    if (!fOutputFile || fOutputFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
        return false;
    }

    InitializeOutput();
    InitializeModules();

    std::cout << "\n Processing events..." << std::endl;
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        if (iEntry % 100 == 0 && fVerbose) {
            std::cout << " Processing event " << iEntry << " / " << nEntries << std::endl;
        }

        fInputTree->GetEntry(iEntry);
        ResetEvent();
        fEvent->eventID = iEntry;

        ExecuteModules();
        FillEvent();
    }

    FinalizeModules();

    fOutputFile->cd();
    fOutputTree->Write();
    fOutputFile->Close();

    std::cout << "\n Processed " << nEntries << " events" << std::endl;

    return true;
}
