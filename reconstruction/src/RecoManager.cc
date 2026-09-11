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

#include <sys/stat.h>

#include <iostream>
#include <set>

namespace {

// Whether two paths name the same file on disk. False when either does not
// exist, which is the usual case for the output.
bool SameFile(const std::string& first, const std::string& second)
{
    struct stat a;
    struct stat b;
    if (::stat(first.c_str(), &a) != 0) return false;
    if (::stat(second.c_str(), &b) != 0) return false;
    return a.st_dev == b.st_dev && a.st_ino == b.st_ino;
}

} // anonymous namespace

RecoManager::RecoManager()
    : fInputFile(nullptr), fOutputFile(nullptr),
      fInputTree(nullptr), fOutputTree(nullptr),
      fEvent(nullptr), fStore(nullptr),
      fMacroParser(nullptr),
      fInputTreeName("Events"), fOutputTreeName("Reco"), fCopyInput(true),
      fVerbose(true)
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

    // When the input was copied to make the output, both point at one file
    if (fInputFile && fInputFile != fOutputFile) {
        fInputFile->Close();
        delete fInputFile;
    }

    if (fOutputFile) {
        fOutputFile->Close();
        delete fOutputFile;
    }
    fInputFile = nullptr;
    fOutputFile = nullptr;

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

    // The factory hands back a module carrying its class default name. The
    // instance name from the macro is the one that identifies it in messages
    // and in the file's schema record, so it wins.
    module->SetName(name);
    module->SetType(type);

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

bool RecoManager::OpenFiles()
{
    if (fCopyInput) {

        // A copy is only a copy if it is a different file. Comparing the
        // strings would miss two spellings of one path, so the file system is
        // asked instead: same device and inode means same file, however it
        // was written or symlinked.
        if (SameFile(fInputFileName, fOutputFileName)) {
            std::cerr << "Error: the output file is the input file. Either give a "
                         "different\n       output name or put "
                         "/reco/global/copyInput false in the macro." << std::endl;
            return false;
        }

        std::cout << "\n Copying " << fInputFileName << " to " << fOutputFileName
                  << ",\n so that the result holds the simulation trees as well "
                     "as the reconstruction..." << std::endl;

        if (!TFile::Cp(fInputFileName.c_str(), fOutputFileName.c_str(), kFALSE)) {
            std::cerr << "Error: could not copy the input file to "
                      << fOutputFileName << std::endl;
            return false;
        }

        fOutputFile = TFile::Open(fOutputFileName.c_str(), "UPDATE");
        if (!fOutputFile || fOutputFile->IsZombie()) {
            std::cerr << "Error: could not reopen the copy: " << fOutputFileName << std::endl;
            return false;
        }

        // Read from the copy and write to it: one open file, and the trees
        // cannot drift apart
        fInputFile = fOutputFile;

        // Re-running reconstruction on a file that already has a Reco tree
        // would leave two cycles of it behind, the older one still readable
        if (fOutputFile->Get(fOutputTreeName.c_str())) {
            std::cout << " Replacing the '" << fOutputTreeName
                      << "' tree already in the file" << std::endl;
            fOutputFile->Delete((fOutputTreeName + ";*").c_str());
        }

    } else {

        fInputFile = TFile::Open(fInputFileName.c_str(), "READ");
        if (!fInputFile || fInputFile->IsZombie()) {
            std::cerr << "Error: Cannot open input file: " << fInputFileName << std::endl;
            return false;
        }

        fOutputFile = new TFile(fOutputFileName.c_str(), "RECREATE");
        if (!fOutputFile || fOutputFile->IsZombie()) {
            std::cerr << "Error: Cannot create output file: " << fOutputFileName << std::endl;
            return false;
        }
    }

    fInputTree = dynamic_cast<TTree*>(fInputFile->Get(fInputTreeName.c_str()));
    if (!fInputTree) {
        std::cerr << "Error: Cannot find tree '" << fInputTreeName << "' in "
                  << fInputFileName << std::endl;
        return false;
    }

    return true;
}

void RecoManager::InitializeOutput()
{
    fOutputFile->cd();
    fOutputTree = new TTree(fOutputTreeName.c_str(), "Reconstruction Output");
    fOutputTree->Branch("RecoEvent", &fEvent);

    fastgarsim::ProductInfo product;
    product.tree = fOutputTreeName;
    product.branch = "RecoEvent";
    product.type = "RecoEvent";
    product.producer = "RecoManager";
    product.producerType = "RecoManager";
    fProducts.push_back(product);
}

void RecoManager::ResetEvent()
{
    *fEvent = RecoEvent();
}

void RecoManager::FillEvent()
{
    fOutputTree->Fill();
}

void RecoManager::WriteSchema()
{
    // Start from what the file actually holds now, so that the record cannot
    // disagree with the file, and label it with what this pass knows
    fastgarsim::ProductSchema schema = fastgarsim::ProductSchema::Describe(fOutputFile);

    // Whatever an earlier pass recorded about its own products, kept
    const fastgarsim::ProductSchema previous = fastgarsim::ProductSchema::Read(fOutputFile);
    schema.Annotate(previous);

    fastgarsim::ProductSchema mine;
    for (const fastgarsim::ProductInfo& product : fProducts) mine.AddProduct(product);
    schema.Annotate(mine);

    for (const fastgarsim::JobInfo& job : previous.Jobs()) schema.AddJob(job);

    fastgarsim::JobInfo job;
    job.stage = "reconstruction";
    job.timestamp = fastgarsim::ProductSchema::Now();
    job.input = fInputFileName;
    if (fMacroParser) {
        job.macro = fMacroParser->GetPath();
        job.macroText = fMacroParser->GetText();
    }
    for (const RecoModule* module : fModules) {
        if (!module->IsEnabled()) continue;
        if (!job.modules.empty()) job.modules += ", ";
        job.modules += module->GetName() + ":" + module->GetType();
    }
    schema.AddJob(job);

    schema.Write(fOutputFile);

    std::cout << "\n Products in " << fOutputFileName << ":" << std::endl;
    schema.Print(std::cout);
}

void RecoManager::InitializeModules()
{
    std::cout << "\n Initializing reconstruction modules..." << std::endl;

    // A module books its own branches in Initialize(), and does not report
    // which. Watching the output tree either side of the call attributes them
    // without every module having to declare them a second time.
    auto branchNames = [this]() {
        std::set<std::string> names;
        TObjArray* branches = fOutputTree->GetListOfBranches();
        for (int i = 0; i < branches->GetEntriesFast(); ++i) {
            names.insert(branches->At(i)->GetName());
        }
        return names;
    };

    for (auto module : fModules) {
        if (!module->IsEnabled()) continue;

        module->SetInputFile(fInputFile);
        module->SetInputTree(fInputTree);
        module->SetOutputTree(fOutputTree);
        module->SetEvent(fEvent);
        module->SetStore(fStore);

        const std::set<std::string> before = branchNames();
        module->Initialize();

        // Flatten the module's configuration into one line for the record
        std::string parameters;
        for (const auto& parameter : module->GetParameters()) {
            if (!parameters.empty()) parameters += "; ";
            parameters += parameter.first + "=" + parameter.second;
        }

        for (const std::string& name : branchNames()) {
            if (before.count(name)) continue;

            fastgarsim::ProductInfo product;
            product.tree = fOutputTreeName;
            product.branch = name;
            product.type = fastgarsim::ProductSchema::BranchType(fOutputTree, name.c_str());
            product.producer = module->GetName();
            product.producerType = module->GetType();
            product.parameters = parameters;
            fProducts.push_back(product);
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

    fInputFileName = inputFile;
    fOutputFileName = outputFile;

    if (fMacroParser) {
        fInputTreeName  = fMacroParser->GetGlobalParameter("inputTreeName", "Events");
        fOutputTreeName = fMacroParser->GetGlobalParameter("outputTreeName", "Reco");
        fCopyInput      = fMacroParser->GetGlobalParameterBool("copyInput", true);
    }

    if (!OpenFiles()) {
        return false;
    }

    Long64_t nEntries = fInputTree->GetEntries();
    std::cout << " Found " << nEntries << " events to process" << std::endl;

    // Verify dependency graph before touching events
    if (!CheckConsistency()) {
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
    fOutputTree->Write(nullptr, TObject::kOverwrite);

    // Last, so that it describes the finished file
    WriteSchema();

    fOutputFile->Close();

    std::cout << "\n Processed " << nEntries << " events" << std::endl;

    return true;
}
