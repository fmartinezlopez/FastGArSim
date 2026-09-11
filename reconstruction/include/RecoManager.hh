//
// RecoManager.hh - Main reconstruction manager class
//

#ifndef RecoManager_h
#define RecoManager_h 1

#include <string>
#include <vector>

#include "ProductSchema.hh"

class TFile;
class TTree;
class RecoModule;
class RecoStore;
class MacroParser;
struct RecoEvent;

class RecoManager {
public:
    RecoManager();
    ~RecoManager();

    // Configuration
    bool LoadMacro(const std::string& macroFile);
    void AddModule(RecoModule* module);

    // Main reconstruction method.
    //
    // By default the output file starts life as a copy of the input, so that
    // the result holds the simulation's Events and Geometry trees as well as
    // the reconstruction's own tree and everything downstream needs one file.
    // Set /reco/global/copyInput false in the macro to write the reconstruction
    // products on their own.
    bool RunReconstruction(const std::string& inputFile, const std::string& outputFile);

private:
    // Module management
    void InitializeModules();
    void ExecuteModules();
    void FinalizeModules();
    void CreateModule(const std::string& name, const std::string& type);

    // Dependency / consistency check.
    // Verifies that every module's declared inputs are satisfied either by a
    // simulation TTree branch or by the output of an earlier module.
    // Returns false (and prints diagnostics) if any required input is missing.
    bool CheckConsistency() const;

    // I/O management
    bool OpenFiles();
    void InitializeOutput();
    void FillEvent();
    void ResetEvent();

    // Record what this pass put in the file, so that a reader can find out
    // what the file holds without being told separately
    void WriteSchema();

    // ROOT I/O
    TFile* fInputFile;
    TFile* fOutputFile;
    TTree* fInputTree;
    TTree* fOutputTree;

    // Event data
    RecoEvent* fEvent;

    // Shared data store (inter-module object passing)
    RecoStore* fStore;

    // Module management
    std::vector<RecoModule*> fModules;
    MacroParser* fMacroParser;

    // Job settings, taken from the macro's /reco/global block
    std::string fInputFileName;
    std::string fOutputFileName;
    std::string fInputTreeName;
    std::string fOutputTreeName;
    bool fCopyInput;

    // One entry per branch this pass added, with the module that added it
    std::vector<fastgarsim::ProductInfo> fProducts;

    // Verbose output
    bool fVerbose;
};

#endif
