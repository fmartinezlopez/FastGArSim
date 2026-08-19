//
// RecoManager.hh - Main reconstruction manager class
//

#ifndef RecoManager_h
#define RecoManager_h 1

#include <string>
#include <vector>

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

    // Main reconstruction method
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
    bool OpenInputFile(const std::string& inputFile);
    void CreateOutputFile(const std::string& outputFile);
    void InitializeOutput();
    void FillEvent();
    void ResetEvent();

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

    // Verbose output
    bool fVerbose;
};

#endif
