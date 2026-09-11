//
// RecoModule.hh - Base class for reconstruction modules
//

#ifndef RecoModule_h
#define RecoModule_h 1

#include <string>
#include <map>
#include <vector>

class TFile;
class TTree;
class RecoStore;
struct RecoEvent;

// Describes one input or output object a module declares.
// key      - Name used in RecoStore (for reco objects) or TTree branch (for sim data)
// typeName - Human-readable type, only used in diagnostics
// optional - If true, missing input is a warning rather than an error
struct ObjectSpec {
    std::string key;
    std::string typeName;
    bool optional = false;
};

class RecoModule {
public:
    RecoModule(const std::string& name);
    virtual ~RecoModule();

    // Pure virtual methods that must be implemented by derived classes
    virtual void Initialize() = 0;
    virtual void Execute() = 0;
    virtual void Finalize() = 0;

    // Declare what this module reads and produces.
    // Inputs may come from the simulation TTree branches or from the RecoStore
    // (i.e. produced by an earlier module). Outputs go into the RecoStore.
    // The base-class implementations return empty vectors (no dependencies).
    virtual std::vector<ObjectSpec> GetInputSpec()  const { return {}; }
    virtual std::vector<ObjectSpec> GetOutputSpec() const { return {}; }

    // Module configuration
    void SetParameter(const std::string& key, const std::string& value);
    std::string GetParameter(const std::string& key, const std::string& defaultValue = "") const;
    int         GetParameterInt(const std::string& key, int defaultValue = 0) const;
    double      GetParameterDouble(const std::string& key, double defaultValue = 0.0) const;
    bool        GetParameterBool(const std::string& key, bool defaultValue = false) const;

    // Every parameter set on this module, for the file's schema record
    const std::map<std::string, std::string>& GetParameters() const { return fParameters; }

    // Identity. The name is the instance name from the macro and the type is
    // the class the factory was asked for; the manager sets both, so a module
    // constructor only has to supply a sensible default name.
    std::string GetName() const { return fName; }
    void SetName(const std::string& name) { fName = name; }
    std::string GetType() const { return fType; }
    void SetType(const std::string& type) { fType = type; }

    bool IsEnabled() const { return fEnabled; }
    void SetEnabled(bool enabled) { fEnabled = enabled; }

    // Access to I/O infrastructure
    void SetInputFile(TFile* file)   { fInputFile  = file;  }
    void SetInputTree(TTree* tree)   { fInputTree  = tree;  }
    void SetOutputTree(TTree* tree)  { fOutputTree = tree;  }
    void SetEvent(RecoEvent* event)  { fEvent      = event; }
    void SetStore(RecoStore* store)  { fStore      = store; }

protected:
    std::string fName;
    std::string fType;
    bool fEnabled;
    std::map<std::string, std::string> fParameters;

    // Pointers to I/O
    TFile*      fInputFile;   // Full input file (for geometry trees, etc.)
    TTree*      fInputTree;
    TTree*      fOutputTree;
    RecoEvent*  fEvent;
    RecoStore*  fStore;

    void Print(const std::string& message) const;
};

#endif
