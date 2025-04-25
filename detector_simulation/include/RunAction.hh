#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"
#include <map>

class G4Run;
class AnalysisManager;

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();
    
    // Methods from base class
    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);
    
    // Method to set output file name
    void SetSaveOutput(G4bool save);
    void SetOutputFileName(G4String& name);

    // Getter methods
    G4bool GetSaveOutput() const { return fSaveOutput; }
    
private:

    // Helper methods
    void DefineCommands();
    
    // Output options
    G4bool fSaveOutput;
    G4String fOutputFileName;
    
    // Analysis manager
    AnalysisManager* fAnalysisManager;

    // Messenger for macro commands
    G4GenericMessenger* fMessenger;
};

#endif