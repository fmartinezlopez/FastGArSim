//
// MacroParser.hh - Parser for reconstruction configuration macros
//

#ifndef MacroParser_h
#define MacroParser_h 1

#include <string>
#include <vector>
#include <map>

class RecoModule;

struct ModuleConfig {
    std::string name;
    std::string type;
    bool enabled;
    std::map<std::string, std::string> parameters;

    ModuleConfig() : enabled(true) {}
};

class MacroParser {
public:
    MacroParser();
    ~MacroParser();

    // Parse a macro file
    bool ParseFile(const std::string& filename);

    // Get parsed module configurations
    const std::vector<ModuleConfig>& GetModuleConfigs() const { return fModuleConfigs; }

    // Get global parameters
    std::string GetGlobalParameter(const std::string& key, const std::string& defaultValue = "") const;
    bool GetGlobalParameterBool(const std::string& key, bool defaultValue) const;

    // The macro that was parsed, kept so that the output file can record how
    // it was produced
    const std::string& GetPath() const { return fPath; }
    const std::string& GetText() const { return fText; }

private:
    void ProcessLine(const std::string& line);
    void ParseModuleCommand(const std::vector<std::string>& tokens);
    void ParseParameterCommand(const std::vector<std::string>& tokens);

    std::vector<std::string> Tokenize(const std::string& line);
    std::string Trim(const std::string& str);

    std::vector<ModuleConfig> fModuleConfigs;
    std::map<std::string, std::string> fGlobalParameters;
    std::string fCurrentModule;
    std::string fPath;
    std::string fText;
};

#endif
