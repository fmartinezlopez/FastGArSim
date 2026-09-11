//
// MacroParser.cc - Implementation of macro parser
//

#include "MacroParser.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

MacroParser::MacroParser()
{
}

MacroParser::~MacroParser()
{
}

bool MacroParser::ParseFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open macro file: " << filename << std::endl;
        return false;
    }

    std::cout << "\n Parsing macro file: " << filename << std::endl;

    fPath = filename;
    fText.clear();

    std::string line;
    int lineNumber = 0;
    while (std::getline(file, line)) {
        lineNumber++;
        // Kept verbatim, so that the output file can say exactly how it was
        // configured rather than only which file was pointed at
        fText += line;
        fText += "\n";
        ProcessLine(line);
    }

    file.close();

    std::cout << " Found " << fModuleConfigs.size() << " module(s) to run" << std::endl;

    return true;
}

void MacroParser::ProcessLine(const std::string& line)
{
    // Trim whitespace
    std::string trimmed = Trim(line);

    // Skip empty lines and comments
    if (trimmed.empty() || trimmed[0] == '#') {
        return;
    }

    // Tokenize the line
    std::vector<std::string> tokens = Tokenize(trimmed);
    if (tokens.empty()) {
        return;
    }

    // Parse commands
    if (tokens[0] == "/reco/addModule") {
        ParseModuleCommand(tokens);
    }
    else if (tokens[0].find("/reco/") == 0) {
        ParseParameterCommand(tokens);
    }
}

void MacroParser::ParseModuleCommand(const std::vector<std::string>& tokens)
{
    // Format: /reco/addModule <ModuleName> <ModuleType>
    // Example: /reco/addModule trackReco TrackRecoModule

    if (tokens.size() < 3) {
        std::cerr << "Warning: Invalid /reco/addModule command" << std::endl;
        return;
    }

    ModuleConfig config;
    config.name = tokens[1];
    config.type = tokens[2];
    config.enabled = true;

    fModuleConfigs.push_back(config);
    fCurrentModule = config.name;

    std::cout << "   Adding module: " << config.name << " (type: " << config.type << ")" << std::endl;
}

void MacroParser::ParseParameterCommand(const std::vector<std::string>& tokens)
{
    // Format: /reco/<module>/<parameter> <value>
    // Example: /reco/trackReco/minHits 5

    if (tokens.size() < 2) {
        return;
    }

    std::string command = tokens[0];
    std::string value = tokens[1];

    // Split command by '/'
    std::vector<std::string> parts;
    std::istringstream iss(command);
    std::string part;
    while (std::getline(iss, part, '/')) {
        if (!part.empty()) {
            parts.push_back(part);
        }
    }

    // Should have at least: reco, module, parameter
    if (parts.size() < 3) {
        return;
    }

    std::string moduleName = parts[1];
    std::string paramName = parts[2];

    // Find the module configuration
    for (auto& config : fModuleConfigs) {
        if (config.name == moduleName) {
            config.parameters[paramName] = value;
            return;
        }
    }

    // If module not found, it might be a global parameter
    if (moduleName == "global") {
        fGlobalParameters[paramName] = value;
    }
}

std::vector<std::string> MacroParser::Tokenize(const std::string& line)
{
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;

    while (iss >> token) {
        tokens.push_back(token);
    }

    return tokens;
}

std::string MacroParser::Trim(const std::string& str)
{
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }

    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::string MacroParser::GetGlobalParameter(const std::string& key, const std::string& defaultValue) const
{
    auto it = fGlobalParameters.find(key);
    if (it != fGlobalParameters.end()) {
        return it->second;
    }
    return defaultValue;
}

bool MacroParser::GetGlobalParameterBool(const std::string& key, bool defaultValue) const
{
    const std::string value = GetGlobalParameter(key, "");
    if (value.empty()) return defaultValue;

    if (value == "true"  || value == "True"  || value == "TRUE"  || value == "1") return true;
    if (value == "false" || value == "False" || value == "FALSE" || value == "0") return false;

    std::cerr << "Warning: '" << value << "' is not a yes or no answer for /reco/global/"
              << key << "; using " << (defaultValue ? "true" : "false") << "." << std::endl;
    return defaultValue;
}
