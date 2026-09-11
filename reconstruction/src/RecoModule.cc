//
// RecoModule.cc - Implementation of base reconstruction module
//

#include "RecoModule.hh"
#include <iostream>
#include <sstream>

RecoModule::RecoModule(const std::string& name)
    : fName(name), fType(name), fEnabled(true),
      fInputFile(nullptr), fInputTree(nullptr),
      fOutputTree(nullptr), fEvent(nullptr), fStore(nullptr)
{
}

RecoModule::~RecoModule()
{
}

void RecoModule::SetParameter(const std::string& key, const std::string& value)
{
    fParameters[key] = value;
}

std::string RecoModule::GetParameter(const std::string& key, const std::string& defaultValue) const
{
    auto it = fParameters.find(key);
    if (it != fParameters.end()) {
        return it->second;
    }
    return defaultValue;
}

int RecoModule::GetParameterInt(const std::string& key, int defaultValue) const
{
    std::string value = GetParameter(key, "");
    if (value.empty()) return defaultValue;

    std::istringstream iss(value);
    int result;
    iss >> result;
    return result;
}

double RecoModule::GetParameterDouble(const std::string& key, double defaultValue) const
{
    std::string value = GetParameter(key, "");
    if (value.empty()) return defaultValue;

    std::istringstream iss(value);
    double result;
    iss >> result;
    return result;
}

bool RecoModule::GetParameterBool(const std::string& key, bool defaultValue) const
{
    std::string value = GetParameter(key, "");
    if (value.empty()) return defaultValue;

    // Handle various boolean representations
    if (value == "true" || value == "True" || value == "TRUE" || value == "1") {
        return true;
    }
    if (value == "false" || value == "False" || value == "FALSE" || value == "0") {
        return false;
    }

    return defaultValue;
}

void RecoModule::Print(const std::string& message) const
{
    std::cout << " [" << fName << "] " << message << std::endl;
}
