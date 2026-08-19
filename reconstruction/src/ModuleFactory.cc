//
// ModuleFactory.cc
//

#include "ModuleFactory.hh"
#include "RecoModule.hh"

ModuleFactory& ModuleFactory::Instance()
{
    static ModuleFactory instance;
    return instance;
}

bool ModuleFactory::Register(const std::string& typeName, Creator creator)
{
    fRegistry[typeName] = std::move(creator);
    return true;
}

RecoModule* ModuleFactory::Create(const std::string& typeName) const
{
    auto it = fRegistry.find(typeName);
    if (it == fRegistry.end()) {
        return nullptr;
    }
    return it->second();
}

bool ModuleFactory::IsRegistered(const std::string& typeName) const
{
    return fRegistry.count(typeName) > 0;
}

std::vector<std::string> ModuleFactory::GetRegisteredTypes() const
{
    std::vector<std::string> types;
    types.reserve(fRegistry.size());
    for (const auto& kv : fRegistry) {
        types.push_back(kv.first);
    }
    return types;
}
