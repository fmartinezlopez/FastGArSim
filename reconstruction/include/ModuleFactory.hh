//
// ModuleFactory.hh - Self-registering factory for reconstruction modules
//
// Each module registers itself via a static initializer in its .cc file:
//
//   namespace {
//     bool kRegistered = ModuleFactory::Instance().Register(
//       "MyModule", []() -> RecoModule* { return new MyModule(); });
//   }
//
// RecoManager then creates modules by type name without knowing about them:
//
//   RecoModule* m = ModuleFactory::Instance().Create("MyModule");
//

#ifndef ModuleFactory_h
#define ModuleFactory_h 1

#include <functional>
#include <map>
#include <string>
#include <vector>

class RecoModule;

class ModuleFactory {
public:
    using Creator = std::function<RecoModule*()>;

    static ModuleFactory& Instance();

    // Register a module type. Returns true so it can drive a static bool init.
    bool Register(const std::string& typeName, Creator creator);

    // Create a module by type name. Returns nullptr if type is unknown.
    RecoModule* Create(const std::string& typeName) const;

    bool IsRegistered(const std::string& typeName) const;

    std::vector<std::string> GetRegisteredTypes() const;

private:
    ModuleFactory() = default;
    std::map<std::string, Creator> fRegistry;
};

#endif
