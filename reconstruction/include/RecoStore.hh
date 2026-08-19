//
// RecoStore.hh - Central data store for inter-module object sharing
//
// Modules register their output objects during Initialize() and downstream
// modules retrieve them. Pointers are stable across events; each module's
// Execute() clears and refills its own containers in place.
//

#ifndef RecoStore_h
#define RecoStore_h 1

#include <any>
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <typeindex>

class RecoStore {
public:
    RecoStore() = default;
    ~RecoStore() = default;

    // Register a stable pointer to a module-owned object (call in Initialize())
    template<typename T>
    void Register(const std::string& key, T* ptr)
    {
        fData[key] = ptr;
        fTypeNames[key] = typeid(T).name();
    }

    // Retrieve a registered pointer (call in Initialize() to cache for Execute())
    template<typename T>
    T* Get(const std::string& key) const
    {
        auto it = fData.find(key);
        if (it == fData.end()) {
            throw std::runtime_error(
                "RecoStore: object '" + key + "' not found. "
                "Check that the producing module runs before this one.");
        }
        try {
            return std::any_cast<T*>(it->second);
        } catch (const std::bad_any_cast&) {
            throw std::runtime_error(
                "RecoStore: object '" + key + "' exists but has wrong type. "
                "Expected: " + typeid(T).name() +
                ", stored: " + fTypeNames.at(key));
        }
    }

    bool Has(const std::string& key) const { return fData.count(key) > 0; }

    std::vector<std::string> GetKeys() const;

private:
    std::unordered_map<std::string, std::any> fData;
    std::unordered_map<std::string, std::string> fTypeNames;
};

#endif
