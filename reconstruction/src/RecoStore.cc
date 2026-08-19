//
// RecoStore.cc
//

#include "RecoStore.hh"

std::vector<std::string> RecoStore::GetKeys() const
{
    std::vector<std::string> keys;
    keys.reserve(fData.size());
    for (const auto& kv : fData) {
        keys.push_back(kv.first);
    }
    return keys;
}
