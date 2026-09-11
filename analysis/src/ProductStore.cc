 /***************************************************************************
 * ProductStore.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   The parts of the product store that do not depend on the product type.
 *
 ***************************************************************************/

#include "ProductStore.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>

#include "TBranch.h"
#include "TDirectory.h"
#include "TObjArray.h"

namespace ana {

void ProductStore::Connect(TTree* tree, const fastgarsim::ProductSchema* schema)
{
    fTree = tree;
    fCurrent = tree ? tree->GetTree() : nullptr;   // the chain's current file
    fSchema = schema;
}

void ProductStore::OnNewTree(TTree* current)
{
    fCurrent = current;
    if (!current) return;

    for (auto& slot : fSlots) {
        const Bool_t present = current->GetBranch(slot->key.c_str()) != nullptr;
        if (present == slot->present) continue;

        slot->present = present;
        if (!present) {
            const TDirectory* file = current->GetDirectory();
            std::cout << "Warning: '" << (file ? file->GetName() : "the next file")
                      << "' has no '" << slot->key << "' branch, so that product "
                      << "is not available for its events." << std::endl;
        }
    }
}

Bool_t ProductStore::Has(const std::string& key) const
{
    const TTree* tree = fCurrent ? fCurrent : fTree;
    return tree && const_cast<TTree*>(tree)->GetBranch(key.c_str()) != nullptr;
}

std::string ProductStore::TypeOf(const std::string& key) const
{
    TTree* tree = fCurrent ? fCurrent : fTree;
    if (!tree) return std::string();
    return fastgarsim::ProductSchema::BranchType(tree, key.c_str());
}

std::vector<std::string> ProductStore::Keys() const
{
    std::vector<std::string> keys;

    TTree* tree = fCurrent ? fCurrent : fTree;
    if (!tree) return keys;

    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntriesFast(); ++i) {
        keys.push_back(branches->At(i)->GetName());
    }
    return keys;
}

std::string ProductStore::Available() const
{
    const std::vector<std::string> keys = Keys();
    if (keys.empty()) return " They hold no reconstruction products at all.";

    std::string message = " Available:";
    for (const std::string& key : keys) message += " " + key;
    return message;
}

void ProductStore::Print(std::ostream& out) const
{
    const std::vector<std::string> keys = Keys();
    if (keys.empty()) {
        out << "No reconstruction products in these files." << std::endl;
        return;
    }

    size_t keyWidth = 0;
    size_t typeWidth = 0;
    for (const std::string& key : keys) {
        keyWidth = std::max(keyWidth, key.size());
        typeWidth = std::max(typeWidth, TypeOf(key).size());
    }

    out << "Reconstruction products:" << std::endl;
    for (const std::string& key : keys) {
        out << "    " << std::left << std::setw(static_cast<int>(keyWidth)) << key
            << "  " << std::setw(static_cast<int>(typeWidth)) << TypeOf(key);

        // The producing module, when the file recorded it
        if (fSchema) {
            if (const fastgarsim::ProductInfo* product = fSchema->Find(key)) {
                if (!product->producer.empty()) {
                    out << "  <- " << product->producer;
                }
            }
        }
        out << "\n";
    }
    out << std::right << std::flush;
}

} // namespace ana
