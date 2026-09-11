//
// ProductSchema.cc - Implementation of the file self-description
//

#include "ProductSchema.hh"

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <set>

#include "TBranch.h"
#include "TBranchElement.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TTree.h"

namespace fastgarsim {

const char* const ProductSchema::kSchemaTree = "Schema";
const char* const ProductSchema::kProvenanceTree = "Provenance";

namespace {

// Branch of std::string, the one type both ROOT and uproot read without a
// dictionary of ours
void StringBranch(TTree* tree, const char* name, std::string* address)
{
    tree->Branch(name, address);
}

bool ReadStringBranch(TTree* tree, const char* name, std::string** address)
{
    if (!tree->GetBranch(name)) return false;
    tree->SetBranchAddress(name, address);
    return true;
}

} // anonymous namespace

/* -------------------------------------------------------------------------- */
/*                                  Building                                  */
/* -------------------------------------------------------------------------- */

bool ProductSchema::IsMetadataTree(const std::string& name)
{
    return name == kSchemaTree || name == kProvenanceTree;
}

void ProductSchema::AddProduct(const ProductInfo& product)
{
    for (ProductInfo& existing : fProducts) {
        if (existing.tree == product.tree && existing.branch == product.branch) {
            existing = product;
            return;
        }
    }
    fProducts.push_back(product);
}

void ProductSchema::AddJob(const JobInfo& job)
{
    fJobs.push_back(job);
}

void ProductSchema::Annotate(const ProductSchema& other)
{
    for (ProductInfo& product : fProducts) {
        const ProductInfo* source = other.Find(product.branch, product.tree);
        if (!source) continue;
        if (!source->producer.empty())     product.producer     = source->producer;
        if (!source->producerType.empty()) product.producerType = source->producerType;
        if (!source->parameters.empty())   product.parameters   = source->parameters;
    }
}

void ProductSchema::RemoveTree(const std::string& tree)
{
    fProducts.erase(std::remove_if(fProducts.begin(), fProducts.end(),
                                   [&tree](const ProductInfo& p) { return p.tree == tree; }),
                    fProducts.end());
}

/* -------------------------------------------------------------------------- */
/*                                   Reading                                  */
/* -------------------------------------------------------------------------- */

const ProductInfo* ProductSchema::Find(const std::string& branch,
                                       const std::string& tree) const
{
    for (const ProductInfo& product : fProducts) {
        if (product.branch != branch) continue;
        if (!tree.empty() && product.tree != tree) continue;
        return &product;
    }
    return nullptr;
}

std::vector<ProductInfo> ProductSchema::InTree(const std::string& tree) const
{
    std::vector<ProductInfo> selected;
    for (const ProductInfo& product : fProducts) {
        if (product.tree == tree) selected.push_back(product);
    }
    return selected;
}

std::vector<std::string> ProductSchema::TreeNames() const
{
    std::vector<std::string> names;
    for (const ProductInfo& product : fProducts) {
        if (std::find(names.begin(), names.end(), product.tree) == names.end()) {
            names.push_back(product.tree);
        }
    }
    return names;
}

std::string ProductSchema::Fingerprint() const
{
    // Sorted, so that the order the products were declared in does not matter
    std::set<std::string> entries;
    for (const ProductInfo& product : fProducts) {
        entries.insert(product.tree + "/" + product.branch + ":" + product.type);
    }

    // FNV-1a, 64 bit
    unsigned long long hash = 1469598103934665603ULL;
    for (const std::string& entry : entries) {
        for (const char c : entry) {
            hash ^= static_cast<unsigned char>(c);
            hash *= 1099511628211ULL;
        }
        hash ^= '\n';
        hash *= 1099511628211ULL;
    }

    char text[17];
    std::snprintf(text, sizeof(text), "%016llx", hash);
    return std::string(text);
}

void ProductSchema::Print() const { Print(std::cout); }

void ProductSchema::Print(std::ostream& out) const
{
    if (fProducts.empty()) {
        out << "No products described." << std::endl;
    } else {
        // Column widths from the content, so the table stays readable
        size_t branchWidth = 6;
        size_t typeWidth = 4;
        for (const ProductInfo& product : fProducts) {
            branchWidth = std::max(branchWidth, product.branch.size());
            typeWidth   = std::max(typeWidth,   product.type.size());
        }

        for (const std::string& tree : TreeNames()) {
            out << "\nTree '" << tree << "'\n";
            for (const ProductInfo& product : InTree(tree)) {
                out << "    " << std::left << std::setw(static_cast<int>(branchWidth))
                    << product.branch << "  "
                    << std::setw(static_cast<int>(typeWidth)) << product.type;
                if (!product.producer.empty()) {
                    out << "  <- " << product.producer;
                    if (!product.producerType.empty()) {
                        out << " (" << product.producerType << ")";
                    }
                }
                out << "\n";
            }
        }
        out << std::right << "\nFingerprint: " << Fingerprint() << std::endl;
    }

    for (const JobInfo& job : fJobs) {
        out << "\n" << job.stage << " pass, " << job.timestamp << "\n";
        if (!job.input.empty())   out << "    input:   " << job.input << "\n";
        if (!job.macro.empty())   out << "    macro:   " << job.macro << "\n";
        if (!job.modules.empty()) out << "    modules: " << job.modules << "\n";
    }
    out << std::flush;
}

/* -------------------------------------------------------------------------- */
/*                                     I/O                                    */
/* -------------------------------------------------------------------------- */

bool ProductSchema::Write(TDirectory* directory) const
{
    if (!directory) return false;

    TDirectory* previous = gDirectory;
    directory->cd();

    {
        TTree schema(kSchemaTree, "FastGArSim product schema");

        std::string tree, branch, type, producer, producerType, parameters;
        StringBranch(&schema, "tree",         &tree);
        StringBranch(&schema, "branch",       &branch);
        StringBranch(&schema, "type",         &type);
        StringBranch(&schema, "producer",     &producer);
        StringBranch(&schema, "producerType", &producerType);
        StringBranch(&schema, "parameters",   &parameters);

        for (const ProductInfo& product : fProducts) {
            tree         = product.tree;
            branch       = product.branch;
            type         = product.type;
            producer     = product.producer;
            producerType = product.producerType;
            parameters   = product.parameters;
            schema.Fill();
        }

        // kOverwrite replaces an earlier cycle rather than adding to it, so a
        // file that has been through reconstruction twice has one schema
        schema.Write(nullptr, TObject::kOverwrite);
    }

    {
        TTree provenance(kProvenanceTree, "FastGArSim processing history");

        std::string stage, timestamp, input, macro, macroText, modules;
        StringBranch(&provenance, "stage",     &stage);
        StringBranch(&provenance, "timestamp", &timestamp);
        StringBranch(&provenance, "input",     &input);
        StringBranch(&provenance, "macro",     &macro);
        StringBranch(&provenance, "macroText", &macroText);
        StringBranch(&provenance, "modules",   &modules);

        for (const JobInfo& job : fJobs) {
            stage     = job.stage;
            timestamp = job.timestamp;
            input     = job.input;
            macro     = job.macro;
            macroText = job.macroText;
            modules   = job.modules;
            provenance.Fill();
        }

        provenance.Write(nullptr, TObject::kOverwrite);
    }

    if (previous) previous->cd();
    return true;
}

ProductSchema ProductSchema::Read(TDirectory* directory)
{
    ProductSchema schema;
    if (!directory) return schema;

    if (TTree* tree = dynamic_cast<TTree*>(directory->Get(kSchemaTree))) {
        std::string treeName, branch, type, producer, producerType, parameters;
        std::string* pTreeName     = &treeName;
        std::string* pBranch       = &branch;
        std::string* pType         = &type;
        std::string* pProducer     = &producer;
        std::string* pProducerType = &producerType;
        std::string* pParameters   = &parameters;

        const bool ok = ReadStringBranch(tree, "tree",   &pTreeName)
                     && ReadStringBranch(tree, "branch", &pBranch)
                     && ReadStringBranch(tree, "type",   &pType);
        ReadStringBranch(tree, "producer",     &pProducer);
        ReadStringBranch(tree, "producerType", &pProducerType);
        ReadStringBranch(tree, "parameters",   &pParameters);

        if (ok) {
            for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
                tree->GetEntry(i);
                ProductInfo product;
                product.tree         = *pTreeName;
                product.branch       = *pBranch;
                product.type         = *pType;
                product.producer     = *pProducer;
                product.producerType = *pProducerType;
                product.parameters   = *pParameters;
                schema.AddProduct(product);
            }
        }
        tree->ResetBranchAddresses();
    }

    if (TTree* tree = dynamic_cast<TTree*>(directory->Get(kProvenanceTree))) {
        std::string stage, timestamp, input, macro, macroText, modules;
        std::string* pStage     = &stage;
        std::string* pTimestamp = &timestamp;
        std::string* pInput     = &input;
        std::string* pMacro     = &macro;
        std::string* pMacroText = &macroText;
        std::string* pModules   = &modules;

        ReadStringBranch(tree, "stage",     &pStage);
        ReadStringBranch(tree, "timestamp", &pTimestamp);
        ReadStringBranch(tree, "input",     &pInput);
        ReadStringBranch(tree, "macro",     &pMacro);
        ReadStringBranch(tree, "macroText", &pMacroText);
        ReadStringBranch(tree, "modules",   &pModules);

        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            JobInfo job;
            job.stage     = *pStage;
            job.timestamp = *pTimestamp;
            job.input     = *pInput;
            job.macro     = *pMacro;
            job.macroText = *pMacroText;
            job.modules   = *pModules;
            schema.AddJob(job);
        }
        tree->ResetBranchAddresses();
    }

    return schema;
}

ProductSchema ProductSchema::Describe(TDirectory* directory)
{
    ProductSchema schema;
    if (!directory) return schema;

    TIter next(directory->GetListOfKeys());
    std::set<std::string> seen;

    while (TKey* key = dynamic_cast<TKey*>(next())) {
        // Only the highest cycle of each name, and only trees
        const std::string name = key->GetName();
        if (!seen.insert(name).second) continue;
        if (IsMetadataTree(name)) continue;

        const TClass* cls = TClass::GetClass(key->GetClassName());
        if (!cls || !cls->InheritsFrom(TTree::Class())) continue;

        TTree* tree = dynamic_cast<TTree*>(directory->Get(name.c_str()));
        if (!tree) continue;

        TObjArray* branches = tree->GetListOfBranches();
        for (int i = 0; i < branches->GetEntriesFast(); ++i) {
            TBranch* branch = dynamic_cast<TBranch*>(branches->At(i));
            if (!branch) continue;

            ProductInfo product;
            product.tree = name;
            product.branch = branch->GetName();
            product.type = BranchType(tree, branch->GetName());
            schema.AddProduct(product);
        }
    }

    return schema;
}

ProductSchema ProductSchema::Load(TDirectory* directory)
{
    ProductSchema schema = Describe(directory);
    const ProductSchema recorded = Read(directory);
    schema.Annotate(recorded);
    for (const JobInfo& job : recorded.Jobs()) schema.AddJob(job);
    return schema;
}

/* -------------------------------------------------------------------------- */
/*                                  Utilities                                 */
/* -------------------------------------------------------------------------- */

std::string ProductSchema::BranchType(TTree* tree, const char* branchName)
{
    if (!tree) return "";

    TBranch* branch = tree->GetBranch(branchName);
    if (!branch) return "";

    if (TBranchElement* element = dynamic_cast<TBranchElement*>(branch)) {
        const char* className = element->GetClassName();
        if (className && *className) return className;
    }

    // A branch of a fundamental type: take the leaf's type
    if (TLeaf* leaf = branch->GetLeaf(branchName)) {
        return leaf->GetTypeName();
    }
    if (branch->GetListOfLeaves()->GetEntriesFast() == 1) {
        TLeaf* first = dynamic_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
        if (first) return first->GetTypeName();
    }

    return "";
}

std::string ProductSchema::Now()
{
    const std::time_t now = std::time(nullptr);
    std::tm parts{};
#ifdef _WIN32
    localtime_s(&parts, &now);
#else
    localtime_r(&now, &parts);
#endif
    char text[32];
    std::strftime(text, sizeof(text), "%Y-%m-%d %H:%M:%S", &parts);
    return std::string(text);
}

} // namespace fastgarsim
