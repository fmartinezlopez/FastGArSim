//
// ProductSchema.hh - Self-description of what a FastGArSim ROOT file contains
//
// The reconstruction is modular: which products a file holds depends on which
// modules were configured, so the set of branches is not known in advance. A
// file therefore carries its own description.
//
// Two trees hold it:
//
//   Schema      one entry per product: the tree and branch it lives in, its
//               C++ type, and which module produced it with what parameters
//   Provenance  one entry per processing pass: when it ran, on what input,
//               with which macro (the macro text is stored verbatim)
//
// Both are written by the reconstruction and read by the analysis, so an
// analysis can be told at start-up that a product it needs is missing rather
// than discovering it mid-loop. They are plain trees of std::string branches,
// so uproot and `root -l` read them without any dictionary.
//
// Files written before this existed, and files from the simulation, have no
// Schema tree. Describe() reconstructs one by walking the trees actually
// present, so every reader can assume a schema exists.
//

#ifndef ProductSchema_hh
#define ProductSchema_hh

#include <iosfwd>
#include <string>
#include <vector>

class TDirectory;
class TTree;

namespace fastgarsim {

// One branch of one tree
struct ProductInfo {
    std::string tree;          // TTree holding it, e.g. "Reco"
    std::string branch;        // Branch name, e.g. "TPCClusters"
    std::string type;          // C++ type, e.g. "vector<digi::TPCCluster>"
    std::string producer;      // Module instance name, empty if not known
    std::string producerType;  // Module class name, empty if not known
    std::string parameters;    // "key=value; key=value", empty if not known
};

// One pass of a processing stage over the data
struct JobInfo {
    std::string stage;      // "reconstruction", "ntuple", ...
    std::string timestamp;  // ISO 8601, local time
    std::string input;      // Input file the pass read
    std::string macro;      // Configuration macro path, if there was one
    std::string macroText;  // Its contents, so the pass can be reproduced
    std::string modules;    // "name:Type, name:Type, ..."
};

class ProductSchema {
public:
    static const char* const kSchemaTree;
    static const char* const kProvenanceTree;

    // True for the trees this class writes, which readers should skip when
    // they enumerate the data trees of a file
    static bool IsMetadataTree(const std::string& name);

    /* ----------------------------- Building ---------------------------- */

    // Adding a product that is already described replaces the description
    void AddProduct(const ProductInfo& product);
    void AddJob(const JobInfo& job);

    // Copy producer, producerType and parameters from `other` onto matching
    // products of this schema, leaving everything else alone. Used to label
    // the branches discovered in a file with the module that wrote them.
    void Annotate(const ProductSchema& other);

    // Drop every product of one tree, e.g. before re-writing that tree
    void RemoveTree(const std::string& tree);

    /* ------------------------------ Reading ---------------------------- */

    const std::vector<ProductInfo>& Products() const { return fProducts; }
    const std::vector<JobInfo>& Jobs() const { return fJobs; }

    bool Empty() const { return fProducts.empty() && fJobs.empty(); }

    // Look a product up by branch name. With `tree` empty the first product
    // of that name in any tree is returned; nullptr when there is none.
    const ProductInfo* Find(const std::string& branch,
                            const std::string& tree = "") const;

    bool Has(const std::string& branch, const std::string& tree = "") const
    { return Find(branch, tree) != nullptr; }

    std::vector<ProductInfo> InTree(const std::string& tree) const;
    std::vector<std::string> TreeNames() const;

    // Stable hash over (tree, branch, type) of every product. Two files with
    // the same fingerprint hold the same products and can be chained.
    std::string Fingerprint() const;

    void Print(std::ostream& out) const;
    void Print() const;

    /* ------------------------------- I/O ------------------------------- */

    // Write the Schema and Provenance trees into `directory`, replacing any
    // that are already there. Returns false if it could not be written.
    bool Write(TDirectory* directory) const;

    // Read the Schema and Provenance trees. The result is empty when the file
    // has none, which is not an error.
    static ProductSchema Read(TDirectory* directory);

    // Build a schema by walking the trees in the file. Producer fields are
    // left empty: the file itself does not say who wrote what.
    static ProductSchema Describe(TDirectory* directory);

    // Describe() the file and then Annotate() it with whatever Read() finds,
    // so the result covers every branch really present while keeping the
    // recorded provenance. This is what readers should use.
    static ProductSchema Load(TDirectory* directory);

    /* ---------------------------- Utilities ---------------------------- */

    // C++ type name of one branch: the class name for an object branch, the
    // leaf type for a branch of a fundamental type.
    static std::string BranchType(TTree* tree, const char* branchName);

    // Local time, formatted as YYYY-MM-DD HH:MM:SS
    static std::string Now();

private:
    std::vector<ProductInfo> fProducts;
    std::vector<JobInfo> fJobs;
};

} // namespace fastgarsim

#endif
