//
// TreeFlattener.hh - Turn a tree of objects into a tree of flat vectors
//
// Nothing here knows about any FastGArSim class. The layout of the output is
// derived from the ROOT dictionaries of whatever the input tree happens to
// contain, so the same code flattens the simulation's root::Event, the
// reconstruction's digi::TPCCluster collections, and any product a module
// added afterwards.
//
// THE RULE
//
// Every branch is walked member by member. A member of a fundamental type
// becomes one output column. A member that is itself a class is walked
// recursively, its name prefixed onto its members'. A member that is a
// collection starts a new *group*: its own set of columns, with one row per
// element, plus an index column saying which row of the parent group each
// element came from. That one rule covers arbitrary nesting -- particles
// inside an event, hits inside a particle, track IDs inside a hit.
//
// NAMING
//
// Columns are named <group>_<member>, the group being the branch name for a
// collection branch and the chain of member names below it:
//
//   Reco tree,   branch TPCClusters   ->  TPCClusters_energy, TPCClusters_x, ...
//                                         TPCClusters_trackIDs, and its
//                                         TPCClusters_trackIDs_parent index
//   Events tree, branch Event         ->  eventID
//                                         particles_pdgCode, ...
//                                         particles_tpcHits_energyDeposit, and
//                                         its particles_tpcHits_parent index
//
// The name of a branch holding a single object is dropped, since it is the
// per-event wrapper and repeating it in every column reads badly. Pass
// keepBranchPrefix to keep it.
//
// TYPES
//
// Integral members become Int_t or Long64_t columns, floating point ones
// Float_t or Double_t, and TString and std::string members std::string. The
// column is a std::vector of that, or the bare scalar for the members of a
// single-object branch, which occur once per entry.
//

#ifndef TreeFlattener_hh
#define TreeFlattener_hh

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Rtypes.h"
#include "TDataType.h"

#include "ProductSchema.hh"

class TClass;
class TTree;
class TVirtualCollectionProxy;

namespace fastgarsim {

// One output column: a std::vector of some fundamental type, or the bare
// scalar for a per-event value. Defined in the implementation, since nothing
// outside it needs to know how a column stores or converts its values.
class FlatColumn;

struct FlattenOptions {
    // How deep to follow nested collections. Depth 1 is the members of the
    // branch itself, 2 its collections' members, and so on. The default
    // reaches the hits of a particle and the track IDs of a hit.
    int maxDepth = 4;

    // Keep the branch name as a prefix even for single-object branches
    bool keepBranchPrefix = false;

    // Branches to flatten. Empty means every branch; exclude wins over
    // include. Names are matched exactly.
    std::vector<std::string> include;
    std::vector<std::string> exclude;

    // Group prefix rewrites, matched on the whole prefix:
    // {"particles_tpcHits", "tpcHit"} turns particles_tpcHits_energyDeposit
    // into tpcHit_energyDeposit. Children inherit the new prefix.
    std::map<std::string, std::string> rename;

    bool verbose = true;
};

class TreeFlattener {
public:
    explicit TreeFlattener(const FlattenOptions& options = FlattenOptions());
    ~TreeFlattener();

    TreeFlattener(const TreeFlattener&) = delete;
    TreeFlattener& operator=(const TreeFlattener&) = delete;

    // Work out the output layout from `input` and attach to its branches.
    // Returns false if nothing could be flattened.
    bool Connect(TTree* input);

    // Book the columns on `output`, which must live in the target file
    void Book(TTree* output);

    // Read one entry of the input and append it to the columns. Returns the
    // number of bytes read, or -1 on failure.
    Long64_t FillEntry(Long64_t entry);

    // Description of the columns booked by Book(), for the file's schema
    std::vector<ProductInfo> Schema(const std::string& outputTreeName) const;

    void PrintLayout(std::ostream& out) const;

    size_t NColumns() const;
    size_t NGroups() const;

private:
    struct Field;
    struct Group;

    // Layout construction
    void AddBranch(TTree* input, const char* name);
    void AddObjectBranch(TTree* input, const char* name, TClass* cls);
    void AddPlainBranch(TTree* input, const char* name, EDataType type);
    void AddMembers(int group, TClass* cls, Long_t offset,
                    const std::string& prefix, int depth);
    int  AddCollection(int parent, const std::string& prefix, Long_t offset,
                       TClass* collectionClass, int depth);
    void AddField(int group, const std::string& name, Long_t offset,
                  EDataType type, int kind);
    int  NewGroup(int parent, const std::string& prefix, bool scalar);

    // Per-entry filling
    void FillCollection(int group, const char* collection, Int_t parentRow);
    void FillObject(int group, const char* object, Int_t parentRow);

    std::string Rename(const std::string& name) const;
    bool Selected(const char* branchName) const;
    // Records a column name, or renames it and says so if it is taken
    std::string Claim(std::string name);

    FlattenOptions fOptions;
    TTree* fInput = nullptr;

    std::vector<std::unique_ptr<Group>> fGroups;

    // Objects the object branches are read into, owned here
    struct BranchObject {
        std::string name;
        TClass* cls = nullptr;
        void* object = nullptr;
        int group = -1;
    };
    // Held by pointer: the addresses are handed to TTree::SetBranchAddress
    // and have to stay put as more branches are added
    std::vector<std::unique_ptr<BranchObject>> fBranchObjects;

    // Branches of a fundamental type, copied straight into a scalar column
    struct PlainBranch {
        std::string name;
        EDataType type = kNoType_t;
        Long64_t buffer = 0;             // wide enough for any of them
        std::unique_ptr<FlatColumn> column;
    };
    std::vector<std::unique_ptr<PlainBranch>> fPlainBranches;

    std::vector<std::string> fColumnNames;  // every name claimed so far
};

} // namespace fastgarsim

#endif
