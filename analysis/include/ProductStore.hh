 /***************************************************************************
 * ProductStore.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Typed access to reconstruction products whose names are not known until
 *   the file is opened.
 *
 *   The reconstruction is modular, so which products a file holds depends on
 *   which modules were run. An analysis says what it needs once, in
 *   BeginJob():
 *
 *       fClusters = Require<std::vector<digi::TPCCluster>>("TPCClusters");
 *       fWaveforms = Optional<std::vector<digi::TPCWaveform>>("TPCWaveforms");
 *
 *   and gets back a handle. A required product that is not in the file, or is
 *   there with the wrong type, ends the job before the first event with a
 *   message naming what the file does hold. An optional one that is missing
 *   gives a handle that is simply never valid, so the code that uses it is
 *   guarded by one `if` and nothing else changes.
 *
 *       if (fWaveforms) {
 *           for (const digi::TPCWaveform& w : *fWaveforms) { ... }
 *       }
 *
 *   Handles stay valid for the whole job: they point at a slot the store
 *   owns, which is refilled by every GetEntry(). A handle also goes invalid
 *   by itself if the chain moves on to a file that lacks its branch, so an
 *   analysis reading a mixed set of files sees "not there" rather than the
 *   previous file's contents.
 *
 ***************************************************************************/

#ifndef ProductStore_hh
#define ProductStore_hh

#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "Rtypes.h"
#include "TClass.h"
#include "TTree.h"

#include "ProductSchema.hh"

namespace ana {

class ProductStore;

namespace detail {

// The store owns one slot per product it reads. The slot address is what
// ROOT writes into, so it has to outlive every handle pointing at it.
struct Slot {
    std::string key;
    std::string type;
    Bool_t present = kFALSE;   // is the branch in the file now open?

    virtual ~Slot() = default;
};

template <class T>
struct TypedSlot : Slot {
    T* object = nullptr;
    ~TypedSlot() override { delete object; }
};

} // namespace detail

/* -------------------------------------------------------------------------- */
/*                                   Handle                                   */
/* -------------------------------------------------------------------------- */

template <class T>
class Handle {
public:
    Handle() = default;

    // Is the product in the file currently being read?
    Bool_t IsValid() const { return fSlot && fSlot->present && fSlot->object; }
    explicit operator bool() const { return IsValid(); }

    // nullptr rather than throwing, for code that would rather test
    const T* Get() const { return IsValid() ? fSlot->object : nullptr; }

    // Throws if the product is not there. Reaching for a product without
    // either requiring it or testing the handle is the mistake being caught.
    const T& operator*() const { return *Checked(); }
    const T* operator->() const { return Checked(); }

    std::string Key() const { return fSlot ? fSlot->key : std::string(); }

private:
    friend class ProductStore;
    explicit Handle(detail::TypedSlot<T>* slot) : fSlot(slot) {}

    const T* Checked() const
    {
        if (!IsValid()) {
            throw std::runtime_error(
                "Product '" + Key() + "' is not available in this file. Ask for "
                "it with Require() to have the job stop at start-up instead, or "
                "test the handle before using it.");
        }
        return fSlot->object;
    }

    detail::TypedSlot<T>* fSlot = nullptr;
};

/* -------------------------------------------------------------------------- */
/*                                    Store                                   */
/* -------------------------------------------------------------------------- */

class ProductStore {
public:
    ProductStore() = default;
    ~ProductStore() = default;

    ProductStore(const ProductStore&) = delete;
    ProductStore& operator=(const ProductStore&) = delete;

    // Attach to the tree the products are read from. `schema`, when given, is
    // only used to make the messages better.
    void Connect(TTree* tree, const fastgarsim::ProductSchema* schema = nullptr);

    // Tell the store the chain has moved on to another file, so that handles
    // for branches the new file lacks stop reporting themselves as valid
    void OnNewTree(TTree* current);

    Bool_t Connected() const { return fTree != nullptr; }

    // What this file holds
    Bool_t Has(const std::string& key) const;
    std::string TypeOf(const std::string& key) const;
    std::vector<std::string> Keys() const;
    void Print(std::ostream& out) const;

    // Ask for a product. `required` records a message in Errors() when it is
    // missing; the caller decides what to do about it.
    template <class T>
    Handle<T> Bind(const std::string& key, Bool_t required);

    const std::vector<std::string>& Errors() const { return fErrors; }
    void ClearErrors() { fErrors.clear(); }

private:
    // Message listing what is available, for a product that is not
    std::string Available() const;

    TTree* fTree = nullptr;      // the chain, which branch addresses are set on
    TTree* fCurrent = nullptr;   // the file's tree, which is asked what it has
    const fastgarsim::ProductSchema* fSchema = nullptr;

    std::vector<std::unique_ptr<detail::Slot>> fSlots;
    std::vector<std::string> fErrors;
};

/* -------------------------------------------------------------------------- */

template <class T>
Handle<T> ProductStore::Bind(const std::string& key, Bool_t required)
{
    // Asking twice hands back the same slot rather than reading the branch
    // into two places
    for (const auto& slot : fSlots) {
        if (slot->key != key) continue;

        auto* typed = dynamic_cast<detail::TypedSlot<T>*>(slot.get());
        if (typed) return Handle<T>(typed);

        fErrors.push_back("product '" + key + "' is already being read as "
                          + slot->type + ", so it cannot also be read as "
                          + std::string(typeid(T).name()));
        return Handle<T>();
    }

    if (!fTree) {
        if (required) {
            fErrors.push_back("product '" + key + "' was required, but these files "
                              "hold no reconstruction tree. Run the reconstruction "
                              "over them first.");
        }
        return Handle<T>();
    }

    TClass* cls = TClass::GetClass(typeid(T));

    if (!fTree->GetBranch(key.c_str())) {
        if (required) {
            fErrors.push_back("product '" + key + "' was required, but it is not in "
                              "these files." + Available());
        }
        return Handle<T>();
    }

    // A branch of the right name holding something else is a mistake worth
    // naming, rather than letting ROOT read one type into another
    const std::string actual =
        fastgarsim::ProductSchema::BranchType(fTree, key.c_str());
    if (cls && !actual.empty() && actual != cls->GetName()) {
        fErrors.push_back("product '" + key + "' is a " + actual
                          + " in these files, but was asked for as a "
                          + cls->GetName() + ".");
        return Handle<T>();
    }

    auto slot = std::unique_ptr<detail::TypedSlot<T>>(new detail::TypedSlot<T>());
    slot->key = key;
    slot->type = cls ? cls->GetName() : actual;
    slot->object = new T();
    slot->present = kTRUE;

    detail::TypedSlot<T>* raw = slot.get();
    fSlots.push_back(std::move(slot));

    fTree->SetBranchAddress(key.c_str(), &raw->object);

    return Handle<T>(raw);
}

} // namespace ana

#endif
