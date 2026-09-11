//
// TreeFlattener.cc - Implementation of the reflection-driven flattener
//
// The layout is worked out once, in Connect(), by walking the ROOT dictionary
// of every branch. After that FillEntry() only follows the offsets it has
// already resolved, so the per-event cost is a member read per column and no
// dictionary lookups.
//

#include "TreeFlattener.hh"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <type_traits>

#include "TBaseClass.h"
#include "TBranch.h"
#include "TClass.h"
#include "TDataMember.h"
#include "TList.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TVirtualCollectionProxy.h"

namespace fastgarsim {

namespace {

// How a member is read out of an object
enum EMemberKind { kNumericMember, kTStringMember, kStdStringMember };

/* -------------------------------------------------------------------------- */
/*                          Reading a member generically                      */
/* -------------------------------------------------------------------------- */

Long64_t ReadIntegral(const void* address, EDataType type)
{
    switch (type) {
        case kBool_t:     return *static_cast<const Bool_t*>(address) ? 1 : 0;
        case kChar_t:     return *static_cast<const Char_t*>(address);
        case kUChar_t:    return *static_cast<const UChar_t*>(address);
        case kShort_t:    return *static_cast<const Short_t*>(address);
        case kUShort_t:   return *static_cast<const UShort_t*>(address);
        case kInt_t:      return *static_cast<const Int_t*>(address);
        case kUInt_t:     return *static_cast<const UInt_t*>(address);
        case kLong_t:     return *static_cast<const Long_t*>(address);
        case kULong_t:    return static_cast<Long64_t>(*static_cast<const ULong_t*>(address));
        case kLong64_t:   return *static_cast<const Long64_t*>(address);
        case kULong64_t:  return static_cast<Long64_t>(*static_cast<const ULong64_t*>(address));
        case kFloat_t:
        case kFloat16_t:  return static_cast<Long64_t>(*static_cast<const Float_t*>(address));
        case kDouble_t:
        case kDouble32_t: return static_cast<Long64_t>(*static_cast<const Double_t*>(address));
        default:          return 0;
    }
}

Double_t ReadReal(const void* address, EDataType type)
{
    switch (type) {
        case kFloat_t:
        case kFloat16_t:  return *static_cast<const Float_t*>(address);
        case kDouble_t:
        case kDouble32_t: return *static_cast<const Double_t*>(address);
        default:          return static_cast<Double_t>(ReadIntegral(address, type));
    }
}

std::string ReadString(const void* address, int kind)
{
    if (kind == kTStringMember) {
        return std::string(static_cast<const TString*>(address)->Data());
    }
    return *static_cast<const std::string*>(address);
}

// Column type a member of the given type is stored in. Keeping the number of
// column types small keeps the number of branch types in the ntuple small,
// without ever narrowing a value.
enum EColumnType { kBoolColumn, kIntColumn, kLongColumn, kFloatColumn,
                   kDoubleColumn, kStringColumn };

EColumnType ColumnTypeFor(EDataType type, int kind)
{
    if (kind != kNumericMember) return kStringColumn;

    switch (type) {
        case kBool_t:
            return kBoolColumn;
        case kChar_t: case kUChar_t: case kShort_t: case kUShort_t:
        case kInt_t:  case kUInt_t:
            return kIntColumn;
        case kLong_t: case kULong_t: case kLong64_t: case kULong64_t:
            return kLongColumn;
        case kFloat_t: case kFloat16_t:
            return kFloatColumn;
        default:
            return kDoubleColumn;
    }
}

bool IsStdString(const TClass* cls)
{
    if (!cls) return false;
    const std::string name = cls->GetName();
    return name == "string" || name == "std::string";
}

bool IsTString(const TClass* cls)
{
    return cls && std::strcmp(cls->GetName(), "TString") == 0;
}

} // anonymous namespace

/* -------------------------------------------------------------------------- */
/*                                   Columns                                  */
/* -------------------------------------------------------------------------- */

class FlatColumn {
public:
    virtual ~FlatColumn() = default;
    virtual void Book(TTree* tree, const std::string& name, bool scalar) = 0;
    virtual void Clear() = 0;
    virtual void Append(const void* address, EDataType type, int kind) = 0;
    virtual void Set(const void* address, EDataType type, int kind) = 0;
    virtual std::string TypeName(bool scalar) const = 0;
};

namespace {

template <class T>
struct NumericColumn : public FlatColumn {
    std::vector<T> values;
    T scalar{};

    static T Convert(const void* address, EDataType type)
    {
        if constexpr (std::is_floating_point<T>::value) {
            return static_cast<T>(ReadReal(address, type));
        } else {
            return static_cast<T>(ReadIntegral(address, type));
        }
    }

    void Book(TTree* tree, const std::string& name, bool asScalar) override
    {
        if (asScalar) tree->Branch(name.c_str(), &scalar);
        else          tree->Branch(name.c_str(), &values);
    }

    void Clear() override { values.clear(); }

    void Append(const void* address, EDataType type, int) override
    { values.push_back(Convert(address, type)); }

    void Set(const void* address, EDataType type, int) override
    { scalar = Convert(address, type); }

    std::string TypeName(bool asScalar) const override
    {
        const char* name = "Double_t";
        if (std::is_same<T, Bool_t>::value)        name = "Bool_t";
        else if (std::is_same<T, Int_t>::value)    name = "Int_t";
        else if (std::is_same<T, Long64_t>::value) name = "Long64_t";
        else if (std::is_same<T, Float_t>::value)  name = "Float_t";
        return asScalar ? std::string(name) : "vector<" + std::string(name) + ">";
    }
};

struct StringColumn : public FlatColumn {
    std::vector<std::string> values;
    std::string scalar;

    void Book(TTree* tree, const std::string& name, bool asScalar) override
    {
        if (asScalar) tree->Branch(name.c_str(), &scalar);
        else          tree->Branch(name.c_str(), &values);
    }

    void Clear() override { values.clear(); }

    void Append(const void* address, EDataType, int kind) override
    { values.push_back(ReadString(address, kind)); }

    void Set(const void* address, EDataType, int kind) override
    { scalar = ReadString(address, kind); }

    std::string TypeName(bool asScalar) const override
    { return asScalar ? "string" : "vector<string>"; }
};

std::unique_ptr<FlatColumn> MakeColumn(EColumnType type)
{
    switch (type) {
        case kBoolColumn:   return std::unique_ptr<FlatColumn>(new NumericColumn<Bool_t>());
        case kIntColumn:    return std::unique_ptr<FlatColumn>(new NumericColumn<Int_t>());
        case kLongColumn:   return std::unique_ptr<FlatColumn>(new NumericColumn<Long64_t>());
        case kFloatColumn:  return std::unique_ptr<FlatColumn>(new NumericColumn<Float_t>());
        case kDoubleColumn: return std::unique_ptr<FlatColumn>(new NumericColumn<Double_t>());
        default:            return std::unique_ptr<FlatColumn>(new StringColumn());
    }
}

} // anonymous namespace

/* -------------------------------------------------------------------------- */
/*                              Fields and groups                             */
/* -------------------------------------------------------------------------- */

struct TreeFlattener::Field {
    std::string name;
    Long_t offset = 0;
    EDataType type = kNoType_t;
    int kind = kNumericMember;
    FlatColumn* column = nullptr;   // owned by the group
};

struct TreeFlattener::Group {
    std::string prefix;
    bool scalar = false;        // one value per entry rather than one per row
    int parent = -1;
    Long_t offset = 0;          // from an element of the parent to this collection

    TClass* collectionClass = nullptr;
    TVirtualCollectionProxy* proxy = nullptr;   // owned
    TClass* valueClass = nullptr;               // null when the elements are numbers
    EDataType valueType = kNoType_t;

    // Where this came from, carried into the schema
    std::string origin;
    std::string originType;

    std::vector<Field> fields;
    std::vector<std::unique_ptr<FlatColumn>> columns;

    std::string parentColumnName;
    std::unique_ptr<FlatColumn> parentColumn;

    std::vector<int> children;
    Int_t rows = 0;

    ~Group() { delete proxy; }
};

/* -------------------------------------------------------------------------- */
/*                             Construction and setup                         */
/* -------------------------------------------------------------------------- */

TreeFlattener::TreeFlattener(const FlattenOptions& options) : fOptions(options) {}

TreeFlattener::~TreeFlattener()
{
    for (auto& branch : fBranchObjects) {
        if (branch->object && branch->cls) branch->cls->Destructor(branch->object);
    }
}

size_t TreeFlattener::NGroups() const { return fGroups.size(); }

size_t TreeFlattener::NColumns() const
{
    size_t n = fPlainBranches.size();
    for (const auto& group : fGroups) {
        n += group->fields.size();
        if (group->parentColumn) n++;
    }
    return n;
}

std::string TreeFlattener::Rename(const std::string& name) const
{
    const auto it = fOptions.rename.find(name);
    return (it == fOptions.rename.end()) ? name : it->second;
}

bool TreeFlattener::Selected(const char* branchName) const
{
    const std::string name = branchName;

    if (std::find(fOptions.exclude.begin(), fOptions.exclude.end(), name)
        != fOptions.exclude.end()) {
        return false;
    }
    if (fOptions.include.empty()) return true;
    return std::find(fOptions.include.begin(), fOptions.include.end(), name)
           != fOptions.include.end();
}

std::string TreeFlattener::Claim(std::string name)
{
    if (std::find(fColumnNames.begin(), fColumnNames.end(), name) == fColumnNames.end()) {
        fColumnNames.push_back(name);
        return name;
    }

    // Two different members flattening onto the same name would silently
    // overwrite each other, so make it unique and say so
    int suffix = 2;
    std::string candidate;
    do {
        candidate = name + "_" + std::to_string(suffix++);
    } while (std::find(fColumnNames.begin(), fColumnNames.end(), candidate)
             != fColumnNames.end());

    std::cerr << "Warning: the column name '" << name << "' is already taken; "
              << "writing this one as '" << candidate << "'." << std::endl;
    fColumnNames.push_back(candidate);
    return candidate;
}

int TreeFlattener::NewGroup(int parent, const std::string& prefix, bool scalar)
{
    fGroups.emplace_back(new Group());
    Group& group = *fGroups.back();
    group.prefix = prefix;
    group.scalar = scalar;
    group.parent = parent;

    if (parent >= 0) {
        group.origin = fGroups[parent]->origin;
        group.originType = fGroups[parent]->originType;
    }

    return static_cast<int>(fGroups.size()) - 1;
}

void TreeFlattener::AddField(int group, const std::string& name, Long_t offset,
                             EDataType type, int kind)
{
    Group& target = *fGroups[group];

    Field field;
    field.name = Claim(name);
    field.offset = offset;
    field.type = type;
    field.kind = kind;

    target.columns.push_back(MakeColumn(ColumnTypeFor(type, kind)));
    field.column = target.columns.back().get();
    target.fields.push_back(field);
}

/* -------------------------------------------------------------------------- */
/*                                Layout discovery                            */
/* -------------------------------------------------------------------------- */

bool TreeFlattener::Connect(TTree* input)
{
    fInput = input;
    if (!input) return false;

    input->SetBranchStatus("*", 1);

    TObjArray* branches = input->GetListOfBranches();
    for (int i = 0; i < branches->GetEntriesFast(); ++i) {
        TBranch* branch = dynamic_cast<TBranch*>(branches->At(i));
        if (!branch || !Selected(branch->GetName())) continue;
        AddBranch(input, branch->GetName());
    }

    return NColumns() > 0;
}

void TreeFlattener::AddBranch(TTree* input, const char* name)
{
    const std::string type = ProductSchema::BranchType(input, name);
    if (type.empty()) {
        std::cerr << "Warning: could not work out the type of branch '" << name
                  << "'; skipping it." << std::endl;
        return;
    }

    if (TClass* cls = TClass::GetClass(type.c_str())) {
        AddObjectBranch(input, name, cls);
        return;
    }

    if (TDataType* dataType = gROOT->GetType(type.c_str())) {
        AddPlainBranch(input, name, static_cast<EDataType>(dataType->GetType()));
        return;
    }

    std::cerr << "Warning: branch '" << name << "' has type '" << type
              << "', which has no dictionary; skipping it." << std::endl;
}

void TreeFlattener::AddPlainBranch(TTree* input, const char* name, EDataType type)
{
    fPlainBranches.emplace_back(new PlainBranch());
    PlainBranch& plain = *fPlainBranches.back();
    plain.name = Claim(Rename(name));
    plain.type = type;
    plain.column = MakeColumn(ColumnTypeFor(type, kNumericMember));

    // Read into the widest buffer we have and convert on the way out. The
    // explicit overload is needed because the buffer's own type is not the
    // branch's, and the deducing overloads would rightly object.
    input->SetBranchAddress(name, &plain.buffer, nullptr, type, kFALSE);
}

void TreeFlattener::AddObjectBranch(TTree* input, const char* name, TClass* cls)
{
    fBranchObjects.emplace_back(new BranchObject());
    BranchObject& branch = *fBranchObjects.back();
    branch.name = name;
    branch.cls = cls;
    branch.object = cls->New();

    if (!branch.object) {
        std::cerr << "Warning: could not create a '" << cls->GetName()
                  << "' to read branch '" << name << "'; skipping it." << std::endl;
        fBranchObjects.pop_back();
        return;
    }

    // The generic overload: a void* alone does not tell ROOT what it points
    // at, so the class and "this is a pointer to a pointer" are passed too
    input->SetBranchAddress(name, &branch.object, cls, kOther_t, kTRUE);

    if (cls->GetCollectionProxy()) {
        // A collection branch: one row per element, named after the branch
        branch.group = AddCollection(-1, Rename(name), 0, cls, 1);
    } else {
        // A single object per entry: its numbers are per-event scalars, and
        // the branch name is dropped from their column names
        const std::string prefix = fOptions.keepBranchPrefix ? Rename(name) : std::string();
        branch.group = NewGroup(-1, prefix, true);
        AddMembers(branch.group, cls, 0, prefix, 1);
    }

    if (branch.group >= 0) {
        fGroups[branch.group]->origin = name;
        fGroups[branch.group]->originType = cls->GetName();
        // Children created above inherited an empty origin, so fix them up
        for (size_t i = static_cast<size_t>(branch.group); i < fGroups.size(); ++i) {
            if (fGroups[i]->origin.empty()) {
                fGroups[i]->origin = name;
                fGroups[i]->originType = cls->GetName();
            }
        }
    }
}

void TreeFlattener::AddMembers(int group, TClass* cls, Long_t offset,
                               const std::string& prefix, int depth)
{
    if (!cls) return;

    // Base classes first, so that inherited members keep their usual order.
    // TObject's bookkeeping members are not data and are skipped.
    if (TList* bases = cls->GetListOfBases()) {
        TIter next(bases);
        while (TBaseClass* base = dynamic_cast<TBaseClass*>(next())) {
            TClass* baseClass = base->GetClassPointer();
            if (!baseClass) continue;
            if (std::strcmp(baseClass->GetName(), "TObject") == 0) continue;
            AddMembers(group, baseClass, offset + base->GetDelta(), prefix, depth);
        }
    }

    TList* members = cls->GetListOfDataMembers();
    if (!members) return;

    TIter next(members);
    while (TDataMember* member = dynamic_cast<TDataMember*>(next())) {

        if (member->Property() & kIsStatic) continue;
        if (!member->IsPersistent()) continue;   // marked //! in the header

        const std::string name = member->GetName();
        const std::string full = prefix.empty() ? name : prefix + "_" + name;
        const Long_t memberOffset = offset + member->GetOffset();

        if (member->IsaPointer()) {
            std::cerr << "Warning: skipping '" << cls->GetName() << "::" << name
                      << "', a pointer member the flattener cannot follow."
                      << std::endl;
            continue;
        }
        if (member->GetArrayDim() > 0) {
            std::cerr << "Warning: skipping '" << cls->GetName() << "::" << name
                      << "', a fixed-size array." << std::endl;
            continue;
        }

        // A number
        TDataType* dataType = member->GetDataType();
        if (dataType && dataType->GetType() != kOther_t
                     && dataType->GetType() != kNoType_t) {
            AddField(group, full, memberOffset,
                     static_cast<EDataType>(dataType->GetType()), kNumericMember);
            continue;
        }

        TClass* memberClass = TClass::GetClass(member->GetTypeName());
        if (!memberClass) {
            std::cerr << "Warning: skipping '" << cls->GetName() << "::" << name
                      << "' of type '" << member->GetTypeName()
                      << "', which has no dictionary." << std::endl;
            continue;
        }

        // A string
        if (IsTString(memberClass)) {
            AddField(group, full, memberOffset, kNoType_t, kTStringMember);
            continue;
        }
        if (IsStdString(memberClass)) {
            AddField(group, full, memberOffset, kNoType_t, kStdStringMember);
            continue;
        }

        // A collection: a group of its own, one row per element
        if (memberClass->GetCollectionProxy()) {
            if (depth >= fOptions.maxDepth) {
                if (fOptions.verbose) {
                    std::cout << "   (not following " << cls->GetName() << "::" << name
                              << ", maxDepth " << fOptions.maxDepth << " reached)"
                              << std::endl;
                }
                continue;
            }
            AddCollection(group, Rename(full), memberOffset, memberClass, depth + 1);
            continue;
        }

        // Any other class: walk into it, prefixing its members
        AddMembers(group, memberClass, memberOffset, full, depth);
    }
}

int TreeFlattener::AddCollection(int parent, const std::string& prefix, Long_t offset,
                                 TClass* collectionClass, int depth)
{
    TVirtualCollectionProxy* prototype = collectionClass->GetCollectionProxy();
    if (!prototype) return -1;

    const int group = NewGroup(parent, prefix, false);
    Group& target = *fGroups[group];
    target.offset = offset;
    target.collectionClass = collectionClass;
    target.proxy = prototype->Generate();   // our own, so nesting is safe
    target.valueClass = prototype->GetValueClass();
    target.valueType = static_cast<EDataType>(prototype->GetType());

    if (parent >= 0) {
        fGroups[parent]->children.push_back(group);

        // An index into the parent is only informative when the parent has
        // more than one row; at the top of an entry there is only ever one
        if (!fGroups[parent]->scalar) {
            target.parentColumnName = Claim(prefix + "_parent");
            target.parentColumn = MakeColumn(kIntColumn);
        }
    }

    if (IsTString(target.valueClass)) {
        AddField(group, prefix, 0, kNoType_t, kTStringMember);
    } else if (IsStdString(target.valueClass)) {
        AddField(group, prefix, 0, kNoType_t, kStdStringMember);
    } else if (target.valueClass) {
        AddMembers(group, target.valueClass, 0, prefix, depth);
    } else if (target.valueType != kNoType_t && target.valueType != kOther_t) {
        // A collection of numbers: one column, named after the member itself
        AddField(group, prefix, 0, target.valueType, kNumericMember);
    } else {
        std::cerr << "Warning: '" << collectionClass->GetName()
                  << "' holds elements the flattener does not understand; "
                  << "'" << prefix << "' will be empty." << std::endl;
    }

    return group;
}

/* -------------------------------------------------------------------------- */
/*                                    Output                                  */
/* -------------------------------------------------------------------------- */

void TreeFlattener::Book(TTree* output)
{
    if (!output) return;

    for (auto& plain : fPlainBranches) {
        plain->column->Book(output, plain->name, true);
    }

    for (auto& group : fGroups) {
        if (group->parentColumn) {
            group->parentColumn->Book(output, group->parentColumnName, false);
        }
        for (const Field& field : group->fields) {
            field.column->Book(output, field.name, group->scalar);
        }
    }
}

std::vector<ProductInfo> TreeFlattener::Schema(const std::string& outputTreeName) const
{
    std::vector<ProductInfo> products;

    for (const auto& plain : fPlainBranches) {
        ProductInfo product;
        product.tree = outputTreeName;
        product.branch = plain->name;
        product.type = plain->column->TypeName(true);
        products.push_back(product);
    }

    for (const auto& group : fGroups) {
        if (group->parentColumn) {
            ProductInfo product;
            product.tree = outputTreeName;
            product.branch = group->parentColumnName;
            product.type = group->parentColumn->TypeName(false);
            product.producer = group->origin;
            product.producerType = group->originType;
            products.push_back(product);
        }
        for (const Field& field : group->fields) {
            ProductInfo product;
            product.tree = outputTreeName;
            product.branch = field.name;
            product.type = field.column->TypeName(group->scalar);
            product.producer = group->origin;
            product.producerType = group->originType;
            products.push_back(product);
        }
    }

    return products;
}

void TreeFlattener::PrintLayout(std::ostream& out) const
{
    out << "Flattening '" << (fInput ? fInput->GetName() : "?") << "' into "
        << NColumns() << " column(s):\n";

    if (!fPlainBranches.empty()) {
        out << "    copied through:";
        for (const auto& plain : fPlainBranches) out << " " << plain->name;
        out << "\n";
    }

    for (const auto& group : fGroups) {
        if (group->fields.empty()) continue;
        out << "    " << (group->scalar ? "per event" : "per row")
            << "  [" << (group->prefix.empty() ? "(branch)" : group->prefix) << "]"
            << "  from " << group->originType;
        if (!group->parentColumnName.empty()) {
            out << ", indexed by " << group->parentColumnName;
        }
        out << "\n";
        for (const Field& field : group->fields) {
            out << "        " << field.name << "  "
                << field.column->TypeName(group->scalar) << "\n";
        }
    }
    out << std::flush;
}

/* -------------------------------------------------------------------------- */
/*                                   Filling                                  */
/* -------------------------------------------------------------------------- */

Long64_t TreeFlattener::FillEntry(Long64_t entry)
{
    if (!fInput) return -1;

    for (auto& group : fGroups) {
        group->rows = 0;
        for (auto& column : group->columns) column->Clear();
        if (group->parentColumn) group->parentColumn->Clear();
    }

    const Long64_t read = fInput->GetEntry(entry);
    if (read < 0) return read;

    for (auto& plain : fPlainBranches) {
        plain->column->Set(&plain->buffer, plain->type, kNumericMember);
    }

    for (auto& branch : fBranchObjects) {
        if (branch->group < 0 || !branch->object) continue;
        const char* object = static_cast<const char*>(branch->object);
        if (fGroups[branch->group]->scalar) {
            FillObject(branch->group, object, 0);
        } else {
            FillCollection(branch->group, object, 0);
        }
    }

    return read;
}

void TreeFlattener::FillCollection(int group, const char* collection, Int_t parentRow)
{
    Group& target = *fGroups[group];
    if (!target.proxy || !collection) return;

    TVirtualCollectionProxy::TPushPop guard(target.proxy,
                                            const_cast<char*>(collection));
    const UInt_t size = target.proxy->Size();
    for (UInt_t i = 0; i < size; ++i) {
        FillObject(group, static_cast<const char*>(target.proxy->At(i)), parentRow);
    }
}

void TreeFlattener::FillObject(int group, const char* object, Int_t parentRow)
{
    Group& target = *fGroups[group];
    if (!object) return;

    const Int_t row = target.scalar ? 0 : target.rows++;

    if (target.parentColumn) {
        target.parentColumn->Append(&parentRow, kInt_t, kNumericMember);
    }

    for (const Field& field : target.fields) {
        const void* address = object + field.offset;
        if (target.scalar) field.column->Set(address, field.type, field.kind);
        else               field.column->Append(address, field.type, field.kind);
    }

    for (const int child : target.children) {
        FillCollection(child, object + fGroups[child]->offset, row);
    }
}

} // namespace fastgarsim
