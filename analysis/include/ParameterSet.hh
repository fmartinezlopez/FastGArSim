 /***************************************************************************
 * ParameterSet.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   The parameters an analysis was configured with, as they were read from
 *   the job macro: a set of name/value pairs that convert themselves to
 *   whatever type the analysis asks for.
 *
 *   An analysis reads them in Configure(), once, before the job starts:
 *
 *     void Configure(const ana::ParameterSet& p) override {
 *         p.Get("truncation", fTruncation);
 *         p.Get("minHits", fMinHits);
 *     }
 *
 *   A name the macro did not set leaves the member alone, so the value it
 *   was declared with is its default -- there is no second place where the
 *   defaults have to be repeated.
 *
 *   Since the names are strings, nothing checks them at compile time, so
 *   both kinds of mistake are caught here instead: a value that cannot be
 *   read as the type asked for counts an error and stops the job, and a name
 *   nothing ever asked for is reported at start-up, which is what catches a
 *   misspelling in the macro.
 *
 ***************************************************************************/

#ifndef ParameterSet_hh
#define ParameterSet_hh

#include <cstddef>
#include <iosfwd>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Rtypes.h"

namespace ana {

class ParameterSet {
public:
    /* ------------------------------- Filling ------------------------------ */

    // Set one parameter, replacing whatever was there before
    void Set(const std::string& name, const std::string& value);

    Bool_t Has(const std::string& name) const;
    size_t Size() const { return fValues.size(); }
    Bool_t IsEmpty() const { return fValues.empty(); }

    /* ------------------------------- Reading ------------------------------ */

    // Read one parameter into `value`.
    //
    //   - the name is there and converts   -- `value` is overwritten, kTRUE
    //   - the name is not there            -- `value` is untouched, kFALSE
    //   - the name is there but the text
    //     cannot be read as this type      -- `value` is untouched, kFALSE,
    //                                         a message, and NErrors() rises
    //
    // Asking for a name is what marks it as used, which is how UnusedKeys()
    // knows what the analysis never looked at.
    Bool_t Get(const std::string& name, std::string& value) const;
    Bool_t Get(const std::string& name, Bool_t& value) const;
    Bool_t Get(const std::string& name, Int_t& value) const;
    Bool_t Get(const std::string& name, Long64_t& value) const;
    Bool_t Get(const std::string& name, Float_t& value) const;
    Bool_t Get(const std::string& name, Double_t& value) const;

    // The same, for a parameter used in an expression rather than stored in a
    // member. `fallback` comes back when the name is absent or unreadable.
    std::string GetString(const std::string& name, const std::string& fallback = "") const;
    Bool_t GetBool(const std::string& name, Bool_t fallback) const;
    Int_t GetInt(const std::string& name, Int_t fallback) const;
    Long64_t GetLong(const std::string& name, Long64_t fallback) const;
    Double_t GetDouble(const std::string& name, Double_t fallback) const;

    // Raw text of a parameter, without conversion and without marking it used
    const std::string& Raw(const std::string& name) const;

    /* ------------------------------ Diagnostics --------------------------- */

    // Values that could not be read as the type the analysis asked for
    Int_t NErrors() const { return fErrors; }

    // Names nothing ever asked for -- in practice, misspellings in the macro
    std::vector<std::string> UnusedKeys() const;

    const std::map<std::string, std::string>& Values() const { return fValues; }

    // One "name  value" line per parameter, in name order
    void Print(std::ostream& out, const std::string& indent = "   ") const;

private:
    // Text of `name`, marked as used; nullptr when it is not there
    const std::string* Lookup(const std::string& name) const;

    std::map<std::string, std::string> fValues;

    // Bookkeeping over a const object: reading a parameter does not change
    // the configuration, only what is known about how it was used
    mutable std::set<std::string> fUsed;
    mutable Int_t fErrors = 0;
};

}  // namespace ana

#endif
