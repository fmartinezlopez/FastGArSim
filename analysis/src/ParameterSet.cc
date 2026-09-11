 /***************************************************************************
 * ParameterSet.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the analysis parameter store and its conversions.
 *
 ***************************************************************************/

#include "ParameterSet.hh"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <ostream>
#include <stdexcept>

namespace ana {

namespace {

const std::string kEmpty;

// Text left over after a conversion, ignoring trailing whitespace. Anything
// there means the value was only partly a number ("10cm", "0.5 0.7"), which
// is a mistake rather than something to take the leading digits of.
Bool_t FullyConsumed(const std::string& text, size_t used)
{
    for (size_t i = used; i < text.size(); ++i) {
        if (!std::isspace(static_cast<unsigned char>(text[i]))) return kFALSE;
    }
    return used > 0;
}

}  // namespace

/* -------------------------------------------------------------------------- */
/*                                   Filling                                  */
/* -------------------------------------------------------------------------- */

void ParameterSet::Set(const std::string& name, const std::string& value)
{
    fValues[name] = value;
}

Bool_t ParameterSet::Has(const std::string& name) const
{
    return fValues.find(name) != fValues.end();
}

const std::string* ParameterSet::Lookup(const std::string& name) const
{
    // Asking for a name counts as using it, whether or not it is there: an
    // analysis that asks for "minHits" and gets the default has still
    // accounted for that name
    fUsed.insert(name);

    const auto entry = fValues.find(name);
    return entry == fValues.end() ? nullptr : &entry->second;
}

const std::string& ParameterSet::Raw(const std::string& name) const
{
    const auto entry = fValues.find(name);
    return entry == fValues.end() ? kEmpty : entry->second;
}

/* -------------------------------------------------------------------------- */
/*                                  Reading                                   */
/* -------------------------------------------------------------------------- */

Bool_t ParameterSet::Get(const std::string& name, std::string& value) const
{
    const std::string* text = Lookup(name);
    if (!text) return kFALSE;

    value = *text;
    return kTRUE;
}

Bool_t ParameterSet::Get(const std::string& name, Bool_t& value) const
{
    const std::string* text = Lookup(name);
    if (!text) return kFALSE;

    std::string lowered = *text;
    std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lowered == "1" || lowered == "true"  || lowered == "yes" || lowered == "on"  ||
        lowered == "ktrue") {
        value = kTRUE;
        return kTRUE;
    }
    if (lowered == "0" || lowered == "false" || lowered == "no"  || lowered == "off" ||
        lowered == "kfalse") {
        value = kFALSE;
        return kTRUE;
    }

    std::cerr << "ParameterSet: '" << name << " " << *text
              << "' is not true or false" << std::endl;
    ++fErrors;
    return kFALSE;
}

Bool_t ParameterSet::Get(const std::string& name, Long64_t& value) const
{
    const std::string* text = Lookup(name);
    if (!text) return kFALSE;

    try {
        size_t used = 0;
        const long long parsed = std::stoll(*text, &used);
        if (!FullyConsumed(*text, used)) throw std::invalid_argument("trailing text");
        value = parsed;
        return kTRUE;
    } catch (const std::exception&) {
        std::cerr << "ParameterSet: '" << name << " " << *text
                  << "' is not a whole number" << std::endl;
        ++fErrors;
        return kFALSE;
    }
}

Bool_t ParameterSet::Get(const std::string& name, Int_t& value) const
{
    Long64_t wide = value;
    if (!Get(name, wide)) return kFALSE;

    value = static_cast<Int_t>(wide);
    return kTRUE;
}

Bool_t ParameterSet::Get(const std::string& name, Double_t& value) const
{
    const std::string* text = Lookup(name);
    if (!text) return kFALSE;

    try {
        size_t used = 0;
        const double parsed = std::stod(*text, &used);
        if (!FullyConsumed(*text, used)) throw std::invalid_argument("trailing text");
        value = parsed;
        return kTRUE;
    } catch (const std::exception&) {
        std::cerr << "ParameterSet: '" << name << " " << *text
                  << "' is not a number" << std::endl;
        ++fErrors;
        return kFALSE;
    }
}

Bool_t ParameterSet::Get(const std::string& name, Float_t& value) const
{
    Double_t wide = value;
    if (!Get(name, wide)) return kFALSE;

    value = static_cast<Float_t>(wide);
    return kTRUE;
}

/* -------------------------------------------------------------------------- */

std::string ParameterSet::GetString(const std::string& name, const std::string& fallback) const
{
    std::string value = fallback;
    Get(name, value);
    return value;
}

Bool_t ParameterSet::GetBool(const std::string& name, Bool_t fallback) const
{
    Bool_t value = fallback;
    Get(name, value);
    return value;
}

Int_t ParameterSet::GetInt(const std::string& name, Int_t fallback) const
{
    Int_t value = fallback;
    Get(name, value);
    return value;
}

Long64_t ParameterSet::GetLong(const std::string& name, Long64_t fallback) const
{
    Long64_t value = fallback;
    Get(name, value);
    return value;
}

Double_t ParameterSet::GetDouble(const std::string& name, Double_t fallback) const
{
    Double_t value = fallback;
    Get(name, value);
    return value;
}

/* -------------------------------------------------------------------------- */
/*                                Diagnostics                                 */
/* -------------------------------------------------------------------------- */

std::vector<std::string> ParameterSet::UnusedKeys() const
{
    std::vector<std::string> unused;

    for (const auto& entry : fValues) {
        if (fUsed.find(entry.first) == fUsed.end()) unused.push_back(entry.first);
    }
    return unused;
}

void ParameterSet::Print(std::ostream& out, const std::string& indent) const
{
    size_t widest = 0;
    for (const auto& entry : fValues) widest = std::max(widest, entry.first.size());

    for (const auto& entry : fValues) {
        out << indent << entry.first
            << std::string(widest - entry.first.size() + 2, ' ')
            << entry.second << "\n";
    }
}

}  // namespace ana
