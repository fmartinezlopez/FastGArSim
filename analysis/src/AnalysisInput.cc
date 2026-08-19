 /***************************************************************************
 * AnalysisInput.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the input-specification handling.
 *
 ***************************************************************************/

#include "AnalysisInput.hh"

#include <cstdlib>
#include <iostream>
#include <set>

#include <glob.h>

#include "TString.h"
#include "TSystem.h"

namespace ana {

const char* const kDefaultXRootDDoor = "root://fndca1.fnal.gov:1094";
const char* const kDefaultPnfsPrefix = "/pnfs/fnal.gov/usr/";

namespace {

const char* const kPnfsRoot = "/pnfs/";

std::string EnvOr(const char* variable, const char* fallback)
{
    const char* value = std::getenv(variable);
    return (value && *value) ? std::string(value) : std::string(fallback);
}

std::string Trim(const std::string& s)
{
    const size_t first = s.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return std::string();
    const size_t last = s.find_last_not_of(" \t\n\r");
    return s.substr(first, last - first + 1);
}

std::vector<std::string> Split(const std::string& s, char delimiter)
{
    std::vector<std::string> parts;
    size_t start = 0;
    while (true) {
        const size_t pos = s.find(delimiter, start);
        if (pos == std::string::npos) {
            parts.push_back(s.substr(start));
            break;
        }
        parts.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return parts;
}

Bool_t HasWildcard(const std::string& s)
{
    return s.find_first_of("*?[") != std::string::npos;
}

// Expand ~ and $VARIABLES, leaving wildcards alone
std::string ExpandShellVariables(const std::string& path)
{
    TString expanded(path.c_str());
    // Returns kTRUE on failure, in which case the original is kept
    if (gSystem->ExpandPathName(expanded)) return path;
    return std::string(expanded.Data());
}

// Shell-style expansion of a single pattern against the local file system
std::vector<std::string> Glob(const std::string& pattern)
{
    std::vector<std::string> matches;

    glob_t results;
    const int status = ::glob(pattern.c_str(), 0, nullptr, &results);

    if (status == 0) {
        matches.reserve(results.gl_pathc);
        for (size_t i = 0; i < results.gl_pathc; ++i) {
            matches.emplace_back(results.gl_pathv[i]);
        }
    } else if (status != GLOB_NOMATCH) {
        std::cerr << "Warning: could not expand the pattern '" << pattern
                  << "' (glob returned " << status << ")" << std::endl;
    }

    ::globfree(&results);
    return matches;
}

} // anonymous namespace

/* -------------------------------------------------------------------------- */
/*                              XRootD settings                               */
/* -------------------------------------------------------------------------- */

std::string XRootDDoor() { return EnvOr("FASTGAR_XROOTD_DOOR", kDefaultXRootDDoor); }
std::string PnfsPrefix() { return EnvOr("FASTGAR_PNFS_PREFIX", kDefaultPnfsPrefix); }

std::string XRootDHost()
{
    const std::string door = XRootDDoor();
    const size_t scheme = door.find("://");
    return (scheme == std::string::npos) ? door : door.substr(scheme + 3);
}

/* -------------------------------------------------------------------------- */
/*                                Path helpers                                */
/* -------------------------------------------------------------------------- */

Bool_t IsPnfsPath(const std::string& path)
{
    return path.compare(0, std::string(kPnfsRoot).size(), kPnfsRoot) == 0;
}

Bool_t IsRemoteURL(const std::string& path)
{
    const size_t scheme = path.find("://");
    return scheme != std::string::npos && scheme > 0;
}

std::string ToXRootD(const std::string& path)
{
    if (IsRemoteURL(path) || !IsPnfsPath(path)) return path;

    // Strip the leading "/pnfs/" and graft on the dCache prefix
    const std::string rest = path.substr(std::string(kPnfsRoot).size());
    return XRootDDoor() + PnfsPrefix() + rest;
}

/* -------------------------------------------------------------------------- */
/*                             Input specification                            */
/* -------------------------------------------------------------------------- */

std::vector<std::string> ExpandInputFiles(const std::string& specification,
                                          Bool_t xrootdForPnfs)
{
    std::vector<std::string> files;
    std::set<std::string> seen;

    // Reading the same file twice would silently double the statistics
    auto append = [&files, &seen](const std::string& path) {
        if (seen.insert(path).second) files.push_back(path);
    };

    for (const std::string& rawEntry : Split(specification, ',')) {

        const std::string entry = Trim(rawEntry);
        if (entry.empty()) continue;

        // A remote path cannot be listed, so it is taken at face value
        if (IsRemoteURL(entry)) {
            if (HasWildcard(entry)) {
                std::cerr << "Warning: '" << entry << "' is a URL, so its wildcards "
                          << "cannot be expanded; passing it through unchanged."
                          << std::endl;
            }
            append(entry);
            continue;
        }

        const std::string path = ExpandShellVariables(entry);

        if (HasWildcard(path)) {
            const std::vector<std::string> matches = Glob(path);
            if (matches.empty()) {
                std::cerr << "Warning: the pattern '" << path << "' matched no files"
                          << std::endl;
            }
            for (const std::string& match : matches) append(match);
            continue;
        }

        // A named /pnfs file will be fetched over XRootD, so it need not be
        // visible through the local mount -- the point of using XRootD is
        // often that it is not. Only wildcards genuinely require the mount,
        // since a remote area cannot be listed.
        if (xrootdForPnfs && IsPnfsPath(path)) {
            append(path);
            continue;
        }

        // AccessPathName is inverted: it returns kFALSE when the path exists
        if (gSystem->AccessPathName(path.c_str())) {
            std::cerr << "Warning: '" << path << "' does not exist; skipping it."
                      << std::endl;
            continue;
        }
        append(path);
    }

    if (files.empty()) {
        std::cerr << "Error: the input specification '" << specification
                  << "' resolved to no files." << std::endl;
        return files;
    }

    // Stream anything living in dCache rather than reading it through the mount
    if (xrootdForPnfs) {
        Int_t nRewritten = 0;
        for (std::string& file : files) {
            const std::string url = ToXRootD(file);
            if (url != file) {
                file = url;
                nRewritten++;
            }
        }
        if (nRewritten > 0) {
            std::cout << "Reading " << nRewritten << " /pnfs file(s) over XRootD via "
                      << XRootDDoor() << std::endl;
        }
    }

    return files;
}

/* -------------------------------------------------------------------------- */
/*                            Failure diagnostics                             */
/* -------------------------------------------------------------------------- */

void ExplainOpenFailure(const std::string& path)
{
    std::cerr << "\nCould not open the input:\n    " << path << "\n" << std::endl;

    if (!IsRemoteURL(path)) {
        std::cerr << "The file exists but ROOT could not read the tree from it. It may "
                     "be\ntruncated or still being written, or the tree name may be "
                     "wrong -- check\nit with:\n\n"
                     "    root -l " << path << "  then  .ls\n" << std::endl;
        return;
    }

    // ROOT has already printed its own diagnosis just above, so point at the
    // most common cause rather than repeating it
    std::cerr <<
        "This is a remote file. If the error above says 'Auth failed' or 'No protocols\n"
        "left to try', the problem is missing or expired credentials, not the file.\n"
        "\n"
        "Get a token and a proxy:\n"
        "\n"
        "    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh\n"
        "    setup_fnal_security\n"
        "\n"
        "or individually:\n"
        "\n"
        "    htgettoken -a htvaultprod.fnal.gov -i dune\n"
        "    voms-proxy-init -rfc -noregen -voms dune:/dune/Role=Analysis\n"
        "\n"
        "Check what you have with 'voms-proxy-info -all' and 'httokendecode', and test\n"
        "the door directly, which isolates the problem from ROOT:\n"
        "\n"
        "    xrdfs " << XRootDHost() << " ls /pnfs/fnal.gov/usr/dune/...\n"
        "\n"
        "A different door or site can be selected with FASTGAR_XROOTD_DOOR and\n"
        "FASTGAR_PNFS_PREFIX.\n"
        << std::endl;
}

} // namespace ana
