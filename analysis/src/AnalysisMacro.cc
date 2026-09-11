 /***************************************************************************
 * AnalysisMacro.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the job macro parser and of the run-time compilation
 *   of an analysis macro.
 *
 ***************************************************************************/

#include "AnalysisMacro.hh"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <system_error>

#include "TInterpreter.h"
#include "TROOT.h"
#include "TSystem.h"

namespace ana {

namespace fs = std::filesystem;

namespace {

const std::string kGlobalPrefix = "/ana/global/";
const std::string kAnaPrefix    = "/ana/";

std::string Trim(const std::string& text)
{
    const size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";

    const size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

// A '#' starts a comment, as in the reconstruction macros, so a value never
// contains one
std::string StripComment(const std::string& line)
{
    const size_t hash = line.find('#');
    return hash == std::string::npos ? line : line.substr(0, hash);
}

}  // namespace

/* -------------------------------------------------------------------------- */
/*                                 Job macro                                  */
/* -------------------------------------------------------------------------- */

const std::vector<std::string>& JobMacroGlobals()
{
    // Order follows the README's table rather than the alphabet: what a job
    // usually sets first
    static const std::vector<std::string> names = {
        "analysis", "input", "output", "genie", "outputTree",
        "maxEvents", "firstEvent", "nReports",
        "requireSim", "printProducts", "printGeometry", "xrootdForPnfs",
        "simTreeName", "recoTreeName", "geoTreeName", "genieTreeName",
        "simBranchName"
    };
    return names;
}

Bool_t ParseJobMacro(const std::string& path, AnalysisConfig& config,
                     std::string* analysisMacro)
{
    std::ifstream macro(path);
    if (!macro.is_open()) {
        std::cerr << "Error: cannot open the job macro '" << path << "'" << std::endl;
        return kFALSE;
    }

    // The global commands are collected first and applied together, so that
    // they are converted by the same code that converts the analysis's own
    // parameters, with the same reporting of a value of the wrong type
    ParameterSet globals;
    Bool_t ok = kTRUE;
    std::string line;

    for (Int_t number = 1; std::getline(macro, line); ++number) {
        const std::string command = Trim(StripComment(line));
        if (command.empty()) continue;

        // Split into the command and everything after it: a value may contain
        // spaces, as a comma-separated list of input files does
        const size_t space = command.find_first_of(" \t");
        const std::string name = command.substr(0, space);
        const std::string value =
            space == std::string::npos ? "" : Trim(command.substr(space + 1));

        if (name.rfind(kAnaPrefix, 0) != 0) {
            std::cerr << "Error: " << path << ":" << number
                      << ": '" << name << "' is not an analysis command"
                      << " (they all start with " << kAnaPrefix << ")" << std::endl;
            ok = kFALSE;
            continue;
        }

        if (value.empty()) {
            std::cerr << "Error: " << path << ":" << number
                      << ": '" << name << "' was given no value" << std::endl;
            ok = kFALSE;
            continue;
        }

        if (name.rfind(kGlobalPrefix, 0) == 0) {
            globals.Set(name.substr(kGlobalPrefix.size()), value);
        } else {
            config.params.Set(name.substr(kAnaPrefix.size()), value);
        }
    }

    /* --------------------------- Apply the globals ------------------------ */

    // Read even when the caller does not want it, so that a job macro naming
    // its analysis is not reported below as an unknown global
    std::string ignored;
    globals.Get("analysis", analysisMacro ? *analysisMacro : ignored);

    globals.Get("input",         config.inputFiles);
    globals.Get("output",        config.outputFile);
    globals.Get("genie",         config.genieFiles);
    globals.Get("outputTree",    config.outputTreeName);
    globals.Get("maxEvents",     config.maxEvents);
    globals.Get("firstEvent",    config.firstEvent);
    globals.Get("nReports",      config.nReports);
    globals.Get("requireSim",    config.requireSim);
    globals.Get("printProducts", config.printProducts);
    globals.Get("printGeometry", config.printGeometry);
    globals.Get("xrootdForPnfs", config.xrootdForPnfs);
    globals.Get("simTreeName",   config.simTreeName);
    globals.Get("recoTreeName",  config.recoTreeName);
    globals.Get("geoTreeName",   config.geoTreeName);
    globals.Get("genieTreeName", config.genieTreeName);
    globals.Get("simBranchName", config.simBranchName);

    // Every name above has now been asked for, so anything left over is a
    // global that does not exist
    for (const std::string& unknown : globals.UnusedKeys()) {
        std::cerr << "Error: " << path << ": there is no " << kGlobalPrefix << unknown
                  << "\n   known names:";
        for (const std::string& known : JobMacroGlobals()) std::cerr << " " << known;
        std::cerr << std::endl;
        ok = kFALSE;
    }

    if (globals.NErrors() > 0) ok = kFALSE;

    return ok;
}

/* -------------------------------------------------------------------------- */
/*                              Analysis macro                                */
/* -------------------------------------------------------------------------- */

namespace {

// The private ACLiC build directory of this process, removed when the process
// ends however it ends. Nothing is written next to the macro, so an analysis
// can live in a read-only directory -- an installed or unpacked build, or a
// shared production area -- and, since the directory starts out empty, ACLiC
// has nothing to reuse and compiles the macro afresh on every run. That is
// what makes an edit between two runs take effect without any bookkeeping.
struct BuildDir {
    std::string path;

    ~BuildDir() { Remove(); }

    Bool_t Create()
    {
        std::error_code ec;
        const fs::path directory =
            fs::temp_directory_path(ec) /
            ("GArAnalysis-" + std::to_string(gSystem->GetPid()));

        if (ec) {
            std::cerr << "Error: no temporary directory to compile in: "
                      << ec.message() << std::endl;
            return kFALSE;
        }

        // A directory left behind by a killed job of the same process number
        // would hand this one a stale library, so start from nothing
        fs::remove_all(directory, ec);
        fs::create_directories(directory, ec);
        if (ec) {
            std::cerr << "Error: cannot create the build directory '"
                      << directory.string() << "': " << ec.message() << std::endl;
            return kFALSE;
        }

        path = directory.string();
        return kTRUE;
    }

    void Remove()
    {
        if (path.empty()) return;

        // The library compiled here is still loaded. Unlinking an open shared
        // object is safe on Linux and macOS -- it stays mapped until the
        // process ends -- which is what lets the job clean up after itself
        // while it is still running.
        std::error_code ec;
        fs::remove_all(path, ec);
        path.clear();
    }
};

BuildDir& TheBuildDir()
{
    static BuildDir directory;
    return directory;
}

}  // namespace

void RemoveBuildDir()
{
    TheBuildDir().Remove();
}

AnalysisBase* CompileAnalysis(const std::string& path,
                              const std::vector<std::string>& includeDirs,
                              const std::vector<std::string>& libraryDirs)
{
    std::error_code ec;
    if (!fs::is_regular_file(path, ec)) {
        std::cerr << "Error: no analysis macro at '" << path << "'" << std::endl;
        return nullptr;
    }

    // Headers first: the analysis macro includes AnalysisBase.hh and the
    // reconstruction data types, and ACLiC has to find them whether or not
    // setup.sh was sourced
    for (const std::string& directory : includeDirs) {
        if (directory.empty() || !fs::is_directory(directory, ec)) continue;
        gSystem->AddIncludePath((" -I\"" + directory + "\"").c_str());
    }

    // Then the libraries, so that the dictionaries of the simulation and
    // reconstruction data types autoload from their .rootmap files. The
    // analysis library itself is not loaded here: this code is part of it, so
    // it is already in the process, and the compiled macro resolves against
    // it when the macro's library is loaded.
    for (const std::string& directory : libraryDirs) {
        if (directory.empty() || !fs::is_directory(directory, ec)) continue;
        gSystem->AddDynamicPath(directory.c_str());
    }

    if (!TheBuildDir().Create()) return nullptr;
    gSystem->SetBuildDir(TheBuildDir().path.c_str(), kTRUE);

    std::cout << "\n Compiling " << path << " ..." << std::endl;

    const auto started = std::chrono::steady_clock::now();

    // '++' forces the compilation, so the macro is rebuilt even if anything
    // else has already produced a library for it
    Int_t error = 0;
    const Int_t loaded = gROOT->LoadMacro((path + "++").c_str(), &error);

    const Double_t seconds =
        std::chrono::duration<Double_t>(std::chrono::steady_clock::now() - started).count();

    if (loaded != 0) {
        std::cerr << "\nError: " << path << " did not compile." << std::endl;
        std::cerr << "   If the errors above are missing headers, the analysis"
                  << " headers are not on ROOT's\n   include path; source the"
                  << " build's setup.sh, or run GArAnalysis from the build"
                  << " directory." << std::endl;
        return nullptr;
    }

    std::cout << "   compiled in " << std::fixed << std::setprecision(1)
              << seconds << " s" << std::defaultfloat << std::endl;

    // ANA_ANALYSIS() in the macro defines this function; calling it is what
    // creates the analysis object
    TInterpreter::EErrorCode status = TInterpreter::kNoError;
    const auto address = gInterpreter->Calc("ana_MakeAnalysis()", &status);

    if (status != TInterpreter::kNoError || address == 0) {
        std::cerr << "\nError: " << path << " compiled but declares no analysis."
                  << "\n   An analysis macro ends with the line"
                  << "\n\n      ANA_ANALYSIS(MyAnalysisClass)\n"
                  << "\n   naming its ana::AnalysisBase subclass, which has to be"
                  << " default-constructible.\n" << std::endl;
        return nullptr;
    }

    return reinterpret_cast<AnalysisBase*>(address);
}

}  // namespace ana
