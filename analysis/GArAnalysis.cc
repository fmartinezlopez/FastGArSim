 /***************************************************************************
 * GArAnalysis.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Single entry point for every FastGArSim analysis.
 *
 *   An analysis is a ROOT macro ending in ANA_ANALYSIS(), and a job is that
 *   macro plus the files to run it over and the parameters it takes:
 *
 *     GArAnalysis -a TruncatedDEDX.C -i 'gun_*.root' -o dedx.root -m dedx.mac
 *
 *   The macro is compiled when the job starts, not when the framework is
 *   built, so an analysis can be written, changed and run without rebuilding
 *   anything -- and a run always uses the file as it stands at that moment.
 *   Its parameters come from the job macro, the same kind of file that
 *   configures GArReconstruction; see AnalysisMacro.hh for its syntax.
 *
 ***************************************************************************/

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"

#include "MacroPath.hh"

#include "AnalysisBase.hh"
#include "AnalysisMacro.hh"

namespace {

void PrintUsage()
{
    std::cout << "\n Usage: GArAnalysis -a <analysis.C> -i <input> [options]\n"
              << "\n Options:"
              << "\n   -a <file.C>   Analysis macro to run; may also be given as the"
              << "\n                 first argument, or by the job macro"
              << "\n   -i <spec>     Input file(s): a path, a wildcard pattern or a"
              << "\n                 comma-separated list of either (required)"
              << "\n   -o <file>     Output ROOT file"
              << "\n   -m <file.mac> Job macro holding the analysis parameters"
              << "\n   -g <spec>     GENIE gst file(s) with the truth record"
              << "\n   -n <N>        Events to process; -1 (the default) is all of them"
              << "\n   -h, --help    Show this message"
              << "\n"
              << "\n The analysis and job macros are looked for in the working"
              << "\n directory, in FASTGARSIM_MACRO_PATH and next to this program."
              << "\n The command line wins over the job macro wherever both set the"
              << "\n same thing.\n" << std::endl;
}

// Everything the command line can say, before the job macro is read
struct CommandLine {
    std::string analysisMacro;
    std::string jobMacro;
    std::string inputFiles;
    std::string outputFile;
    std::string genieFiles;
    Long64_t maxEvents = -1;
    Bool_t maxEventsGiven = kFALSE;
};

// Anything that draws writes its canvases to file rather than opening a
// window, as there is no display attached to a batch job. Set
// FASTGARSIM_NO_BATCH for an analysis that really is meant to put something
// on the screen.
void SetBatch()
{
    if (std::getenv("FASTGARSIM_NO_BATCH") == nullptr) gROOT->SetBatch(kTRUE);
}

}  // namespace

int main(int argc, char** argv)
{
    std::cout << "\n==================================================" << std::endl;
    std::cout << "   GArAnalysis - Analysis Framework" << std::endl;
    std::cout << "==================================================" << std::endl;

    /* ---------------------------- Command line ---------------------------- */

    CommandLine given;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];

        const Bool_t hasValue = (i + 1 < argc);

        if (argument == "-h" || argument == "--help") {
            PrintUsage();
            return 0;
        } else if (argument == "-a" && hasValue) {
            given.analysisMacro = argv[++i];
        } else if (argument == "-i" && hasValue) {
            given.inputFiles = argv[++i];
        } else if (argument == "-o" && hasValue) {
            given.outputFile = argv[++i];
        } else if (argument == "-m" && hasValue) {
            given.jobMacro = argv[++i];
        } else if (argument == "-g" && hasValue) {
            given.genieFiles = argv[++i];
        } else if (argument == "-n" && hasValue) {
            given.maxEvents = std::atoll(argv[++i]);
            given.maxEventsGiven = kTRUE;
        } else if (!argument.empty() && argument[0] != '-' && given.analysisMacro.empty()) {
            // `GArAnalysis TruncatedDEDX.C -i ...` reads more naturally than
            // the flag does, and there is nothing else a bare argument in
            // front could mean
            given.analysisMacro = argument;
        } else {
            std::cerr << "\nUnknown or incomplete option: " << argument << std::endl;
            PrintUsage();
            return 1;
        }
    }

    /* ------------------------------ Job macro ----------------------------- */

    ana::AnalysisConfig config;
    std::string analysisMacro;
    std::string jobMacro;

    if (!given.jobMacro.empty()) {
        jobMacro = fastgarsim::FindMacro(given.jobMacro, argv[0]);

        if (!ana::ParseJobMacro(jobMacro, config, &analysisMacro)) {
            std::cerr << "\nError: the job macro '" << given.jobMacro
                      << "' could not be used; nothing was run.\n" << std::endl;
            return 1;
        }
    }

    // The command line is read after the macro, so that it wins
    if (!given.analysisMacro.empty()) analysisMacro = given.analysisMacro;
    if (!given.inputFiles.empty())    config.inputFiles = given.inputFiles;
    if (!given.outputFile.empty())    config.outputFile = given.outputFile;
    if (!given.genieFiles.empty())    config.genieFiles = given.genieFiles;
    if (given.maxEventsGiven)         config.maxEvents  = given.maxEvents;

    if (analysisMacro.empty()) {
        std::cerr << "\nError: no analysis to run. Name one with -a, or with"
                  << " /ana/global/analysis in the job macro." << std::endl;
        PrintUsage();
        return 1;
    }

    if (config.inputFiles.empty()) {
        std::cerr << "\nError: no input files. Name them with -i, or with"
                  << " /ana/global/input in the job macro." << std::endl;
        PrintUsage();
        return 1;
    }

    const std::string resolvedAnalysis = fastgarsim::FindMacro(analysisMacro, argv[0]);

    if (!std::filesystem::is_regular_file(resolvedAnalysis)) {
        std::cerr << "\nError: no analysis macro '" << analysisMacro << "'." << std::endl;
        std::cerr << "Looked in:" << std::endl;
        for (const std::string& directory : fastgarsim::MacroSearchPath(argv[0])) {
            std::cerr << "   " << directory << std::endl;
        }
        std::cerr << std::endl;
        return 1;
    }

    std::cout << "\n Configuration:" << std::endl;
    std::cout << "   Analysis:    " << resolvedAnalysis << std::endl;
    std::cout << "   Input:       " << config.inputFiles << std::endl;
    std::cout << "   Output:      "
              << (config.outputFile.empty() ? "(none)" : config.outputFile) << std::endl;
    if (!config.genieFiles.empty()) {
        std::cout << "   GENIE:       " << config.genieFiles << std::endl;
    }
    std::cout << "   Job macro:   "
              << (jobMacro.empty() ? "(none, defaults throughout)" : jobMacro) << std::endl;

    if (!config.params.IsEmpty()) {
        std::cout << "\n Analysis parameters:" << std::endl;
        config.params.Print(std::cout);
    }

    /* ------------------------- Compile and run it ------------------------- */

    SetBatch();

    // Where the analysis headers and libraries are, for ACLiC: where this
    // build put them, and where an installed or unpacked copy keeps them
    // relative to the program. Whatever the environment already says is
    // searched as well, so a build laid out differently still works once its
    // setup.sh has been sourced.
    const std::string executableDir = fastgarsim::ExecutableDir(argv[0]);

    const std::vector<std::string> includeDirs = {
        FASTGARSIM_ANALYSIS_INCLUDE_DIR,
        FASTGARSIM_COMMON_INCLUDE_DIR,
        executableDir + "/include",
        executableDir + "/../include"
    };

    const std::vector<std::string> libraryDirs = {
        FASTGARSIM_ANALYSIS_LIB_DIR,
        FASTGARSIM_COMMON_LIB_DIR,
        executableDir,
        executableDir + "/../lib"
    };

    ana::AnalysisBase* analysis =
        ana::CompileAnalysis(resolvedAnalysis, includeDirs, libraryDirs);
    if (!analysis) {
        ana::RemoveBuildDir();
        return 1;
    }

    const Bool_t succeeded = analysis->Execute(config);

    // The output is written and closed by Execute(), so by here the compiled
    // macro has done everything it was going to do
    delete analysis;
    ana::RemoveBuildDir();

    if (!succeeded) {
        std::cerr << "\n Analysis failed!\n" << std::endl;
        return 1;
    }

    std::cout << "\n==================================================" << std::endl;
    std::cout << "   GArAnalysis finished" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    return 0;
}
