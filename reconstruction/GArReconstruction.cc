//
// GArReconstruction - Main executable for reconstruction algorithms
//

#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>

#include "RecoDataTypes.hh"
#include "RecoManager.hh"
#include "MacroPath.hh"

void PrintUsage() {
    std::cout << "\n Usage: GArReconstruction [options]" << std::endl;
    std::cout << "\n Options:" << std::endl;
    std::cout << "   -i <file>     Input ROOT file from simulation (required)" << std::endl;
    std::cout << "   -o <file>     Output ROOT file (default: the input with _reco before .root)" << std::endl;
    std::cout << "   -m <file>     Macro file for reconstruction configuration (required)" << std::endl;
    std::cout << "   -h, --help    Show this help message" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char** argv) {

    std::cout << "\n==================================================" << std::endl;
    std::cout << "   GArReconstruction - Reconstruction Framework" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Parse command line arguments
    std::string inputFile = "";
    std::string outputFile = "";
    std::string macroFile = "";

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            PrintUsage();
            return 0;
        }
        else if (arg == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        }
        else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        }
        else if (arg == "-m" && i + 1 < argc) {
            macroFile = argv[++i];
        }
        else {
            std::cerr << "Unknown option: " << arg << std::endl;
            PrintUsage();
            return 1;
        }
    }

    // Check if required arguments were provided
    if (inputFile.empty()) {
        std::cerr << "\nError: Input file required!" << std::endl;
        PrintUsage();
        return 1;
    }

    if (macroFile.empty()) {
        std::cerr << "\nError: Macro file required!" << std::endl;
        PrintUsage();
        return 1;
    }

    // The output is a copy of the input with the reconstruction added, so name
    // it after the input rather than after the program
    if (outputFile.empty()) {
        const size_t extension = inputFile.rfind(".root");
        outputFile = (extension == std::string::npos)
                   ? inputFile + "_reco.root"
                   : inputFile.substr(0, extension) + "_reco.root";
    }

    // The macros are copied next to the executable, so look for them there as
    // well as in the working directory; see common/include/MacroPath.hh
    const std::string resolvedMacro = fastgarsim::FindMacro(macroFile, argv[0]);

    std::cout << "\n Configuration:" << std::endl;
    std::cout << "   Input file:  " << inputFile << std::endl;
    std::cout << "   Output file: " << outputFile << std::endl;
    std::cout << "   Macro file:  " << resolvedMacro << std::endl;
    std::cout << std::endl;

    // Initialize reconstruction manager
    RecoManager* recoManager = new RecoManager();

    // Load reconstruction configuration from macro
    if (!recoManager->LoadMacro(resolvedMacro)) {
        std::cerr << "\nError: Failed to load macro file '" << macroFile << "'!" << std::endl;
        std::cerr << "Looked in:" << std::endl;
        for (const std::string& directory : fastgarsim::MacroSearchPath(argv[0])) {
            std::cerr << "   " << directory << std::endl;
        }
        delete recoManager;
        return 1;
    }

    // Run reconstruction
    std::cout << "\n Starting reconstruction..." << std::endl;
    bool success = recoManager->RunReconstruction(inputFile, outputFile);

    if (success) {
        std::cout << "\n Reconstruction completed successfully!" << std::endl;
        std::cout << "   Output saved to: " << outputFile << std::endl;
    } else {
        std::cerr << "\n Reconstruction failed!" << std::endl;
        delete recoManager;
        return 1;
    }

    // Cleanup
    delete recoManager;

    std::cout << "\n==================================================" << std::endl;
    std::cout << "   GArReconstruction finished" << std::endl;
    std::cout << "==================================================\n" << std::endl;

    return 0;
}
