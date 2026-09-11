 /***************************************************************************
 * DumpSchema.C
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Print what a FastGArSim ROOT file contains: every tree and branch with
 *   its C++ type, which module produced it, and the processing passes the
 *   file has been through. Use it to find out what an analysis can ask for
 *   before writing the analysis.
 *
 * Usage:
 *   DumpSchema <file> [showMacros]
 *
 *   showMacros also prints the configuration macro of each pass, verbatim.
 *
 ***************************************************************************/

R__LOAD_LIBRARY(libGArCommon)

#include <iostream>

#include "TFile.h"
#include "TSystem.h"

#include "ProductSchema.hh"

void DumpSchema(const char* fileName, Bool_t showMacros = kFALSE)
{
#ifdef __CLING__
    gSystem->Load("libSimDataDict");
    gSystem->Load("libDigiDataDict");
#endif

    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open " << fileName << std::endl;
        return;
    }

    const fastgarsim::ProductSchema schema = fastgarsim::ProductSchema::Load(file);

    std::cout << "\n" << fileName << std::endl;
    schema.Print(std::cout);

    if (showMacros) {
        for (const fastgarsim::JobInfo& job : schema.Jobs()) {
            if (job.macroText.empty()) continue;
            std::cout << "\n--- " << job.stage << " macro (" << job.macro << ") ---\n"
                      << job.macroText << std::endl;
        }
    }

    std::cout << std::endl;
    file->Close();
    delete file;
}
