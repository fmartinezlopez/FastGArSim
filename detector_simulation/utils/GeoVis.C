 /***************************************************************************
 * GeoVis.C
 * 
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 * 
 * Created: November 2025
 * 
 * Description:
 *   Draw geometry from GDML file with OpenGL
 *   Adapted from original by Brian Rebel
 * 
 ***************************************************************************/

#include <iostream>
#include <string>

#include "TGeoManager.h"
#include "TObjArray.h"

void GeoVis(const char* fileName, const char* volName="World_log", Bool_t checkOverlaps=true, Bool_t writeROOT=false) {
  
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  // Get file basename
  std::string basename(fileName);
  size_t lastdot = basename.find_last_of(".");
  if (lastdot != std::string::npos) {
      basename = basename.substr(0, lastdot);
  }
  
  // Import GDML file
  TGeoManager::Import(fileName);

  // Check for overlaps
  if (checkOverlaps) {
    gGeoManager->CheckOverlaps(0.01);
    gGeoManager->PrintOverlaps();
  }

  // Inspect volumes list
  TObjArray *vollist = gGeoManager->GetListOfVolumes();
  TIter next(vollist);
  TObject *obj = next();
  while (obj) {
    TString volName = obj->GetName();
    std::cout << volName << std::endl;
    obj = next();
  }

  // Draw geometry
  gGeoManager->FindVolumeFast(volName)->Draw("ogl");

  // Optionaly write geometry as ROOT file
  if (writeROOT) {
    std::string outFileName = basename + ".root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    gGeoManager->Write();
    outFile->Close();
  }

}
