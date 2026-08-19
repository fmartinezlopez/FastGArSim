//
// RecoDataTypesLinkDef.h - ROOT dictionary linkdef for reconstruction data types
//

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ struct RecoTrack+;
#pragma link C++ struct RecoCluster+;
#pragma link C++ struct RecoVertex+;
#pragma link C++ struct RecoEvent+;

#pragma link C++ class std::vector<RecoTrack>+;
#pragma link C++ class std::vector<RecoCluster>+;
#pragma link C++ class std::vector<RecoVertex>+;

#endif
