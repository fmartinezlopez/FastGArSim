//
// RecoDataTypes.hh - Data structures for reconstruction output
//

#ifndef RecoDataTypes_h
#define RecoDataTypes_h 1

#include <vector>

// Reconstructed track structure
struct RecoTrack {
    int trackID;
    double startX, startY, startZ;
    double endX, endY, endZ;
    double momentum;
    double theta, phi;
    double length;
    int nHits;
    double chi2;
    int pdgCode;  // Particle hypothesis

    RecoTrack() : trackID(0), startX(0), startY(0), startZ(0),
                  endX(0), endY(0), endZ(0), momentum(0),
                  theta(0), phi(0), length(0), nHits(0),
                  chi2(0), pdgCode(0) {}
};

// Reconstructed cluster structure
struct RecoCluster {
    int clusterID;
    double energy;
    double centerX, centerY, centerZ;
    double timeAvg;
    int nHits;
    double rmsX, rmsY, rmsZ;

    RecoCluster() : clusterID(0), energy(0), centerX(0), centerY(0), centerZ(0),
                    timeAvg(0), nHits(0), rmsX(0), rmsY(0), rmsZ(0) {}
};

// Reconstructed vertex structure
struct RecoVertex {
    int vertexID;
    double x, y, z;
    double time;
    int nTracks;
    std::vector<int> trackIDs;
    double chi2;

    RecoVertex() : vertexID(0), x(0), y(0), z(0), time(0),
                   nTracks(0), chi2(0) {}
};

// Event-level reconstruction information
struct RecoEvent {
    int eventID;
    double totalEnergy;
    int nTracks;
    int nClusters;
    int nVertices;

    RecoEvent() : eventID(0), totalEnergy(0), nTracks(0),
                  nClusters(0), nVertices(0) {}
};

#endif
