//
// TPCClusterModule.hh - Grouping of TPC hits into clusters
//
// Runs DBSCAN over the reconstructed hit positions, with a uniform grid of
// cells the size of the neighbourhood radius so that neighbour lookups stay
// local. Hits that end up in no cluster are noise and are dropped.
//
// Each cluster is summarized by its charge-weighted centroid and spread, the
// principal axis of the hits it contains, and its extent along that axis --
// enough to seed a track fit or to feed a dE/dx estimate.
//
// Input:
//   "TPCHits"     <- std::vector<digi::TPCHit>
// Output:
//   "TPCClusters" -> std::vector<digi::TPCCluster>  (RecoStore + TTree)
//

#ifndef TPCClusterModule_h
#define TPCClusterModule_h 1

#include "RecoModule.hh"
#include "DigiDataTypes.hh"

#include <unordered_map>
#include <vector>

class TPCClusterModule : public RecoModule {
public:
    TPCClusterModule();
    virtual ~TPCClusterModule();

    virtual void Initialize() override;
    virtual void Execute()    override;
    virtual void Finalize()   override;

    virtual std::vector<ObjectSpec> GetInputSpec()  const override;
    virtual std::vector<ObjectSpec> GetOutputSpec() const override;

private:
    // Uniform grid over the hits, cells one neighbourhood radius across
    using CellMap = std::unordered_map<long long, std::vector<int>>;

    void BuildGrid(CellMap& grid) const;
    void Neighbours(int hitIndex, const CellMap& grid,
                    std::vector<int>& found) const;

    std::vector<std::vector<int>> RunDBSCAN() const;
    digi::TPCCluster MakeCluster(int clusterID,
                                 const std::vector<int>& members) const;

    // --- I/O ---
    const std::vector<digi::TPCHit>* fHits;
    std::vector<digi::TPCCluster>*   fClusters;

    // --- Parameters ---
    double fEpsilon;      // Neighbourhood radius [cm]
    int    fMinPoints;    // Neighbours needed for a hit to be a core hit
    int    fMinHits;      // Smallest cluster kept
    bool   fSeparatePlanes;  // Keep the two drift volumes' hits apart

    // --- Running totals for the end-of-job summary ---
    long fNEvents, fNClusters, fNClustered, fNHits;
};

#endif
