//
// TPCClusterModule.cc - Grouping of TPC hits into clusters
//

#include "TPCClusterModule.hh"
#include "ModuleFactory.hh"
#include "RecoStore.hh"

#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>

// Self-registration
namespace {
    bool kRegistered = ModuleFactory::Instance().Register(
        "TPCClusterModule",
        []() -> RecoModule* { return new TPCClusterModule(); });

    // Iterations of the power method used to find the principal axis
    constexpr int kPowerIterations = 32;
}

// ---------------------------------------------------------------------------
// Constructor / Destructor
// ---------------------------------------------------------------------------

TPCClusterModule::TPCClusterModule()
    : RecoModule("TPCCluster"),
      fHits(nullptr), fClusters(nullptr),
      fEpsilon(2.0), fMinPoints(3), fMinHits(3), fSeparatePlanes(true),
      fNEvents(0), fNClusters(0), fNClustered(0), fNHits(0)
{
}

TPCClusterModule::~TPCClusterModule()
{
}

// ---------------------------------------------------------------------------
// Module spec
// ---------------------------------------------------------------------------

std::vector<ObjectSpec> TPCClusterModule::GetInputSpec() const
{
    return {{"TPCHits", "std::vector<digi::TPCHit>"}};
}

std::vector<ObjectSpec> TPCClusterModule::GetOutputSpec() const
{
    return {{"TPCClusters", "std::vector<digi::TPCCluster>"}};
}

// ---------------------------------------------------------------------------
// Initialize
// ---------------------------------------------------------------------------

void TPCClusterModule::Initialize()
{
    Print("Initializing TPC clustering module");

    fEpsilon        = GetParameterDouble("epsilon",        2.0);
    fMinPoints      = GetParameterInt   ("minPoints",      3);
    fMinHits        = GetParameterInt   ("minHits",        3);
    fSeparatePlanes = GetParameterBool  ("separatePlanes", true);

    fHits = fStore->Get<std::vector<digi::TPCHit>>("TPCHits");

    Print("Clustering parameters:");
    std::cout << "   Neighbourhood radius [cm]: " << fEpsilon   << "\n"
              << "   Minimum neighbours: "        << fMinPoints << "\n"
              << "   Minimum cluster size: "      << fMinHits   << "\n"
              << "   Planes clustered separately: "
              << (fSeparatePlanes ? "yes" : "no") << std::endl;

    fClusters = new std::vector<digi::TPCCluster>();
    fStore->Register("TPCClusters", fClusters);

    if (fOutputTree) {
        fOutputTree->Branch("TPCClusters", &fClusters);
    }
}

// ---------------------------------------------------------------------------
// Execute
// ---------------------------------------------------------------------------

void TPCClusterModule::Execute()
{
    fClusters->clear();
    if (!fHits) return;

    int clusterID = 0;
    for (const auto& members : RunDBSCAN()) {
        if (static_cast<int>(members.size()) < fMinHits) continue;
        fClusters->push_back(MakeCluster(clusterID++, members));
        fNClustered += members.size();
    }

    fNEvents++;
    fNHits     += fHits->size();
    fNClusters += fClusters->size();
}

// ---------------------------------------------------------------------------
// Finalize
// ---------------------------------------------------------------------------

void TPCClusterModule::Finalize()
{
    Print("TPC clustering summary:");
    std::cout << "   Events processed: " << fNEvents   << "\n"
              << "   Hits seen: "        << fNHits     << "\n"
              << "   Clusters found: "   << fNClusters << std::endl;
    if (fNHits > 0) {
        std::cout << "   Fraction of hits clustered: "
                  << double(fNClustered) / fNHits << std::endl;
    }
    if (fNEvents > 0) {
        std::cout << "   Clusters per event: "
                  << double(fNClusters) / fNEvents << std::endl;
    }

    delete fClusters;
    fClusters = nullptr;
}

// ---------------------------------------------------------------------------
// Neighbour lookup
// ---------------------------------------------------------------------------

namespace {
    // Pack three cell indices into one key. The ranges involved (a few hundred
    // cells per axis for a 250 cm TPC) fit comfortably in 21 bits each.
    inline long long CellKey(int i, int j, int k)
    {
        return ((long long)(i + (1 << 20)) << 42)
             | ((long long)(j + (1 << 20)) << 21)
             |  (long long)(k + (1 << 20));
    }
}

void TPCClusterModule::BuildGrid(CellMap& grid) const
{
    for (size_t i = 0; i < fHits->size(); ++i) {
        const auto& hit = (*fHits)[i];
        grid[CellKey(static_cast<int>(std::floor(hit.x / fEpsilon)),
                     static_cast<int>(std::floor(hit.y / fEpsilon)),
                     static_cast<int>(std::floor(hit.z / fEpsilon)))]
            .push_back(static_cast<int>(i));
    }
}

void TPCClusterModule::Neighbours(int hitIndex, const CellMap& grid,
                                  std::vector<int>& found) const
{
    found.clear();

    const auto& hit = (*fHits)[hitIndex];
    const int ci = static_cast<int>(std::floor(hit.x / fEpsilon));
    const int cj = static_cast<int>(std::floor(hit.y / fEpsilon));
    const int ck = static_cast<int>(std::floor(hit.z / fEpsilon));
    const double eps2 = fEpsilon * fEpsilon;

    // Cells are one radius across, so the neighbourhood is inside the 27 cells
    // touching this one
    for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
            for (int dk = -1; dk <= 1; ++dk) {
                auto cell = grid.find(CellKey(ci + di, cj + dj, ck + dk));
                if (cell == grid.end()) continue;

                for (int other : cell->second) {
                    const auto& candidate = (*fHits)[other];
                    if (fSeparatePlanes && candidate.plane != hit.plane) continue;

                    const double dx = candidate.x - hit.x;
                    const double dy = candidate.y - hit.y;
                    const double dz = candidate.z - hit.z;
                    if (dx*dx + dy*dy + dz*dz <= eps2) found.push_back(other);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// DBSCAN
// ---------------------------------------------------------------------------

std::vector<std::vector<int>> TPCClusterModule::RunDBSCAN() const
{
    std::vector<std::vector<int>> clusters;
    const int n = static_cast<int>(fHits->size());
    if (n == 0) return clusters;

    CellMap grid;
    BuildGrid(grid);

    enum State { kUnvisited, kNoise, kAssigned };
    std::vector<State> state(n, kUnvisited);

    std::vector<int> neighbours, more, members;

    for (int i = 0; i < n; ++i) {
        if (state[i] != kUnvisited) continue;

        Neighbours(i, grid, neighbours);
        if (static_cast<int>(neighbours.size()) < fMinPoints) {
            state[i] = kNoise;   // May still be picked up as a border hit
            continue;
        }

        // Grow the cluster out from this core hit
        members.clear();
        members.push_back(i);
        state[i] = kAssigned;

        std::vector<int> queue = neighbours;
        for (size_t q = 0; q < queue.size(); ++q) {
            const int hit = queue[q];
            if (state[hit] == kAssigned) continue;

            const bool wasNoise = (state[hit] == kNoise);
            state[hit] = kAssigned;
            members.push_back(hit);

            // Border hits join the cluster but do not extend it
            if (wasNoise) continue;

            Neighbours(hit, grid, more);
            if (static_cast<int>(more.size()) >= fMinPoints) {
                for (int candidate : more) {
                    if (state[candidate] != kAssigned) queue.push_back(candidate);
                }
            }
        }

        clusters.push_back(members);
    }

    return clusters;
}

// ---------------------------------------------------------------------------
// Cluster properties
// ---------------------------------------------------------------------------

digi::TPCCluster TPCClusterModule::MakeCluster(int clusterID,
                                               const std::vector<int>& members) const
{
    digi::TPCCluster cluster;
    cluster.clusterID = clusterID;
    cluster.nHits     = static_cast<int>(members.size());
    cluster.hitIndices.assign(members.begin(), members.end());

    // --- Charge-weighted centroid ---
    double weightSum = 0;
    double cx = 0, cy = 0, cz = 0;
    for (int index : members) {
        const auto& hit = (*fHits)[index];
        const double w  = std::max(0.0f, hit.charge);

        weightSum += w;
        cx += w * hit.x;
        cy += w * hit.y;
        cz += w * hit.z;

        cluster.charge     += hit.charge;
        cluster.energy     += hit.energy;
        cluster.trueEnergy += hit.trueEnergy;
    }

    // Fall back on an unweighted centroid if nothing carried charge
    if (weightSum <= 0) {
        weightSum = members.size();
        cx = cy = cz = 0;
        for (int index : members) {
            cx += (*fHits)[index].x;
            cy += (*fHits)[index].y;
            cz += (*fHits)[index].z;
        }
    }
    cluster.x = cx / weightSum;
    cluster.y = cy / weightSum;
    cluster.z = cz / weightSum;

    // --- Spread and covariance about the centroid ---
    double covariance[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double sumWeights = 0;
    for (int index : members) {
        const auto& hit = (*fHits)[index];
        const double w  = (hit.charge > 0) ? hit.charge : 1.0;
        const double d[3] = {hit.x - cluster.x, hit.y - cluster.y, hit.z - cluster.z};

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) covariance[a][b] += w * d[a] * d[b];
        }
        sumWeights += w;
    }
    if (sumWeights > 0) {
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) covariance[a][b] /= sumWeights;
        }
    }
    cluster.rmsX = std::sqrt(std::max(0.0, covariance[0][0]));
    cluster.rmsY = std::sqrt(std::max(0.0, covariance[1][1]));
    cluster.rmsZ = std::sqrt(std::max(0.0, covariance[2][2]));

    // --- Principal axis, by power iteration on the covariance matrix ---
    double axis[3] = {cluster.rmsX + 1e-9, cluster.rmsY, cluster.rmsZ};
    for (int iteration = 0; iteration < kPowerIterations; ++iteration) {
        double next[3] = {0, 0, 0};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) next[a] += covariance[a][b] * axis[b];
        }
        const double norm = std::sqrt(next[0]*next[0] + next[1]*next[1] + next[2]*next[2]);
        if (norm <= 0) break;
        for (int a = 0; a < 3; ++a) axis[a] = next[a] / norm;
    }
    cluster.dirX = axis[0];
    cluster.dirY = axis[1];
    cluster.dirZ = axis[2];

    // --- Extent along the axis, and the hits at either end of it ---
    double minProjection = 0, maxProjection = 0;
    int    minIndex = members.front(), maxIndex = members.front();
    double transverse = 0;
    bool   first = true;

    for (int index : members) {
        const auto& hit = (*fHits)[index];
        const double d[3] = {hit.x - cluster.x, hit.y - cluster.y, hit.z - cluster.z};
        const double along = d[0]*axis[0] + d[1]*axis[1] + d[2]*axis[2];

        const double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        transverse += std::max(0.0, d2 - along * along);

        if (first || along < minProjection) { minProjection = along; minIndex = index; }
        if (first || along > maxProjection) { maxProjection = along; maxIndex = index; }
        first = false;
    }

    cluster.length = maxProjection - minProjection;
    cluster.width  = std::sqrt(transverse / members.size());

    cluster.startX = (*fHits)[minIndex].x;
    cluster.startY = (*fHits)[minIndex].y;
    cluster.startZ = (*fHits)[minIndex].z;
    cluster.endX   = (*fHits)[maxIndex].x;
    cluster.endY   = (*fHits)[maxIndex].y;
    cluster.endZ   = (*fHits)[maxIndex].z;

    // --- Which readout plane, and which tracks, the hits came from ---
    cluster.plane = (*fHits)[members.front()].plane;
    std::vector<std::pair<int, float>> tracks;
    double trackTotal = 0;

    for (int index : members) {
        const auto& hit = (*fHits)[index];
        if (hit.plane != cluster.plane) cluster.plane = -1;

        for (size_t t = 0; t < hit.trackIDs.size(); ++t) {
            const float energy = hit.trueEnergy * hit.trackFractions[t];
            trackTotal += energy;

            auto found = std::find_if(tracks.begin(), tracks.end(),
                [&](const std::pair<int, float>& entry) {
                    return entry.first == hit.trackIDs[t]; });
            if (found != tracks.end()) found->second += energy;
            else tracks.emplace_back(hit.trackIDs[t], energy);
        }
    }

    std::sort(tracks.begin(), tracks.end(),
              [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
                  return a.second > b.second; });
    for (const auto& track : tracks) {
        cluster.trackIDs.push_back(track.first);
        cluster.trackFractions.push_back(trackTotal > 0 ? track.second / trackTotal : 0);
    }

    return cluster;
}
