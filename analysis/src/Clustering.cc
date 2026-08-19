 /***************************************************************************
 * Clustering.cc
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Implementation of the grid-accelerated DBSCAN clustering.
 *
 ***************************************************************************/

#include "Clustering.hh"

#include <cmath>
#include <functional>
#include <unordered_map>

namespace ana {

namespace {

struct GridKey {
    Int_t x, y, z;
    Bool_t operator==(const GridKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct KeyHash {
    size_t operator()(const GridKey& k) const {
        return ((51 + std::hash<Int_t>()(k.x)) * 51 + std::hash<Int_t>()(k.y)) * 51
               + std::hash<Int_t>()(k.z);
    }
};

using Grid = std::unordered_map<GridKey, std::vector<Int_t>, KeyHash>;

GridKey CellOf(const TVector3& p, Double_t eps)
{
    return GridKey{ Int_t(std::floor(p.X() / eps)),
                    Int_t(std::floor(p.Y() / eps)),
                    Int_t(std::floor(p.Z() / eps)) };
}

// Assign points to a 3D grid of cells of side eps, for fast neighbour lookup
Grid BuildGrid(const std::vector<TVector3>& points, Double_t eps)
{
    Grid grid;
    for (size_t i = 0; i < points.size(); ++i) {
        grid[CellOf(points[i], eps)].push_back(static_cast<Int_t>(i));
    }
    return grid;
}

// Retrieve the neighbours of point `idx` using the grid
std::vector<Int_t> RegionQuery(const std::vector<TVector3>& points,
                               const Grid& grid,
                               Int_t idx,
                               Double_t eps)
{
    std::vector<Int_t> neighbors;
    const TVector3& p = points[idx];
    const GridKey cell = CellOf(p, eps);

    for (Int_t dx = -1; dx <= 1; ++dx) {
        for (Int_t dy = -1; dy <= 1; ++dy) {
            for (Int_t dz = -1; dz <= 1; ++dz) {
                auto it = grid.find(GridKey{cell.x + dx, cell.y + dy, cell.z + dz});
                if (it == grid.end()) continue;
                for (Int_t j : it->second) {
                    if ((points[j] - p).Mag() <= eps) neighbors.push_back(j);
                }
            }
        }
    }
    return neighbors;
}

} // anonymous namespace

/* -------------------------------------------------------------------------- */
/*                               Cluster helpers                              */
/* -------------------------------------------------------------------------- */

TVector3 Cluster::Centroid(const std::vector<TVector3>& points) const
{
    TVector3 sum(0., 0., 0.);
    if (indices.empty()) return sum;
    for (Int_t i : indices) sum += points[i];
    return sum * (1.0 / indices.size());
}

TVector3 Cluster::Centroid(const std::vector<TVector3>& points,
                           const std::vector<Float_t>& weights) const
{
    TVector3 sum(0., 0., 0.);
    Double_t total = 0.;
    for (Int_t i : indices) {
        const Double_t w = weights[i];
        sum += w * points[i];
        total += w;
    }
    if (total <= 0.) return Centroid(points);
    return sum * (1.0 / total);
}

Double_t Cluster::TotalWeight(const std::vector<Float_t>& weights) const
{
    Double_t total = 0.;
    for (Int_t i : indices) total += weights[i];
    return total;
}

Double_t Cluster::Radius(const std::vector<TVector3>& points) const
{
    const TVector3 centre = Centroid(points);
    Double_t maxDistance = 0.;
    for (Int_t i : indices) {
        const Double_t d = (points[i] - centre).Mag();
        if (d > maxDistance) maxDistance = d;
    }
    return maxDistance;
}

/* -------------------------------------------------------------------------- */
/*                                   DBSCAN                                   */
/* -------------------------------------------------------------------------- */

std::vector<Cluster> DBSCAN3D(const std::vector<TVector3>& points,
                              Double_t eps,
                              Int_t minPts)
{
    std::vector<Cluster> clusters;
    if (points.empty() || eps <= 0.) return clusters;

    const Int_t n = static_cast<Int_t>(points.size());
    std::vector<Int_t> labels(n, -1);  // -1 = unvisited, -2 = noise
    Int_t clusterID = 0;

    const Grid grid = BuildGrid(points, eps);

    for (Int_t i = 0; i < n; ++i) {
        if (labels[i] != -1) continue;

        auto neighbors = RegionQuery(points, grid, i, eps);
        if (static_cast<Int_t>(neighbors.size()) < minPts) {
            labels[i] = -2;  // noise
            continue;
        }

        Cluster cluster;
        cluster.indices.push_back(i);
        labels[i] = clusterID;

        std::vector<Int_t> seeds = neighbors;
        for (size_t j = 0; j < seeds.size(); ++j) {
            const Int_t idx = seeds[j];

            if (labels[idx] == -2) {
                // Previously flagged as noise, now a border point of this cluster
                labels[idx] = clusterID;
                cluster.indices.push_back(idx);
                continue;
            }
            if (labels[idx] != -1) continue;

            labels[idx] = clusterID;
            cluster.indices.push_back(idx);

            auto subNeighbors = RegionQuery(points, grid, idx, eps);
            if (static_cast<Int_t>(subNeighbors.size()) >= minPts) {
                seeds.insert(seeds.end(), subNeighbors.begin(), subNeighbors.end());
            }
        }

        clusters.push_back(std::move(cluster));
        clusterID++;
    }

    return clusters;
}

} // namespace ana
