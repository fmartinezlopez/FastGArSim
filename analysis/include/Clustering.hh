 /***************************************************************************
 * Clustering.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Density-based (DBSCAN) clustering of 3D points, with a uniform grid for
 *   fast neighbour lookup. Used to group calorimeter hits into clusters.
 *
 ***************************************************************************/

#ifndef Clustering_hh
#define Clustering_hh

#include <vector>

#include "Rtypes.h"
#include "TVector3.h"

namespace ana {

/* -------------------------------------------------------------------------- */
/*                            Clustering algorithm                            */
/* -------------------------------------------------------------------------- */

struct Cluster {
    // Indices into the point collection that was clustered
    std::vector<Int_t> indices;

    size_t Size() const { return indices.size(); }

    // Unweighted centroid of the cluster
    TVector3 Centroid(const std::vector<TVector3>& points) const;

    // Centroid weighted by `weights` (e.g. hit energy deposits). Falls back to
    // the unweighted centroid if the total weight is zero.
    TVector3 Centroid(const std::vector<TVector3>& points,
                      const std::vector<Float_t>& weights) const;

    // Sum of `weights` over the cluster members
    Double_t TotalWeight(const std::vector<Float_t>& weights) const;

    // Largest distance between any member and the centroid
    Double_t Radius(const std::vector<TVector3>& points) const;
};

// DBSCAN clustering of `points`.
//   eps    -- neighbourhood radius, in the same units as the points (cm)
//   minPts -- minimum number of neighbours for a point to be a core point
// Points that end up in no cluster are treated as noise and simply omitted
// from the result.
std::vector<Cluster> DBSCAN3D(const std::vector<TVector3>& points,
                              Double_t eps,
                              Int_t minPts);

} // namespace ana

#endif
