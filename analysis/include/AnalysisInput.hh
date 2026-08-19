 /***************************************************************************
 * AnalysisInput.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   Turning an input specification into the concrete list of files an
 *   analysis should read: shell wildcards, comma-separated lists, and
 *   rewriting dCache /pnfs paths as XRootD URLs.
 *
 ***************************************************************************/

#ifndef AnalysisInput_hh
#define AnalysisInput_hh

#include <string>
#include <vector>

#include "Rtypes.h"

namespace ana {

/* -------------------------------------------------------------------------- */
/*                              XRootD settings                               */
/* -------------------------------------------------------------------------- */

// FNAL dCache defaults. A /pnfs/<rest> path is served as
//     <door><prefix><rest>
// i.e. /pnfs/dune/scratch/users/me/f.root becomes
//     root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/scratch/users/me/f.root
//
// Both can be overridden at run time with the FASTGAR_XROOTD_DOOR and
// FASTGAR_PNFS_PREFIX environment variables, for a different door or site.
extern const char* const kDefaultXRootDDoor;
extern const char* const kDefaultPnfsPrefix;

std::string XRootDDoor();
std::string PnfsPrefix();

// The door with its scheme stripped, e.g. "fndca1.fnal.gov:1094", as xrdfs
// and friends expect it
std::string XRootDHost();

/* -------------------------------------------------------------------------- */
/*                                Path helpers                                */
/* -------------------------------------------------------------------------- */

// Whether `path` is under the /pnfs dCache area
Bool_t IsPnfsPath(const std::string& path);

// Whether `path` is already a URL (root://, xroot://, http://, …) rather than
// a local file-system path
Bool_t IsRemoteURL(const std::string& path);

// Rewrite a /pnfs path as an XRootD URL. Anything else is returned unchanged,
// so this is safe to apply to a whole list.
std::string ToXRootD(const std::string& path);

/* -------------------------------------------------------------------------- */
/*                             Input specification                            */
/* -------------------------------------------------------------------------- */

// Expand an input specification into the list of files to read.
//
// The specification is a comma-separated list of entries, each of which may
// be a plain path or a shell wildcard pattern (*, ?, [...]). Patterns are
// expanded against the local file system, which covers NFS-mounted /pnfs
// areas; ~ and $VARIABLES are expanded first. Entries that are already URLs
// are passed through untouched, since a remote path cannot be listed this
// way. Duplicates are removed, so a file matched by two patterns is still
// read once.
//
// When `xrootdForPnfs` is true, every resolved /pnfs path is rewritten as an
// XRootD URL, so the files are streamed rather than read through the mount.
//
// Returns an empty vector, having said why, if nothing matched.
std::vector<std::string> ExpandInputFiles(const std::string& specification,
                                          Bool_t xrootdForPnfs = kTRUE);

// Print guidance for a file that ROOT could not open, tailored to whether it
// is a remote URL or a local path. Call this after ROOT has already reported
// its own error, so that an XRootD authentication failure is not mistaken for
// a missing file.
void ExplainOpenFailure(const std::string& path);

} // namespace ana

#endif
