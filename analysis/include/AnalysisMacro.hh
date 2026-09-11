 /***************************************************************************
 * AnalysisMacro.hh
 *
 * Author: Francisco Martinez Lopez
 * Email: frmart@iu.edu
 *
 * Description:
 *   The two things GArAnalysis does before an analysis can run: read the job
 *   macro that configures it, and compile the analysis macro that defines it.
 *
 *   JOB MACRO
 *   A job macro is the analysis counterpart of the reconstruction's .mac
 *   files, and is read the same way: one command per line, '#' starts a
 *   comment, blank lines are ignored.
 *
 *     # What to run over
 *     /ana/global/input       'gun_*.root'
 *     /ana/global/output      dedx.root
 *     /ana/global/outputTree  DEDX
 *     /ana/global/maxEvents   1000
 *
 *     # What this particular analysis makes of it
 *     /ana/truncation         0.6
 *     /ana/minHits            10
 *
 *   Everything under /ana/global/ configures the job itself and is written
 *   into AnalysisConfig, so those names are fixed and checked here. Anything
 *   else under /ana/ belongs to the analysis, which reads it in Configure();
 *   those names are whatever that analysis chose to read.
 *
 *   ANALYSIS MACRO
 *   The analysis itself is an ordinary ROOT macro ending in ANA_ANALYSIS(),
 *   compiled with ACLiC when the job starts rather than when the framework is
 *   built. It is compiled afresh every time, in a private build directory
 *   that is thrown away afterwards, so a job always runs the code that is in
 *   the file at that moment and nothing is left behind next to the macro.
 *
 ***************************************************************************/

#ifndef AnalysisMacro_hh
#define AnalysisMacro_hh

#include <string>
#include <vector>

#include "Rtypes.h"

#include "AnalysisBase.hh"

namespace ana {

/* -------------------------------------------------------------------------- */
/*                                 Job macro                                  */
/* -------------------------------------------------------------------------- */

// Read `path` into `config`: the /ana/global/ commands fill its fields and
// everything else under /ana/ goes into config.params for the analysis to
// read in Configure(). When `analysisMacro` is given it receives
// /ana/global/analysis, so that one file can describe a whole job.
//
// Returns kFALSE, having said what is wrong, for a macro that cannot be read,
// a command that is not understood, a /ana/global/ name that does not exist
// or a value that cannot be read as the type that name requires. Whatever was
// parsed before the error is still in `config`, which is of no use to a
// caller that stops -- as GArAnalysis does.
Bool_t ParseJobMacro(const std::string& path,
                     AnalysisConfig& config,
                     std::string* analysisMacro = nullptr);

// The /ana/global/ names ParseJobMacro understands, in the order they are
// documented. For error messages and for --help.
const std::vector<std::string>& JobMacroGlobals();

/* -------------------------------------------------------------------------- */
/*                              Analysis macro                                */
/* -------------------------------------------------------------------------- */

// Compile `path` with ACLiC and return a new instance of the analysis class
// its ANA_ANALYSIS() line names. The caller owns the object.
//
// `includeDirs` and `libraryDirs` are added to the compiler's header search
// path and to ROOT's library search path, so that the analysis headers are
// found and the data type dictionaries autoload, without the environment
// having been set up first.
//
// Returns nullptr, having said what is wrong, when the macro does not exist,
// does not compile, or compiles but declares no analysis.
AnalysisBase* CompileAnalysis(const std::string& path,
                              const std::vector<std::string>& includeDirs = {},
                              const std::vector<std::string>& libraryDirs = {});

// Remove the private build directory CompileAnalysis() compiled into. It is
// removed at exit in any case; call this once the job's output is written and
// closed to have it gone at a point of your choosing. Doing nothing twice is
// safe.
void RemoveBuildDir();

}  // namespace ana

#endif
