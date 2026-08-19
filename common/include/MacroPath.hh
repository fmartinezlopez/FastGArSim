//
// MacroPath.hh - locating the run macros
//
// The macros are copied next to their executable at build time and refer to
// each other by paths relative to it (/control/execute macros/init.mac). If
// they were only ever resolved against the working directory, the programs
// would have to be started from their own build directory. Both GArSimulation
// and GArReconstruction therefore search, in order:
//
//   1. the working directory, so a local macros/ still wins and nothing that
//      already works changes behaviour;
//   2. FASTGARSIM_MACRO_PATH, exported by the generated setup.sh;
//   3. the directory holding the executable and its macros/ subdirectory,
//      which is what makes a relocated, installed or tarballed build work --
//      the grid jobs unpack exactly such a copy.
//
#ifndef MacroPath_hh
#define MacroPath_hh

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

namespace fastgarsim {

//---------------------------------------------------------------------------
// Directory holding the running executable, resolved from argv[0]. When the
// program was found through PATH there is no directory component, so PATH is
// searched the same way the shell did.
//---------------------------------------------------------------------------
inline std::string ExecutableDir(const char* argv0)
{
    namespace fs = std::filesystem;

    fs::path exe(argv0 != nullptr ? argv0 : "");

    if (!exe.has_parent_path()) {
        if (const char* path = std::getenv("PATH")) {
            std::stringstream candidates(path);
            std::string dir;
            while (std::getline(candidates, dir, ':')) {
                if (dir.empty()) continue;
                std::error_code ec;
                const fs::path candidate = fs::path(dir) / exe;
                if (fs::is_regular_file(candidate, ec)) {
                    exe = candidate;
                    break;
                }
            }
        }
    }

    std::error_code ec;
    fs::path resolved = fs::weakly_canonical(exe, ec);
    if (ec) resolved = exe;

    return resolved.has_parent_path() ? resolved.parent_path().string() : ".";
}

//---------------------------------------------------------------------------
// The directories to look in, most specific first
//---------------------------------------------------------------------------
inline std::vector<std::string> MacroSearchPath(const char* argv0)
{
    std::vector<std::string> directories;
    directories.push_back(".");

    if (const char* fromEnv = std::getenv("FASTGARSIM_MACRO_PATH")) {
        std::stringstream entries(fromEnv);
        std::string entry;
        while (std::getline(entries, entry, ':')) {
            if (!entry.empty()) directories.push_back(entry);
        }
    }

    const std::string exeDir = ExecutableDir(argv0);
    directories.push_back(exeDir);
    directories.push_back(exeDir + "/macros");

    return directories;
}

// Colon separated, for handing to Geant4's macro search path
inline std::string MacroSearchPathString(const char* argv0)
{
    std::string joined;
    for (const std::string& directory : MacroSearchPath(argv0)) {
        if (!joined.empty()) joined += ":";
        joined += directory;
    }
    return joined;
}

//---------------------------------------------------------------------------
// First readable match for name, or name unchanged when nothing matches, so
// that the caller reports the name the user actually typed.
//---------------------------------------------------------------------------
inline std::string FindMacro(const std::string& name, const char* argv0)
{
    namespace fs = std::filesystem;

    if (name.empty()) return name;

    std::error_code ec;
    const fs::path asGiven(name);

    // An absolute path, or one that already resolves, is taken as it stands
    if (asGiven.is_absolute() || fs::is_regular_file(asGiven, ec)) {
        return name;
    }

    for (const std::string& directory : MacroSearchPath(argv0)) {
        const fs::path candidate = fs::path(directory) / asGiven;
        if (fs::is_regular_file(candidate, ec)) {
            return candidate.string();
        }
    }

    return name;
}

}  // namespace fastgarsim

#endif
