//
// AppOptions.hh - machinery behind the FastGArSim command line tools
//
// Every ROOT macro in the tool directories is also built as an executable. The
// executable is generated: CMake writes a one-line main() that #includes the
// macro and hands its entry function to app::Run() below. The macro itself is
// never modified and keeps working interpreted in ROOT.
//
// The awkward part is that a macro's entry function has default arguments, and
// default arguments are not part of a function's type. Taking a pointer to
// MakeNtuple would give a five-argument function with no way to recover the
// defaults. So the generated main() wraps the call in a
// generic lambda instead: inside it the call is written out normally, which is
// what makes the compiler apply the defaults. app::Run() then asks, at compile
// time, which argument counts that lambda actually accepts, and calls it with
// however many the user supplied.
//
#ifndef AppOptions_hh
#define AppOptions_hh

#include "TROOT.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

namespace app {

//---------------------------------------------------------------------------
// One command line argument, which converts itself to whatever type the macro
// declared for that parameter. This is what lets one generated main() drive
// macros with different signatures.
//
// Only the types the macros actually use are supported. A parameter of some
// other class type (a TString, say) would need two user-defined conversions
// and so fails to compile, which is the right answer: exclude that macro.
//---------------------------------------------------------------------------
class Arg {
public:
    explicit Arg(const char* text) : fText(text) {}

    operator const char*()  const { return fText; }
    operator std::string()  const { return fText; }
    operator double()       const { return std::atof(fText); }
    operator float()        const { return static_cast<float>(std::atof(fText)); }
    operator int()          const { return static_cast<int>(std::atoll(fText)); }
    operator long()         const { return static_cast<long>(std::atoll(fText)); }
    operator long long()    const { return std::atoll(fText); }
    operator unsigned int() const { return static_cast<unsigned int>(std::atoll(fText)); }
    operator bool()         const { return AsBool(); }

private:
    bool AsBool() const
    {
        const std::string value(fText);
        if (value == "1" || value == "true"  || value == "yes" || value == "kTRUE")  return true;
        if (value == "0" || value == "false" || value == "no"  || value == "kFALSE") return false;

        std::cerr << "Warning: cannot read '" << value
                  << "' as true or false, taking it as true" << std::endl;
        return true;
    }

    const char* fText;
};

//---------------------------------------------------------------------------
// Compile-time questions about the wrapped entry function
//---------------------------------------------------------------------------
namespace detail {

// N copies of Arg, for asking whether the macro takes N arguments
template <std::size_t> using ArgAt = Arg;

template <class F, std::size_t... I>
constexpr bool AcceptsImpl(std::index_sequence<I...>)
{
    return std::is_invocable_v<F&, ArgAt<I>...>;
}

template <class F, std::size_t N>
constexpr bool Accepts()
{
    return AcceptsImpl<F>(std::make_index_sequence<N>{});
}

template <class F, std::size_t... I>
bool CallWith(F& f, char** argv, std::index_sequence<I...>)
{
    if constexpr (std::is_invocable_v<F&, ArgAt<I>...>) {
        f(Arg(argv[I + 1])...);
        return true;
    } else {
        (void)f; (void)argv;
        return false;
    }
}

template <class F, std::size_t... N>
bool Dispatch(F& f, int n, char** argv, std::index_sequence<N...>)
{
    bool called = false;
    ((n == static_cast<int>(N) &&
      (called = CallWith(f, argv, std::make_index_sequence<N>{}))) || ...);
    return called;
}

// Lowest and highest argument counts the macro can be called with
template <class F, std::size_t... N>
void ArityRange(std::index_sequence<N...>, int& lowest, int& highest)
{
    const bool accepted[] = { Accepts<F, N>()... };

    lowest = -1;
    highest = -1;
    for (int i = 0; i < static_cast<int>(sizeof...(N)); ++i) {
        if (accepted[i]) {
            if (lowest < 0) lowest = i;
            highest = i;
        }
    }
}

}  // namespace detail

//---------------------------------------------------------------------------
// Helpers
//---------------------------------------------------------------------------
inline bool WantsHelp(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            return true;
        }
    }
    return false;
}

// No display is attached to a command line tool, so anything that draws writes
// its canvases to file rather than opening a window. Set FASTGARSIM_NO_BATCH to
// override, for a macro that really is meant to put something on the screen.
inline void Batch()
{
    if (std::getenv("FASTGARSIM_NO_BATCH") == nullptr) {
        gROOT->SetBatch(kTRUE);
    }
}

inline void PrintUsage(const char* name, const char* description, int lowest, int highest)
{
    std::cout << '\n' << description << "\n\nUsage: " << name;

    for (int i = 1; i <= highest; ++i) {
        std::cout << (i <= lowest ? " <arg" : " [arg") << i << (i <= lowest ? '>' : ']');
    }
    std::cout << "\n\n";

    if (highest > lowest) {
        std::cout << "Takes " << lowest << " to " << highest << " arguments; the optional "
                  << "ones default to\nwhatever " << name << ".C declares. "
                  << "See the package README for what each one means.\n" << std::endl;
    } else {
        std::cout << "Takes exactly " << lowest << " argument" << (lowest == 1 ? "" : "s")
                  << ". See the package README for what they mean.\n" << std::endl;
    }
}

//---------------------------------------------------------------------------
// Entry point used by the generated main()
//---------------------------------------------------------------------------
template <class F>
int Run(F f, int argc, char** argv, const char* name, const char* description)
{
    // One more than the largest number of parameters any macro declares
    constexpr std::size_t kSlots = 17;
    using Slots = std::make_index_sequence<kSlots>;

    int lowest = -1;
    int highest = -1;
    detail::ArityRange<F>(Slots{}, lowest, highest);

    if (WantsHelp(argc, argv)) {
        PrintUsage(name, description, lowest, highest);
        return 0;
    }

    const int given = argc - 1;
    if (given < lowest || given > highest) {
        std::cerr << '\n' << name << ": got " << given << " argument"
                  << (given == 1 ? "" : "s") << ", expected "
                  << lowest << " to " << highest << std::endl;
        PrintUsage(name, description, lowest, highest);
        return 1;
    }

    Batch();

    if (!detail::Dispatch(f, given, argv, Slots{})) {
        std::cerr << name << ": cannot be called with " << given << " arguments"
                  << std::endl;
        return 1;
    }
    return 0;
}

}  // namespace app

//---------------------------------------------------------------------------
// Written by the generated main(). FN is the macro's entry function; the
// generic lambda is what preserves its default arguments.
//---------------------------------------------------------------------------
#define FASTGARSIM_TOOL(FN, DESCRIPTION)                                       \
    int main(int argc, char** argv)                                            \
    {                                                                          \
        auto invoke = [](auto&&... a) -> decltype(FN(a...)) { return FN(a...); }; \
        return app::Run(invoke, argc, argv, #FN, DESCRIPTION);                 \
    }

#endif
