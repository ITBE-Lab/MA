
#ifdef __GNUC__
// under gnu EXPORTED is not needed to do anything
#define EXPORTED
#elif _MSC_VER
// under msc we need to export or import a function according to weather we are building the dll or
// using it
#ifdef EXPORT
#define EXPORTED __declspec( dllexport )
#else
#define EXPORTED __declspec( dllimport )
#endif
#endif