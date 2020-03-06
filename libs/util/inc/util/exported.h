
#ifdef __GNUC__
// under gnu DLL_PORT is not needed to do anything
#define DLL_PORT( libName )
#elif _MSC_VER


// under msc we need to export or import a function according to weather we are building the dll or
// using it
#define DLL_PORT( libName ) __declspec( DLL_PORT_##libName )

#endif // _MSC_VER