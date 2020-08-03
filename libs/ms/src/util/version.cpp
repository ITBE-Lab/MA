#include "ms/util/version.h"

const std::string DLL_PORT(MS) sLibMaVersion = MA_VERSION;

#ifdef WITH_PYTHON
const bool DLL_PORT(MS) bLibMaWithPython = true;
#else
const bool DLL_PORT(MS) bLibMaWithPython = false;
#endif
