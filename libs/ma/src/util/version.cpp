#include "version.h"

const std::string EXPORTED sLibMaVersion = MA_VERSION;

#ifdef WITH_PYTHON
const bool EXPORTED bLibMaWithPython = true;
#else
const bool EXPORTED bLibMaWithPython = false;
#endif