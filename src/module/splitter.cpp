/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "module/splitter.h"

using namespace libMA;


#ifdef WITH_PYTHON
void exportSplitter( )
{
    // export the Lock class
    exportModule<Lock<Container>>( "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( "UnLock" );
} // function
#endif
