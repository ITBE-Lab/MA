/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "module/splitter.h"

using namespace libMA;


#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportSplitter( )
{
    // export the Lock class
    exportModule<Lock<Container>>( "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( "UnLock" );


    // export the TupleGet class
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 0>>( "GetFirstQuery" );
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 1>>( "GetSecondQuery" );
} // function
#else
void exportSplitter( py::module& rxPyModuleId )
{
    // export the Lock class
    exportModule<Lock<Container>>( rxPyModuleId, "Lock" );

    // export the UnLock class
    exportModule<UnLock<Container>, std::shared_ptr<BasePledge>>( rxPyModuleId, "UnLock" );


    // export the TupleGet class
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 0>>( rxPyModuleId, "GetFirstQuery" );
    exportModule<TupleGet<ContainerVector<std::shared_ptr<NucSeq>>, 1>>( rxPyModuleId, "GetSecondQuery" );
} // function
#endif
#endif
