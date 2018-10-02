/**
 * @file splitter.cpp
 * @author Markus Schmidt
 */
#include "module/splitter.h"

using namespace libMA;


#ifdef WITH_PYTHON
void exportSplitter( )
{
    // export the Splitter class
    boost::python::class_<Splitter, boost::python::bases<Module>, std::shared_ptr<Splitter>>(
        "Splitter", boost::python::init<std::shared_ptr<Pledge>>( )
        // make sure that the given pledge is not deallocated
        // before the splitter module
        //[boost::python::with_custodian_and_ward<1,2>()]
    );

    boost::python::implicitly_convertible<std::shared_ptr<Splitter>, std::shared_ptr<Module>>( );
    // export the Splitter class
    boost::python::class_<Collector, boost::python::bases<Module>, std::shared_ptr<Collector>>(
        "Collector", boost::python::init<std::shared_ptr<Container>>( ) )
        .def_readwrite( "content", &Collector::pVec );

    boost::python::implicitly_convertible<std::shared_ptr<Collector>, std::shared_ptr<Module>>( );

    // export the Lock class
    boost::python::class_<Lock, boost::python::bases<Module>, std::shared_ptr<Lock>>(
        "Lock", boost::python::init<std::shared_ptr<Container>>( ) );

    boost::python::implicitly_convertible<std::shared_ptr<Lock>, std::shared_ptr<Module>>( );

    // export the UnLock class
    boost::python::class_<UnLock, boost::python::bases<Module>, std::shared_ptr<UnLock>>(
        "UnLock", boost::python::init<std::shared_ptr<Pledge>>( ) );

    boost::python::implicitly_convertible<std::shared_ptr<UnLock>, std::shared_ptr<Module>>( );

    // export the ReadSplitter class
    boost::python::class_<ReadSplitter, boost::python::bases<Module>, std::shared_ptr<ReadSplitter>>( "ReadSplitter" );
    boost::python::implicitly_convertible<std::shared_ptr<ReadSplitter>, std::shared_ptr<Module>>( );

} // function
#endif
