/**
 * @file cyclic_queue_module.h
 * @brief Implements a modules that interact with a queue of containers.
 * @author Markus Schmidt
 */
#include "container/cyclic_queue_container.h"
#include "module/module.h"

namespace libMA
{

template <typename ContentType> class QueuePicker : public Module<ContentType, true, CyclicQueue<ContentType>>
{
  public:
    QueuePicker( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<ContentType> execute( std::shared_ptr<CyclicQueue<ContentType>> pQueue )
    {
        return pQueue->pop( );
    } // method
}; // class

// works with UnLock
template <typename ContentType>
class QueuePlacer : public Module<ContentType, false, ContentType, CyclicQueue<ContentType>>
{
  public:
    QueuePlacer( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<ContentType>
    execute( std::shared_ptr<ContentType> pContainer, std::shared_ptr<CyclicQueue<ContentType>> pQueue )
    {
        if( pContainer->eof( ) )
            pQueue->informThatContainerIsFinished( );
        else
            pQueue->push( pContainer );
        return pContainer;
    } // method
}; // class

template <typename ContentType>
std::shared_ptr<QueuePicker<ContentType>> getQueuePicker( const ParameterSetManager& rParameters,
                                                          std::shared_ptr<CyclicQueue<ContentType>> )
{
    return std::make_shared<QueuePicker<ContentType>>( rParameters );
} // function

template <typename ContentType>
std::shared_ptr<QueuePlacer<ContentType>> getQueuePlacer( const ParameterSetManager& rParameters,
                                                          std::shared_ptr<CyclicQueue<ContentType>> )
{
    return std::make_shared<QueuePlacer<ContentType>>( rParameters );
} // function


#ifdef WITH_PYTHON
template <typename ContentType, typename... ReturnTypes>
void exportCyclicQueue( py::module& rxPyModuleId, const std::string& sNamePrefix )
{
    py::class_<CyclicQueue<ContentType>, Container, std::shared_ptr<CyclicQueue<ContentType>>>(
        rxPyModuleId, std::string( sNamePrefix ).append( "Queue" ).c_str( ) )
        .def( py::init<std::shared_ptr<ContainerVector<std::shared_ptr<ContentType>>>>( ) )
        .def( py::init<>( ) )
        .def( "pop", &CyclicQueue<ContentType>::pop )
        .def( "push", &CyclicQueue<ContentType>::push )
        .def( "add", &CyclicQueue<ContentType>::add )
        .def( "inform_finished", &CyclicQueue<ContentType>::informThatContainerIsFinished );

    exportModule<QueuePicker<ContentType>>( rxPyModuleId, std::string( sNamePrefix ).append( "Picker" ).c_str( ) );
    exportModule<QueuePlacer<ContentType>>( rxPyModuleId, std::string( sNamePrefix ).append( "Placer" ).c_str( ) );
} // function

#endif

} // namespace libMA