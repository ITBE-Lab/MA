/**
 * @file cyclic_queue_module.h
 * @brief Implements a modules that interact with a queue of containers.
 * @author Markus Schmidt
 */
#include "container/cyclic_queue_container.h"
#include "module/module.h"

namespace libMS
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
template <typename ReturnType, typename ContentType>
class QueuePlacer : public Module<ReturnType, false, ReturnType, ContentType, CyclicQueue<ContentType>>
{
  public:
    QueuePlacer( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<ReturnType> execute( std::shared_ptr<ReturnType> pRet,
                                         std::shared_ptr<ContentType>
                                             pContainer,
                                         std::shared_ptr<CyclicQueue<ContentType>>
                                             pQueue )
    {
        if( pContainer->eof( ) )
        {
            pContainer->close( );
            pQueue->informThatContainerIsFinished( );
        } // if
        else
            pQueue->push( pContainer );
        return pRet;
    } // method
}; // class

template <typename ContentType>
std::shared_ptr<QueuePicker<ContentType>> getQueuePicker( const ParameterSetManager& rParameters,
                                                          std::shared_ptr<CyclicQueue<ContentType>> )
{
    return std::make_shared<QueuePicker<ContentType>>( rParameters );
} // function

template <typename ReturnType, typename ContentType>
std::shared_ptr<QueuePlacer<ReturnType, ContentType>> getQueuePlacer( const ParameterSetManager& rParameters,
                                                                      std::shared_ptr<CyclicQueue<ContentType>> )
{
    return std::make_shared<QueuePlacer<ReturnType, ContentType>>( rParameters );
} // function


#ifdef WITH_PYTHON

template <size_t uiN, typename ContentType>
void exportCyclicQueueHelper( py::module& rxPyModuleId, const std::string& sNamePrefix,
                              const std::vector<std::string>& vReturnTypeNames )
{} // function

template <size_t uiN, typename ContentType, typename ReturnType, typename... ReturnTypes>
void exportCyclicQueueHelper( py::module& rxPyModuleId, const std::string& sNamePrefix,
                              const std::vector<std::string>& vReturnTypeNames )
{
    exportModule<QueuePlacer<ReturnType, ContentType>>(
        rxPyModuleId,
        std::string( sNamePrefix ).append( uiN == 0 ? "" : vReturnTypeNames[ uiN - 1 ] ).append( "Placer" ).c_str( ) );
    exportCyclicQueueHelper<uiN + 1, ContentType, ReturnTypes...>( rxPyModuleId, sNamePrefix, vReturnTypeNames );
} // function

/** @brief expose a cyclic queue to python
 * @details
 * This exposes:
 * - CyclicQueue<ContentType> under the name sNamePrefix + "Queue"
 * - QueuePicker<ContentType> under the name sNamePrefix + "Picker"
 * - QueuePlacer<ContentType, ContentType> under the name sNamePrefix + "Placer"
 * - QueuePlacer<ReturnTypes[i], ContentType> under the name sNamePrefix + vReturnTypeNames[i] + "Placer"
 *   where i is in the range of 0 to size(ReturnTypes...)
 * vReturnTypeNames must be the same size as ReturnTypes..., both can be 0 in that case vReturnTypeNames default
 * value can be used.
 * This is purely a convenience function.
 * Example usage:
    exportCyclicQueue<FileStream, NucSeq>( rxPyModuleId, "File", {"NucSeq"} );
    exportCyclicQueue<PairedFileStream, PairedReadsContainer>( rxPyModuleId, "PairedFile", {"NucSeq"} );
 * In both cases one extra placer is exposed: the FileNucSeqPlacer and the PairedFileNucSeqPlacer
 */
template <typename ContentType, typename... ReturnTypes>
void exportCyclicQueue( py::module& rxPyModuleId, const std::string& sNamePrefix,
                        const std::vector<std::string>& vReturnTypeNames = {} )
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
    exportCyclicQueueHelper<0, ContentType, ContentType, ReturnTypes...>( rxPyModuleId, sNamePrefix, vReturnTypeNames );
} // function

#endif

} // namespace libMS