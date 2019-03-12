/**
 * @file hashMapSeeding.h
 * @brief Computes the mapping quality of Alignments.
 * @author Markus Schmidt
 */
#pragma once

#include "container/alignment.h"
#include "module/module.h"

namespace libMA
{
/**
 * @brief 
 * @ingroup module
 * @details
 */
class HashMapSeeding : public Module<Seeds, false, NucSeq, NucSeq>
{
  public:
    const size_t uiSeedSize;

    HashMapSeeding( const ParameterSetManager& rParameters )
        : uiSeedSize( rParameters.getSelected( )->xSecSeedSize->get( ) )
    {} // constructor

    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 );

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the HashMapSeeding @ref Module "module" to python.
 * @ingroup export
 */
void exportHashMapSeeding( py::module& rxPyModuleId );
#endif
