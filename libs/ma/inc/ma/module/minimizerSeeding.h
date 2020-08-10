/**
 * @file minimizerSeeding.h
 * @brief links libMA to the minimiap 2 code...
 * @author Markus Schmidt
 */
#ifndef MINIMIZER_SEEDING_H
#define MINIMIZER_SEEDING_H

#include "ma/container/minimizer_index.h"
#include "ma/container/pack.h"
#include "ma/container/seed.h"
#include "ms/module/module.h"
#include "util/system.h"

#ifdef WITH_ZLIB
namespace libMA
{

/**
 * @brief Computes a maximally covering set of seeds.
 * @details
 * Can use either the extension scheme by Li et Al. or ours.
 * @ingroup module
 */
class MinimizerSeeding : public libMS::Module<Seeds, false, minimizer::Index, NucSeq, Pack>
{
  public:
    /**
     * @brief Initialize a MinimizerSeeding Module
     * @details
     */
    MinimizerSeeding( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<Seeds> execute( std::shared_ptr<minimizer::Index> pMMIndex, std::shared_ptr<NucSeq> pQuery,
                                    std::shared_ptr<Pack> pPack )
    {
        pQuery->vTranslateToCharacterForm( );
        const char* sSeq = (const char*)pQuery->pxSequenceRef;
        const int iSize = (int)pQuery->length( );
        auto pRet = pMMIndex->seed_one( sSeq, iSize, pPack );
        pQuery->vTranslateToNumericForm( );
        return pRet;
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief exports the Segmentation @ref Module "module" to python.
 * @ingroup export
 */
void exportMinimizerSeeding( libMS::SubmoduleOrganizer& rxPyModuleId );
#endif
#endif

#endif