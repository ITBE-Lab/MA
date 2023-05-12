/**
 * @file coverageUpdater.h
 * @brief 
 * @author Markus Schmidt
 */
#pragma once

#include "ma/container/seed.h"
#include "ms/module/get_inserter_container_module.h"
#include "msv/container/sv_db/tables/coverage.h"
#include "util/geom.h"

using namespace libMA;
using namespace libMS;

namespace libMSV
{
/**
 * @brief Insertion of coverage into the database.
 */
template <typename DBCon>
class CoverageUpdaterContainer
    : public InserterContainer<DBCon, AbstractInserterContainer, CoverageTable, Seeds>
{
  public:
    using ParentType =
        InserterContainer<DBCon, libMS::AbstractInserterContainer, CoverageTable, Seeds>;
    const size_t uiBinSize;


    CoverageUpdaterContainer( std::shared_ptr<PoolContainer<DBCon>> pPool, int64_t iId,
                           std::shared_ptr<SharedInserterProfiler> pSharedProfiler )
        : ParentType::InserterContainer( pPool, iId, pSharedProfiler ),
        uiBinSize(pGlobalParams->xCoverageBinSize->get())
    {} // constructur

  protected:
    virtual size_t insert_override( std::shared_ptr<Seeds> pSeeds )
    {
        std::set<int64_t> xCoveredBinIds;
        for( Seed& rSeed : *pSeeds )
        {
            const nucSeqIndex uiStartBinId = rSeed.start_ref_cons_rev() / uiBinSize;
            const nucSeqIndex uiEndBinId = 1 + (rSeed.end_ref_cons_rev() - 1) / uiBinSize;
            
            for(nucSeqIndex uiX = uiStartBinId; uiX < uiEndBinId; uiX++)
                xCoveredBinIds.insert(uiX);
        } // for
    
        for(int64_t uiBinId : xCoveredBinIds)
            ParentType::pInserter->incCoverage( ParentType::iId /* <- id of the sv jump run */, uiBinId, (uint32_t)1 );
        return xCoveredBinIds.size( );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetCoverageUpdaterContainerModule =
    GetInserterContainerModule<CoverageUpdaterContainer, DBCon, DBConInit, CoverageTable>;


template <typename DBCon> using CoverageUpdaterModule = InserterModule<CoverageUpdaterContainer<DBCon>>;


} // namespace libMSV

#ifdef WITH_PYTHON
/// @brief used to expose libMSV::CoverageUpdaterContainer to python
void exportCoverageUpdater( libMS::SubmoduleOrganizer& xOrganizer );
#endif
