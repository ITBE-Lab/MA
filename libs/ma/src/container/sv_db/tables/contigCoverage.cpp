
#include "container/sv_db/tables/contigCoverage.h"
#include "container/sv_db/svDb.h"

using namespace libMA;

ContigCovTable::CovInserter::CovInserter( int64_t iSequencerId, std::shared_ptr<Pack> pPack,
                                          std::shared_ptr<SV_DB> pDb )
    : iSequencerId( iSequencerId ), pDb( pDb ), pPack( pPack ), vNumNts( pPack->uiNumContigs( ) )
{
    if( iSequencerId != -1 )
        pDb->pContigCovTable->resetCount( iSequencerId );
} // constructor

void ContigCovTable::CovInserter::commit( )
{
    if( iSequencerId == -1 )
        return;
    std::lock_guard<std::mutex> xGuard( *pDb->pWriteLock );
    for( size_t uiI = 0; uiI < vNumNts.size( ); uiI++ )
        if( vNumNts[ uiI ] > 0 )
        {
            pDb->pContigCovTable->incrementNt( iSequencerId, uiI, vNumNts[ uiI ] );
            vNumNts[ uiI ] = 0;
        } // if
} // function