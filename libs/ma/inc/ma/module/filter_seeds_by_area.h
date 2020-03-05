/**
 * @file filter_seeds_by_area.h
 * @brief implements FilterSeedsByArea
 * @author Markus Schmidt
 */
#include "ms/module/module.h"
#include "ma/container/segment.h"
namespace libMS
{

/**
 * @brief ?
 */
class FilterSeedsByArea : public Module<Seeds, false, SegmentVector, FMIndex, NucSeq>
{
  public:
    nucSeqIndex uiStart;
    nucSeqIndex uiSize;
    nucSeqIndex uiMaxAmbiguity;
    nucSeqIndex uiSeedSize;

    FilterSeedsByArea( const ParameterSetManager& rParameters, nucSeqIndex uiStart, nucSeqIndex uiSize )
        : uiStart( uiStart ),
          uiSize( uiSize ),
          uiMaxAmbiguity( rParameters.getSelected( )->xMaxAmbiguitySv->get( ) ),
          uiSeedSize( rParameters.getSelected( )->xMinSeedSizeSV->get( ) )
    {} // constructor

    typename std::shared_ptr<Seeds> execute( std::shared_ptr<SegmentVector> pSegments,
                                             std::shared_ptr<FMIndex> pFmIndex, std::shared_ptr<NucSeq> pQuery )
    {
        std::shared_ptr<Seeds> pRet = std::make_shared<Seeds>( );
        pSegments->forEachSeed( *pFmIndex, pQuery->length( ), uiMaxAmbiguity, uiSeedSize, true,
                                [&]( Seed& s ) {
                                    s.uiDelta = (nucSeqIndex)pQuery->iId;
                                    if( s.start_ref( ) <= uiStart + uiSize && s.end_ref( ) >= uiStart )
                                        pRet->push_back( s );
                                    return true;
                                } // lambda
        ); // for each
        return pRet;
    } // method
}; // class
} // namespace libMS