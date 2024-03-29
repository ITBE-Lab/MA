/**
 * @file alignment.cpp
 * @author Markus Schmidt
 */
#include "ma/container/alignment.h"
#include "ms/util/parameter.h"
#include "ms/util/pybind11.h"
using namespace libMA;


void DLL_PORT( MA ) Alignment::append( MatchType type, nucSeqIndex size )
{
#if DEBUG_LEVEL >= 2
    // get a copy of the alignment for later comparison in case something goes wrong
    std::vector<std::pair<MatchType, nucSeqIndex>> vCopyOfData( data.begin( ), data.end( ) );
    // const char vTranslate[5] = {'S', '=', 'X', 'I', 'D'};
    // std::cout << vTranslate[type] << size << std::endl;
#endif
    /*
     * we are storing in a compressed format
     */
    if( size == 0 )
        return;
    // adjust the score of the alignment
    if( type == MatchType::seed || type == MatchType::match )
    {
        iScore += pGlobalParams->iMatch->get( ) * size;
        uiEndOnRef += size;
        uiEndOnQuery += size;
    } // if
    else if( type == MatchType::missmatch )
    {
        // pGlobalParams->iMissMatch->get() is a penalty not a score
        iScore -= pGlobalParams->iMissMatch->get( ) * size;
        uiEndOnRef += size;
        uiEndOnQuery += size;
    } // else if
    else if( type == MatchType::insertion || type == MatchType::deletion )
    {
        // add the sizes first so we do not have to revert them
        if( type == MatchType::insertion )
            uiEndOnQuery += size;
        else
            uiEndOnRef += size;
        // pGlobalParams->iGap->get() & pGlobalParams->iExtend->get() is a penalty not a score
        if( data.size( ) != 0 && ( data.back( ).first == type ) )
        {
            // revert last alignment but memorize the length
            size += data.back( ).second;
            uiLength -= data.back( ).second;
            if( pGlobalParams->iExtend->get( ) * data.back( ).second + pGlobalParams->iGap->get( ) <
                (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
                iScore += pGlobalParams->iExtend->get( ) * data.back( ).second + pGlobalParams->iGap->get( );
            else
                iScore += (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
            data.pop_back( );
        } // if
        // add the penalty for this indel (plus the possibly removed last one...)
        if( pGlobalParams->iExtend->get( ) * size + pGlobalParams->iGap->get( ) <
            (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
            iScore -= pGlobalParams->iExtend->get( ) * size + pGlobalParams->iGap->get( );
        else
            iScore -= (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
    } // else if

    /*
     * we are storing in a compressed format
     * since it actually makes quite a lot of things easier
     * same thing here: just check weather the last symbol is the same as the inserted one
     *      if so just add the amount of new symbols
     *      else make a new entry with the correct amount of symbols
     */
    if( data.size( ) != 0 && data.back( ).first == type )
        data.back( ).second += size;
    else
        data.push_back( std::make_pair( type, size ) );
    uiLength += size;

    DEBUG_2( if( reCalcScore( ) != iScore ) {
        std::cerr << "WARNING set wrong score in append name: " << xStats.sName << " actual score: " << reCalcScore( )
                  << " score: " << iScore << std::endl;
        for( auto tup : vCopyOfData )
            std::cout << std::get<0>( tup ) << ":" << std::get<1>( tup ) << " ";
        std::cout << std::endl;
        for( auto tup : data )
            std::cout << std::get<0>( tup ) << ":" << std::get<1>( tup ) << " ";
        std::cout << std::endl;
        assert( false );
    } // if
             ) // DEBUG
    DEBUG( nucSeqIndex uiCheck = 0; for( auto xTup
                                         : data ) uiCheck += xTup.second;
           if( uiCheck != uiLength ) {
               std::cout << "Alignment length check failed: " << uiCheck << " != " << uiLength << std::endl;
               assert( false );
           } // if
           ) // DEBUG
} // function

unsigned int DLL_PORT( MA ) Alignment::localscore( ) const
{
    unsigned int uiMaxScore = 0;
    int64_t iScoreCurr = 0;
    for( unsigned int index = 0; index < data.size( ); index++ )
    {
        switch( data[ index ].first )
        {
            case MatchType::deletion:
            case MatchType::insertion:
                if( pGlobalParams->iExtend->get( ) * data[ index ].second + pGlobalParams->iGap->get( ) <
                    (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
                    iScoreCurr -= pGlobalParams->iExtend->get( ) * data[ index ].second + pGlobalParams->iGap->get( );
                else
                    iScoreCurr -= (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
                break;
            case MatchType::missmatch:
                iScoreCurr -= pGlobalParams->iMissMatch->get( ) * data[ index ].second;
                break;
            case MatchType::match:
            case MatchType::seed:
                iScoreCurr += pGlobalParams->iMatch->get( ) * data[ index ].second;
                break;
        } // switch
        if( iScoreCurr < 0 )
            iScoreCurr = 0;
        if( uiMaxScore < (unsigned int)iScoreCurr )
            uiMaxScore = (unsigned int)iScoreCurr;
    } // for
    return uiMaxScore;
} // function

void DLL_PORT( MA ) Alignment::makeLocal( )
{
    if( uiLength == 0 )
        return;
    std::vector<int> vScores;
    int64_t iMaxScore = 0;
    nucSeqIndex iMaxStart = 0;
    nucSeqIndex iMaxEnd = data.size( );
    nucSeqIndex iLastStart = 0;
    int64_t iScoreCurr = 0;
    /*
     * the purpose of this loop is to set iMaxStart & iMaxEnd correctly.
     * this is done using an approach similar to SW backtracking:
     * - run from the beginning to the end of the alignment
     * - always keep track of the current score
     * - always keep track of the last index where the score was zero
     * - if the current score is the largest so far encountered score do the following:
     *      - overwrite iMaxStart with the last index where the score was zero
     *      - overwrite iMaxEnd with the current index
     * once the loop finished iMaxStart & iMaxEnd are set correctly.
     */
    for( unsigned int index = 0; index < data.size( ); index++ )
    {
        switch( data[ index ].first )
        {
            case MatchType::deletion:
            case MatchType::insertion:
                if( pGlobalParams->iExtend->get( ) * data[ index ].second + pGlobalParams->iGap->get( ) <
                    (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
                    iScoreCurr -= pGlobalParams->iExtend->get( ) * data[ index ].second + pGlobalParams->iGap->get( );
                else
                    iScoreCurr -= (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
                break;
            case MatchType::missmatch:
                iScoreCurr -= pGlobalParams->iMissMatch->get( ) * data[ index ].second;
                break;
            case MatchType::match:
            case MatchType::seed:
                iScoreCurr += pGlobalParams->iMatch->get( ) * data[ index ].second;
                break;
        } // switch
        if( iScoreCurr < 0 )
        {
            iScoreCurr = 0;
            iLastStart = index + 1;
        } // if
        DEBUG_2( std::cout << data[ index ].first << "," << data[ index ].second << " (" << iScoreCurr << ") | "; )
        if( iScoreCurr >= iMaxScore )
        {
            iMaxScore = iScoreCurr;
            iMaxStart = iLastStart;
            iMaxEnd = index + 1;
        } // if
    } // for
    DEBUG_2( std::cout << std::endl; std::cout << iMaxStart << " " << iMaxEnd << std::endl; )
    // adjust the begin/end on ref/query according to the area that will be erased
    if( iMaxStart <= iMaxEnd )
    {
        for( nucSeqIndex index = 0; index < iMaxStart; index++ )
            switch( data[ index ].first )
            {
                case MatchType::deletion:
                    uiBeginOnRef += data[ index ].second;
                    break;
                case MatchType::insertion:
                    uiBeginOnQuery += data[ index ].second;
                    break;
                default:
                    uiBeginOnRef += data[ index ].second;
                    uiBeginOnQuery += data[ index ].second;
                    break;
            } // switch
        for( nucSeqIndex index = iMaxEnd; index < data.size( ); index++ )
            switch( data[ index ].first )
            {
                case MatchType::deletion:
                    uiEndOnRef -= data[ index ].second;
                    break;
                case MatchType::insertion:
                    uiEndOnQuery -= data[ index ].second;
                    break;
                default:
                    uiEndOnRef -= data[ index ].second;
                    uiEndOnQuery -= data[ index ].second;
                    break;
            } // switch
    } // if
    // erase everything before and after
    if( iMaxEnd < data.size( ) )
        data.erase( data.begin( ) + iMaxEnd, data.end( ) );
    if( iMaxStart > 0 )
        data.erase( data.begin( ), data.begin( ) + iMaxStart );
    // adjust score accordingly
    iScore = iMaxScore;
    // adjust length variable accordingly
    uiLength = 0;
    for( unsigned int index = 0; index < data.size( ); index++ )
        uiLength += data[ index ].second;
    DEBUG_2( if( uiEndOnRef < uiBeginOnRef ) {
        std::cout << "---" << std::endl;
        for( unsigned int index = 0; index < data.size( ); index++ )
            std::cout << data[ index ].first << "," << data[ index ].second << std::endl;
        exit( 0 );
    } )
    DEBUG( if( reCalcScore( ) != iScore ) std::cerr << "WARNING set wrong score or removed wrong elements in makeLocal"
                                                    << std::endl; )
} // function

void DLL_PORT( MA ) Alignment::removeDangeling( )
{
#if DEBUG_LEVEL >= 1
    // get a copy of the alignment for later comparison in case something goes wrong
    std::vector<std::tuple<MatchType, nucSeqIndex>> vCopyOfData( data.begin( ), data.end( ) );
#endif

    if( data.empty( ) )
        return;
    while( data.front( ).first == MatchType::deletion || data.front( ).first == MatchType::insertion )
    {
        if( data.front( ).first == MatchType::deletion ) // deletion
            uiBeginOnRef += data.front( ).second;
        else // insertion
            uiBeginOnQuery += data.front( ).second;
        if( pGlobalParams->iGap->get( ) + pGlobalParams->iExtend->get( ) * data.front( ).second <
            (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
            iScore += pGlobalParams->iGap->get( ) + pGlobalParams->iExtend->get( ) * data.front( ).second;
        else
            iScore += (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
        uiLength -= data.front( ).second;
        data.erase( data.begin( ), data.begin( ) + 1 );
    } // if
    while( data.back( ).first == MatchType::deletion || data.back( ).first == MatchType::insertion )
    {
        if( data.back( ).first == MatchType::deletion ) // deletion
            uiEndOnRef -= data.back( ).second;
        else // insertion
            uiEndOnQuery -= data.back( ).second;
        if( pGlobalParams->iGap->get( ) + pGlobalParams->iExtend->get( ) * data.back( ).second <
            (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
            iScore += pGlobalParams->iGap->get( ) + pGlobalParams->iExtend->get( ) * data.back( ).second;
        else
            iScore += (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
        uiLength -= data.back( ).second;
        data.pop_back( );
    } // if
    //iScore = reCalcScore();
    DEBUG( if( reCalcScore( ) != iScore ) {
        std::cerr << "WARNING set wrong score or removed wrong elements in remove dangeling" << std::endl;
        for( auto tup : vCopyOfData )
            std::cout << std::get<0>( tup ) << ":" << std::get<1>( tup ) << " ";
        std::cout << std::endl;
        for( auto tup : data )
            std::cout << std::get<0>( tup ) << ":" << std::get<1>( tup ) << " ";
        std::cout << std::endl;
    } // if

           nucSeqIndex uiCheck = 0;
           for( auto xTup
                : data ) uiCheck += xTup.second;
           if( uiCheck != uiLength ) {
               std::cout << "Alignment length check failed: " << uiCheck << " != " << uiLength << std::endl;
               assert( false );
           } // if
           ) // DEBUG
} // function

int64_t Alignment::reCalcScore( ) const
{
    int64_t iScore = 0;
    //int64_t iSVScore = -pGlobalParams->uiSVPenalty->get( );
    for( unsigned int index = 0; index < data.size( ); index++ )
    {
        //if(iScore < iSVScore)
        //    iScore = iSVScore;
        switch( std::get<0>( data[ index ] ) )
        {
            case MatchType::deletion:
            case MatchType::insertion:
                if( pGlobalParams->iExtend->get( ) * data[ index ].second + pGlobalParams->iGap->get( ) <
                    (nucSeqIndex)pGlobalParams->uiSVPenalty->get( ) )
                    iScore -= pGlobalParams->iExtend->get( ) * data[ index ].second + pGlobalParams->iGap->get( );
                else
                    iScore -= (nucSeqIndex)pGlobalParams->uiSVPenalty->get( );
                break;
            case MatchType::missmatch:
                iScore -= pGlobalParams->iMissMatch->get( ) * data[ index ].second;
                break;
            case MatchType::match:
            case MatchType::seed:
                iScore += pGlobalParams->iMatch->get( ) * data[ index ].second;
                break;
        } // switch
        //if(iScore - pGlobalParams->uiSVPenalty->get( ) > iSVScore)
        //    iSVScore = iScore - pGlobalParams->uiSVPenalty->get( );
    }
    return iScore;
} // function


#ifdef WITH_PYTHON
#include <pybind11/stl.h>

void exportAlignment( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<Alignment, libMS::Container, std::shared_ptr<Alignment>>( xOrganizer.container(), "Alignment" )
        .def( "at", &Alignment::at )
        .def( py::init<nucSeqIndex>( ) )
        .def( py::init<nucSeqIndex, nucSeqIndex>( ) )
        .def( py::init<nucSeqIndex, nucSeqIndex, nucSeqIndex, nucSeqIndex>( ) )
        .def( "__getitem__", &Alignment::at )
        .def( "append", &Alignment::append_boost1 )
        .def( "append", &Alignment::append_boost2 )
        .def( "begin_on_ref", &Alignment::beginOnRef )
        .def( "end_on_ref", &Alignment::endOnRef )
        .def( "__len__", &Alignment::length )
        .def( "length", &Alignment::length )
        .def( "seed_coverage", &Alignment::seedCoverage )
        .def( "get_score", &Alignment::score )
        .def( "get_local_score", &Alignment::localscore )
        .def( "num_seeds", &Alignment::getNumSeeds )
        .def( "num_by_seeds", &Alignment::getNumBySeeds )
        .def( "make_local", &Alignment::makeLocal )
        .def( "extract", &Alignment::extract )
        .def( "overlap", &Alignment::overlap )
        .def( "cigarString", &Alignment::cigarString )
        .def( "getSamFlag", &Alignment::getSamFlag )
        .def( "getContig", &Alignment::getContig )
        .def( "getSamPosition", &Alignment::getSamPosition )
        .def( "getQuerySequence", &Alignment::getQuerySequence )
        .def( "to_seeds", &Alignment::toSeeds )
        .def( "set_other", &Alignment::setOther )
        .def_readonly( "stats", &Alignment::xStats )
        .def_readonly( "data", &Alignment::data )
        .def_readwrite( "begin_on_query", &Alignment::uiBeginOnQuery )
        .def_readwrite( "end_on_query", &Alignment::uiEndOnQuery )
        .def_readwrite( "begin_on_ref", &Alignment::uiBeginOnRef )
        .def_readwrite( "end_on_ref", &Alignment::uiEndOnRef )
        .def_readwrite( "mapping_quality", &Alignment::fMappingQuality )
        .def_readwrite( "secondary", &Alignment::bSecondary )
        .def_readwrite( "supplementary", &Alignment::bSupplementary )
            DEBUG(.def_readwrite( "vGapsScatter", &Alignment::vGapsScatter ) );


    // export the matchType enum
    py::enum_<MatchType>( xOrganizer.util(), "MatchType" )
        .value( "match", MatchType::match )
        .value( "seed", MatchType::seed )
        .value( "missmatch", MatchType::missmatch )
        .value( "insertion", MatchType::insertion )
        .value( "deletion", MatchType::deletion )
        .export_values( );

    py::bind_vector<std::vector<MatchType>>( xOrganizer._util(), "MatchTypeVector", "docstr" );

    py::bind_vector_ext<libMS::ContainerVector<std::shared_ptr<Alignment>>,
                        libMS::Container,
                        std::shared_ptr<libMS::ContainerVector<std::shared_ptr<Alignment>>>>(
        xOrganizer.container(), "AlignmentVector", "docstr" );

    xOrganizer.container().def("toAlignmentVector", [](std::shared_ptr<libMS::Container> pPtr){
        return std::dynamic_pointer_cast<libMS::ContainerVector<std::shared_ptr<Alignment>>>(pPtr);
    });

    // tell boost python that pointers of these classes can be converted implicitly
    py::implicitly_convertible<libMS::Container, libMS::ContainerVector<std::shared_ptr<Alignment>>>( );

} // function
#endif