/**
 * @file seedFilters.h
 * @brief Implements several filters for seeds
 * @author Markus Schmidt
 */
#pragma once

namespace libMA
{

#if 1
/**
 * @brief Extends seeds at the end to create maximally extened seeds
 * @ingroup module
 */
class SeedExtender : public libMS::Module<Seeds, false, Seeds, NucSeq, Pack>
{
    template <typename Func1_t, typename Func2_t>
    static void extendSeedHelper( Seed& rSeed, Func1_t&& fGetNucQuery, Func2_t&& fGetNucRef, size_t uiQSize,
                                  size_t uiRSize )
    {
        size_t uiForw = 1;
        if( rSeed.bOnForwStrand )
            while( uiForw <= rSeed.start( ) && //
                   uiForw <= rSeed.start_ref( ) && //
                   fGetNucQuery( rSeed.start( ) - uiForw ) == fGetNucRef( rSeed.start_ref( ) - uiForw ) )
                uiForw++;
        else
            while( uiForw <= rSeed.start( ) && //
                   uiForw + rSeed.start_ref( ) < uiRSize && //
                   fGetNucQuery( rSeed.start( ) - uiForw ) == 3 - fGetNucRef( rSeed.start_ref( ) + uiForw ) )
                uiForw++;
        uiForw--; // uiForward is one too high after loop
        rSeed.iSize += uiForw;
        rSeed.iStart -= uiForw;
        if( rSeed.bOnForwStrand )
            rSeed.uiPosOnReference -= uiForw;
        else
            rSeed.uiPosOnReference += uiForw;

        size_t uiBackw = 0;
        if( rSeed.bOnForwStrand )
            while( uiBackw + rSeed.end( ) < uiQSize && //
                   uiBackw + rSeed.end_ref( ) < uiRSize && //
                   fGetNucQuery( rSeed.end( ) + uiBackw ) == fGetNucRef( rSeed.end_ref( ) + uiBackw ) )
                uiBackw++;
        else
            while( uiBackw + rSeed.end( ) < uiQSize && //
                   rSeed.start_ref( ) >= uiBackw + rSeed.size( ) && //
                   fGetNucQuery( rSeed.end( ) + uiBackw ) ==
                       3 - fGetNucRef( rSeed.start_ref( ) - rSeed.size( ) - uiBackw ) )
                uiBackw++;
        rSeed.iSize += uiBackw;
    }

  public:
    SeedExtender( const ParameterSetManager& rParameters )
    {} // default constructor

    static void extendSeed( Seed& rSeed, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRef )
    {
        extendSeedHelper(
            rSeed, [&]( size_t uiI ) { return pQuery->pxSequenceRef[ uiI ]; },
            [&]( size_t uiI ) { return pRef->vExtract( uiI ); }, pQuery->length( ),
            pRef->uiUnpackedSizeForwardPlusReverse( ) );
    }

    static void extendSeed( Seed& rSeed, std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
    {
        extendSeedHelper(
            rSeed, [&]( size_t uiI ) { return pQ1->pxSequenceRef[ uiI ]; },
            [&]( size_t uiI ) { return pQ2->pxSequenceRef[ uiI ]; }, pQ1->length( ), pQ2->length( ) );
    }

    static void extendSeed( Seed& rSeed, NucSeq& rQ1, NucSeq& rQ2 )
    {
        extendSeedHelper(
            rSeed, [&]( size_t uiI ) { return rQ1.pxSequenceRef[ uiI ]; },
            [&]( size_t uiI ) { return rQ2.pxSequenceRef[ uiI ]; }, rQ1.length( ), rQ2.length( ) );
    }

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRef )
    {
        for( auto& rSeed : *pSeeds )
            extendSeed( rSeed, pQuery, pRef );
        return pSeeds;
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> extend( std::vector<std::shared_ptr<libMA::Seeds>> vIn,
                                                               std::vector<std::shared_ptr<NucSeq>>
                                                                   vQueries,
                                                               std::shared_ptr<Pack>
                                                                   pRef )
    {
        if( vIn.size( ) != vQueries.size( ) )
            throw std::runtime_error( "vIn and vQueries have different lenghts" );
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ], vQueries[ uiI ], pRef ) );
        return vRet;
    } // method
}; // class
#endif

/**
 * @brief Combines overlapping seeds
 * @ingroup module
 * @details
 * Uses a n log(n) algorithm to combine overlapping seeds
 * @todo move me to another file
 */
class SeedLumping : public libMS::Module<Seeds, false, Seeds, NucSeq, Pack>
{
  public:
    /*
     * currently seeds on the reverse strand are pointing to the top left instead of the top right.
     * The following algorithm requires all seeds to point to the top right
     * therefore we mirror the reverse strand seeds to the reverse complement strand (and back once were done)
     *
     * We do not need to know the actual size of the reference in order to do that
     * we can just use the maximal possible reference size.
     */
    static int64_t getDelta( const Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return rSeed.start_ref( ) - rSeed.start( );
        else
            return rSeed.start_ref( ) + rSeed.start( );
    } // method

  private:
    template <typename Func1_t, typename Func2_t>
    std::shared_ptr<Seeds> execute_helper( std::shared_ptr<Seeds> pIn, Func1_t&& fExtendSeedFunc,
                                           Func2_t&& fExtendSeedRightFunc )
    {
        if( pIn->size( ) == 0 )
            return std::make_shared<Seeds>( );

        std::vector<std::pair<Seed*, int64_t>> vSortedSeeds;
        vSortedSeeds.reserve( pIn->size( ) );
        for( auto& rSeed : *pIn )
            vSortedSeeds.push_back( std::make_pair( &rSeed, getDelta( rSeed ) ) );

        std::sort( //
            vSortedSeeds.begin( ),
            vSortedSeeds.end( ),
            [&]( const std::pair<Seed*, int64_t>& rA, const std::pair<Seed*, int64_t>& rB ) //
            { //
                if( rA.first->bOnForwStrand != rB.first->bOnForwStrand )
                    return rA.first->bOnForwStrand;
                if( rA.second != rB.second )
                    return rA.second < rB.second;
                return rA.first->start( ) < rB.first->start( );
            } // lambda
        ); // sort call

        auto pRet = std::make_shared<Seeds>( );
        pRet->reserve( vSortedSeeds.size( ) );

        Seed* pLast = vSortedSeeds.front( ).first;
        int64_t iLastDelta = vSortedSeeds.front( ).second;

        for( size_t iI = 1; iI < vSortedSeeds.size( ); iI++ )
        {
            Seed* pSeed = vSortedSeeds[ iI ].first;
            int64_t iCurrDelta = vSortedSeeds[ iI ].second;
            if( pLast->bOnForwStrand == pSeed->bOnForwStrand && iLastDelta == iCurrDelta )
            {
                fExtendSeedRightFunc( *pLast, *pSeed );
                if( pLast->end( ) >= pSeed->start( ) )
                {
                    if( pSeed->end( ) > pLast->end( ) )
                        pLast->iSize = pSeed->end( ) - pLast->start( );
                    assert( pLast->end( ) >= pSeed->end( ) );
                    assert( pLast->end_ref( ) >= pSeed->end_ref( ) );
                    // do not add *pSeed, to pRet since it has been merged with pLast.
                    continue;
                } // if
            } // if

            fExtendSeedFunc( *pLast );
            pRet->push_back( *pLast );

            iLastDelta = iCurrDelta;
            pLast = pSeed;
        } // for

        // add remaining seed
        fExtendSeedFunc( *pLast );
        pRet->push_back( *pLast );

        return pRet;
    }

  public:
    SeedLumping( const ParameterSetManager& rParameters )
    {} // default constructor
    SeedLumping( )
    {} // default constructor


    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pIn, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRef )
    {
        return execute_helper(
            pIn, [&]( Seed& rSeed ) { SeedExtender::extendSeed( rSeed, pQuery, pRef ); },
            [&]( Seed& rLast, Seed& rSeed ) {
                size_t uiBackw = 0;
                if( rLast.bOnForwStrand )
                    while( rLast.end( ) + uiBackw < rSeed.start( ) && pQuery->pxSequenceRef[ rLast.end( ) + uiBackw ] ==
                                                                          pRef->vExtract( rLast.end_ref( ) + uiBackw ) )
                        uiBackw++;
                else
                    while( rLast.end( ) + uiBackw < rSeed.start( ) &&
                           pQuery->pxSequenceRef[ rLast.end( ) + uiBackw ] ==
                               3 - pRef->vExtract( rLast.start_ref( ) - rLast.size( ) - uiBackw ) )
                        uiBackw++;
                rLast.iSize += uiBackw;
            } );
    } // method

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pIn, NucSeq& rQ1, NucSeq& rQ2 )
    {
        return execute_helper(
            pIn, [&]( Seed& rSeed ) { SeedExtender::extendSeed( rSeed, rQ1, rQ2 ); },
            [&]( Seed& rLast, Seed& rSeed ) {
                size_t uiBackw = 0;
                if( rLast.bOnForwStrand )
                    while( rLast.end( ) + uiBackw < rSeed.start( ) &&
                           rQ1.pxSequenceRef[ rLast.end( ) + uiBackw ] ==
                               rQ2.pxSequenceRef[ rLast.end_ref( ) + uiBackw ] )
                        uiBackw++;
                else
                    while( rLast.end( ) + uiBackw < rSeed.start( ) &&
                           rQ1.pxSequenceRef[ rLast.end( ) + uiBackw ] ==
                               3 - rQ2.pxSequenceRef[ rLast.start_ref( ) - rLast.size( ) - uiBackw ] )
                        uiBackw++;
                rLast.iSize += uiBackw;
            } );
    } // method

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pIn, std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
    {
        return execute( pIn, *pQ1, *pQ2 );
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> lump( std::vector<std::shared_ptr<libMA::Seeds>> vIn,
                                                             std::vector<std::shared_ptr<NucSeq>>
                                                                 vQueries,
                                                             std::shared_ptr<Pack>
                                                                 pRef )
    {
        if( vIn.size( ) != vQueries.size( ) )
            throw std::runtime_error( "vIn and vQueries have different lenghts" );
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        vRet.reserve( vIn.size( ) );
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ], vQueries[ uiI ], pRef ) );
        return vRet;
    } // method
}; // class


/**
 * @brief filters out exact duplicates among seeds by sorting and comparing them.
 * @details
 * This is implemented to show the runtime advantage of removeing duplicates first.
 */
class SortRemoveDuplicates : public libMS::Module<Seeds, false, Seeds>
{
  public:
    SortRemoveDuplicates( const ParameterSetManager& rParameters )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {
        std::sort( pSeeds->begin( ), pSeeds->end( ), []( const Seed& rA, const Seed& rB ) {
            if( rA.bOnForwStrand != rB.bOnForwStrand )
                return rA.bOnForwStrand;
            if( rA.start_ref( ) != rB.start_ref( ) )
                return rA.start_ref( ) < rB.start_ref( );
            if( rA.start( ) != rB.start( ) )
                return rA.start( ) < rB.start( );
            return rA.size( ) < rB.size( );
        } );
        Seed* pLast = nullptr;
        auto pRet = std::make_shared<Seeds>( );
        pRet->reserve( pSeeds->size( ) );
        for( auto& rSeed : *pSeeds )
        {
            if( pLast == nullptr || pLast->start_ref( ) != rSeed.start_ref( ) || pLast->start( ) != rSeed.start( ) ||
                pLast->size( ) != rSeed.size( ) )
                pRet->push_back( rSeed );
            pLast = &rSeed;
        } // for
        return pRet;
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> filter( std::vector<std::shared_ptr<libMA::Seeds>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ] ) );
        return vRet;
    } // method
}; // class

/**
 * @brief Filters a set of seeds until there are only unique seeds left
 * @details
 * allows up to n missmatches
 * implemented very inefficiently at the moment.
 * @ingroup module
 */
class FilterToUnique : public libMS::Module<Seeds, false, Seeds, NucSeq, NucSeq>
{
  public:
    size_t uiNumMissmatchesAllowed = 3;

    FilterToUnique( const ParameterSetManager& rParameters )
    {} // default constructor
    FilterToUnique( )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef )
    {
        auto pRet = std::make_shared<Seeds>( );
        for( auto& rSeed : *pSeeds )
        {
            bool bFoundMatch = false;
            for( nucSeqIndex uiI = 0; uiI < pRef->length( ) - rSeed.size( ); uiI++ )
            {
                size_t uiCurrNumMM = 0;
                if( uiI != rSeed.start_ref( ) )
                {
                    for( nucSeqIndex uiJ = 0; uiJ < rSeed.size( ); uiJ++ )
                        if( ( *pRef )[ uiI + uiJ ] != ( *pQuery )[ rSeed.start( ) + uiJ ] )
                            uiCurrNumMM++;
                    if( uiCurrNumMM <= uiNumMissmatchesAllowed )
                    {
                        bFoundMatch = true;
                        break;
                    } // if
                }
            } // for
            if( !bFoundMatch )
                pRet->push_back( rSeed );
        } // for
        return pRet;
    } // method
}; // class

#if 1
/**
 * @brief Filters a set of maximally extended seeds down to SMEMs
 * @ingroup module
 */
class MaxExtendedToSMEM : public libMS::Module<Seeds, false, Seeds>
{
  public:
    MaxExtendedToSMEM( const ParameterSetManager& rParameters )
    {} // constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {
        std::sort( //
            pSeeds->begin( ),
            pSeeds->end( ),
            []( const Seed& rA, const Seed& rB ) //
            { //
                if( rA.start( ) == rB.start( ) )
                {
                    if( rA.size( ) == rB.size( ) )
                        return rA.start_ref( ) < rB.start_ref( );
                    return rA.size( ) > rB.size( );
                } // if
                return rA.start( ) < rB.start( );
            } // lambda
        ); // sort call

        nucSeqIndex uiMaxSeenPos = 0;

        auto pRet = std::make_shared<Seeds>( );

        for( auto& rSeed : *pSeeds )
        {
            if( rSeed.end( ) > uiMaxSeenPos )
                pRet->push_back( rSeed );
            // allow seeds that have exactly the same query interval, but filter out duplicates here
            else if( rSeed.end( ) == uiMaxSeenPos && rSeed.start( ) == pRet->back( ).start( ) &&
                     rSeed.start_ref( ) != pRet->back( ).start_ref( ) )
                pRet->push_back( rSeed );
            uiMaxSeenPos = std::max( rSeed.end( ), uiMaxSeenPos );
        } // for
        return pRet;
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> filter( std::vector<std::shared_ptr<libMA::Seeds>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ] ) );
        return vRet;
    } // method
}; // class
#endif

/**
 * @brief Filters a set of maximally extended seeds down to SMEMs
 * @ingroup module
 */
class MinLength : public libMS::Module<Seeds, false, Seeds>
{
    size_t uiMinLen;

  public:
    MinLength( const ParameterSetManager& rParameters, size_t uiMinLen ) : uiMinLen( uiMinLen )
    {} // constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {
        pSeeds->erase( std::remove_if( pSeeds->begin( ), pSeeds->end( ),
                                       [&]( const Seed& rSeed ) { return rSeed.size( ) < uiMinLen; } ),
                       pSeeds->end( ) );
        // std::cout << pSeeds->size( ) << std::endl;
        return pSeeds;
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> filter( std::vector<std::shared_ptr<libMA::Seeds>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ] ) );
        return vRet;
    } // method
}; // class


#if 1
/**
 * @brief Filters a set of maximally extended seeds down to MaxSpanning
 * @ingroup module
 */
class MaxExtendedToMaxSpanning : public libMS::Module<Seeds, false, Seeds>
{
    struct SeedSmallerComp
    {
        bool operator( )( const Seed* pA, const Seed* pB )
        {
            if( pA->size( ) != pB->size( ) )
                return pA->size( ) < pB->size( );
            if( pA->start( ) != pB->start( ) )
                return pA->start( ) < pB->start( );
            return pA->start_ref( ) < pB->start_ref( );
        } // operator
    }; // struct

  public: //
    MaxExtendedToMaxSpanning( const ParameterSetManager& rParameters )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {

        interval_tree::IntervalTree<nucSeqIndex, Seed*>::interval_vector vIntervals;
        std::vector<nucSeqIndex> vStartVec;
        for( auto& rSeed : *pSeeds )
        {
            vIntervals.emplace_back( rSeed.start( ), rSeed.end( ) - 1, &rSeed );
            vStartVec.push_back( rSeed.start( ) );
        } // for
        interval_tree::IntervalTree<nucSeqIndex, Seed*> xTree( std::move( vIntervals ) );
        std::sort( vStartVec.begin( ), vStartVec.end( ) );

        auto pRet = std::make_shared<Seeds>( );

        nucSeqIndex uiX = 0;

        while( true )
        {
            std::vector<Seed*> vOverlaps;
            xTree.visit_overlapping( uiX, uiX,
                                     [&]( const interval_tree::IntervalTree<nucSeqIndex, Seed*>::interval& rInterval ) {
                                         vOverlaps.push_back( rInterval.value );
                                     } ); // visit_overlapping call
            if( vOverlaps.size( ) == 0 )
            {
                // check for non overlapping seed to the right
                auto pIt = std::lower_bound( vStartVec.begin( ), vStartVec.end( ), uiX );

                if( pIt == vStartVec.end( ) ) // we still have not found any intervals...
                    break;
                else // we are not at the end of the query yet; there is just a gap between seeds
                    uiX = *pIt; // jump to the end of this gap
                // from here the loop will just repeat so we will use the visit_overlapping call from above
                // to fill vOverlaps instead of rewriting the same code here...
            } // if
            else
            {
                std::make_heap( vOverlaps.begin( ), vOverlaps.end( ), SeedSmallerComp( ) );
                pRet->push_back( *vOverlaps.front( ) );
                uiX = vOverlaps.front( )->end( );
                // next two lines: pop heap
                std::pop_heap( vOverlaps.begin( ), vOverlaps.end( ), SeedSmallerComp( ) );
                vOverlaps.pop_back( );
                while( vOverlaps.size( ) > 0 && vOverlaps.front( )->size( ) == pRet->back( ).size( ) )
                {
                    if( pRet->back( ).start_ref( ) != vOverlaps.front( )->start_ref( ) )
                        pRet->push_back( *vOverlaps.front( ) );
                    // next two lines: pop heap
                    std::pop_heap( vOverlaps.begin( ), vOverlaps.end( ), SeedSmallerComp( ) );
                    vOverlaps.pop_back( );
                } // while
            } // else
        } // while

        return pRet;
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> filter( std::vector<std::shared_ptr<libMA::Seeds>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ] ) );
        return vRet;
    } // method
}; // class
#endif

/**
 * @brief Filters out seeds that are overlapping more than x nt (filters the smaller seed)
 * @ingroup module
 */
class FilterOverlappingSeeds : public libMS::Module<Seeds, false, Seeds>
{
    nucSeqIndex uiMaxOverlap;

  public:
    FilterOverlappingSeeds( const ParameterSetManager& rParameters, nucSeqIndex uiMaxOverlap )
        : uiMaxOverlap( uiMaxOverlap )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {
        std::sort( pSeeds->begin( ), pSeeds->end( ),
                   []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );

        std::shared_ptr<Seeds> pRet;
        pRet->reserve( pSeeds->size( ) );


        for( auto& rS : *pSeeds )
        {
            if( pRet->size( ) == 0 || rS.start( ) > pRet->back( ).end( ) + uiMaxOverlap )
                pRet->push_back( rS );
            else if( pRet->back( ).size( ) > rS.size( ) )
            {
                pRet->pop_back( );
                pRet->push_back( rS );
            }
        } // for

        return pRet;
    } // method

    virtual std::vector<std::shared_ptr<libMA::Seeds>> filter( std::vector<std::shared_ptr<libMA::Seeds>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( size_t uiI = 0; uiI < vIn.size( ); uiI++ )
            vRet.push_back( execute( vIn[ uiI ] ) );
        return vRet;
    } // method
}; // class

/**
 * @brief Filters a set of seeds removing all parlindrome seeds
 * @details
 * parlindrome seeds are two seeds that are crossing and on opposite strands.
 * @ingroup module
 */
class ParlindromeFilter : public libMS::Module<Seeds, false, Seeds>
{
    static int64_t startX( Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return rSeed.start_ref( );
        else
            return rSeed.start_ref( ) - rSeed.size( ) + 1;
    } // method

    static int64_t endX( Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return rSeed.end_ref( ) - 1;
        else
            return rSeed.start_ref( );
    } // method

    static int64_t startY( Seed& rSeed )
    {
        return rSeed.start( );
    } // method

    static int64_t endY( Seed& rSeed )
    {
        return rSeed.end( ) - 1;
    } // method

    static int64_t rotatedStartX( Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return startX( rSeed ) + startY( rSeed );
        else
            return endX( rSeed ) + startY( rSeed );
    } // method

    static int64_t rotatedEndX( Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return endX( rSeed ) + endY( rSeed );
        else
            return rotatedStartX( rSeed );
    } // method

    static int64_t rotatedStartY( Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return startX( rSeed ) - startY( rSeed );
        else
            return startX( rSeed ) - endY( rSeed );
    } // method

    static int64_t rotatedEndY( Seed& rSeed )
    {
        if( rSeed.bOnForwStrand )
            return rotatedStartY( rSeed );
        else
            return endX( rSeed ) - startY( rSeed );
    } // method

    /**
     * @brief sorting order for the linesweep
     * @details
     * Sorts after the rotatedStartX of seeds.
     * Also sorts all forward strand seeds before the reverse strand seeds.
     */
    static bool seedIsSmallerSort( Seed& rA, Seed& rB )
    {
        if( rA.bOnForwStrand == rB.bOnForwStrand )
            return rotatedStartX( rA ) < rotatedStartX( rB );
        return rA.bOnForwStrand;
    } // method

    /**
     * @brief sorting order for the heap
     * @details
     * Sorts after the rotatedEndX of seeds.
     */
    static bool seedIsSmallerHeap( Seed* pA, Seed* pB )
    {
        return rotatedEndX( *pA ) < rotatedEndX( *pB );
    } // method

    /**
     * @brief sorting order for the set
     * @details
     * sorts after the rotatedStartY.
     * for equal rotatedStartY seeds are sorter after their memory positions (so that we can use std::find).
     */
    static bool seedIsSmallerMap( Seed* pA, Seed* pB )
    {
        if( rotatedStartY( *pA ) == rotatedStartY( *pB ) && pA->bOnForwStrand == pB->bOnForwStrand )
        {
            return *pA < *pB;
        } // if
        return rotatedStartY( *pA ) < rotatedStartY( *pB );
    } // method

  public:
    std::shared_ptr<Seeds> pParlindromes = nullptr;

    ParlindromeFilter( const ParameterSetManager& rParameters )
    {} // default constructor

    /**
     * @brief call this to keep all parlindrom seeds
     * @details
     * calling this will discard all currently kept parlindromes and set up a datastructor to keep all future ones.
     * @note Once this is called the module is no longer threadsave!
     */
    void keepParlindromes( )
    {
        pParlindromes = std::make_shared<Seeds>( );
    } // method

/**
 * @brief
 * removes all parlindrom seeds
 * @details
 * Seeds that indicate parlindromes must be of opposite strands and cross each other.
 * This is a linesweep:
 * all seeds are rotated rightwards 45 degrees:
 *  - forward strand seeds are then horizontal lines
 *  - reverse strand seeds are then vertical lines
 * we iterate over the rotated x-axis (i.e. sort seeds by their start position on the x-axis). There are 3 cases
 * depending on wether we stop at the start of a forward strand seed, end of a forward strand seed, or position of a
 * reverse strand seed (vertical line: start & end are the same)
 *  - start of forward strand seed: remember that the seed is "open" and that we have not seen a longer crossing
 *    seed
 *  - end of forward strand seed: remove the seed form the set of "open" seeds; if we have never found a longer
 *    crossing seed we can append this seed to the output list
 *  - reverse strand seed: check the "open" seeds for crossing seeds; If this seed has no longer crossing seed
 *    append it directly to the output list. Set the marker for shorter crossing seeds in the "open" list
 *    accordingly.
 * The "open" set is implemented as a binary search tree ordered by the y-position of open seeds (forward strand
 * seeds are horizontal lines -> only one y-position). This allows for the efficient search of overlapping seeds on
 * the y-axis.
 * The end of forward strand seed case is implemented via a heap that holds all open seeds ordered by their
 * endpoints on the x-axis.
 *
 * @todo O( n^2 ) with respect to the amount of parlindromes (can this be improved?)
 *       O( n log n ) with respect to the amount of seeds.
 */
#define DEBUG_CODE 0
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {
        auto pRet = std::make_shared<Seeds>( );
        std::sort( pSeeds->begin( ), pSeeds->end( ), &seedIsSmallerSort );
#if DEBUG_CODE
        for( auto& rSeed : *pSeeds )
            std::cout << ( rSeed.bOnForwStrand ? "f " : "r " ) << rSeed.start( ) << "," << rSeed.start_ref( ) << ","
                      << rSeed.size( ) << " " << rotatedStartX( rSeed ) << " - " << rotatedEndX( rSeed ) << "   "
                      << rotatedStartY( rSeed ) << " - " << rotatedEndY( rSeed ) << std::endl;
        std::cout << std::endl;
#endif
        // forward strand seeds are at the start of the array
        auto xForwStartIt = pSeeds->begin( );
        // reverse strand seeds are at the end of the array
        auto xRevIt = std::lower_bound( pSeeds->begin( ), pSeeds->end( ), 0,
                                        []( Seed& rSeed, int ) { return rSeed.bOnForwStrand; } );
        // check if there are no forward strand or no reverse strand seeds
        if( xRevIt == pSeeds->begin( ) || xRevIt == pSeeds->end( ) )
        {
#if DEBUG_CODE
            std::cout << "XXX" << std::endl;
#endif
            // if so just return all seeds
            return pSeeds;
        }
        assert( !xRevIt->bOnForwStrand );
        assert( xForwStartIt->bOnForwStrand );
        /*
         * This binary search tree holds the open seeds so that we can perform a lookup over them.
         * All seeds in here must be on the forward strand.
         * Therefore, these seeds have only one rotated y position (i.e. they are a horizontal line)
         * and so we can use a simple set as search datastructure.
         * The boolean value holds wether we have found a longer crossing seed or not.
         */
        std::map<Seed*, bool, decltype( &seedIsSmallerMap )> xIntervals( &seedIsSmallerMap );
        /**
         * holds the open seeds in order of their endpoints on the rotated x-axis, so that we know when to remove them
         * form xIntervals.
         */
        std::vector<Seed*> xHeap;
        /*
         * Iterate with the iterators (xForwStartIt and xRevIt) and heap xHeap at the same time.
         * The line sweep is over the rotated x-axis.
         * - xForwStartIt: inserts a new seed into xIntervals once we pass the start of a forward strand seed on x-axis
         * - xHeap.front(): removes a seed from xIntervals once we pass the end of such a forward strand seed on x-axis
         * - xRevIt: checks for parlindrome seeds via xIntervals (only iterator on reverse strand seeds)
         */
        while( xRevIt != pSeeds->end( ) )
        {
            // if xForwStartIt points to the seed that appears first on the rotated x-axis
            if( xForwStartIt->bOnForwStrand &&
                ( xHeap.size( ) == 0 || rotatedStartX( *xForwStartIt ) <= rotatedEndX( *xHeap.front( ) ) ) &&
                rotatedStartX( *xForwStartIt ) <= rotatedStartX( *xRevIt ) )
            {
                // check if there was a previous seed
                if( xForwStartIt != pSeeds->begin( ) )
                {
                    // the seed pointed to by xForwStartIt might be equal to a seed that we already passed
                    if( *xForwStartIt == *( xForwStartIt - 1 ) )
                    {
                        // if so, ignore this duplicate
                        xForwStartIt++;
                        continue;
                    } // if
                } // if
#if DEBUG_CODE
                std::cout << "AAA " << ( xForwStartIt->bOnForwStrand ? "f " : "r " ) << xForwStartIt->start( ) << ","
                          << xForwStartIt->start_ref( ) << "," << xForwStartIt->size( ) << std::endl;
#endif
                assert( xIntervals.find( &*xForwStartIt ) == xIntervals.end( ) );
                xIntervals[&*xForwStartIt ] = true; // insert
                xHeap.push_back( &*xForwStartIt );
                std::push_heap( xHeap.begin( ), xHeap.end( ), &seedIsSmallerHeap );
                xForwStartIt++;
            } // if
            // if xRevIt points to the seed that appears first on the rotated x-axis
            else if( xHeap.size( ) == 0 || rotatedStartX( *xRevIt ) <= rotatedEndX( *xHeap.front( ) ) )
            {
#if DEBUG_CODE
                std::cout << "BBB " << ( xRevIt->bOnForwStrand ? "f " : "r " ) << xRevIt->start( ) << ","
                          << xRevIt->start_ref( ) << "," << xRevIt->size( ) << std::endl;
                std::cout << xIntervals.size( ) << std::endl;
#endif
                auto xCheckIt = xIntervals.lower_bound( &*xRevIt );
                bool bRevItIsLargest = true;
                while( xCheckIt != xIntervals.end( ) && rotatedEndY( *xCheckIt->first ) <= rotatedEndY( *xRevIt ) )
                {
#if DEBUG_CODE
                    std::cout << "DDD " << ( xCheckIt->first->bOnForwStrand ? "f " : "r " ) << xCheckIt->first->start( )
                              << "," << xCheckIt->first->start_ref( ) << "," << xCheckIt->first->size( ) << std::endl;
#endif
                    if( xRevIt->size( ) <= xCheckIt->first->size( ) )
                        bRevItIsLargest = false;
                    if( xCheckIt->first->size( ) <= xRevIt->size( ) )
                        xCheckIt->second = false;
                    xCheckIt++;
                } // while
#if DEBUG_CODE
                std::cout << ( bRevItIsLargest ? "true " : "false " ) << std::endl;
#endif
                if( bRevItIsLargest )
                    pRet->push_back( *xRevIt );
                else if( pParlindromes != nullptr )
                    pParlindromes->push_back( *xRevIt );
                xRevIt++;
            } // else
            // if xHeap.front() points to the seed that appears first on the rotated x-axis
            else
            {
#if DEBUG_CODE
                std::cout << "CCC " << ( xHeap.front( )->bOnForwStrand ? "f " : "r " ) << xHeap.front( )->start( )
                          << "," << xHeap.front( )->start_ref( ) << "," << xHeap.front( )->size( ) << std::endl;
#endif
                auto xFind = xIntervals.find( xHeap.front( ) );
                assert( xFind != xIntervals.end( ) );
                if( xFind->second )
                    pRet->push_back( *xHeap.front( ) );
                else if( pParlindromes != nullptr )
                    pParlindromes->push_back( *xHeap.front( ) );
                xIntervals.erase( xFind );
                std::pop_heap( xHeap.begin( ), xHeap.end( ), &seedIsSmallerHeap );
                xHeap.pop_back( );
            } // if
        } // while
        // add the forward strand seeds we have not visited yet to pRet
        while( xForwStartIt->bOnForwStrand )
        {
            pRet->push_back( *xForwStartIt );
            xForwStartIt++;
        } // while
        // add the seeds in xIntervals that are no parlindromes to pRet
        for( auto xPair : xIntervals )
        {
            if( xPair.second )
                pRet->push_back( *xPair.first );
            else if( pParlindromes != nullptr )
                pParlindromes->push_back( *xPair.first );
        } // for

#if DEBUG_CODE
        size_t uiCntForw = 0;
        size_t uiCntRev = 0;
        for( auto& rSeed : *pSeeds )
            if( rSeed.bOnForwStrand )
                uiCntForw++;
            else
                uiCntRev++;
        std::cout << std::flush;
        if( uiCntForw == 11 && uiCntRev == 8 )
            raise( SIGINT );
#endif

        return pRet;
    } // method
}; // class


} // namespace libMA
