/**
 * @file seedFilters.h
 * @brief Implements several filters for seeds
 * @author Markus Schmidt
 */
#pragma once


#include "IntervalTree.h"
#include "stripOfConsideration.h"

namespace libMA
{


struct SeedCmp
{
    bool operator( )( const Seed& rA, const Seed& rB ) const
    {
        if( rA.start( ) != rB.start( ) )
            return rA.start( ) < rB.start( );

        if( rA.start_ref( ) != rB.start_ref( ) )
            return rA.start_ref( ) < rB.start_ref( );

        if( rA.bOnForwStrand != rB.bOnForwStrand )
            return rA.bOnForwStrand;

        return rA.size( ) < rB.size( );
    }
}; // struct

/// @brief this class exists mereley to expose the return value of execute_helper_py to python
class HelperRetVal
{
  public:
    std::shared_ptr<Seeds> pSeeds;
    std::shared_ptr<Seeds> pRemovedSeeds;
    std::vector<size_t> vLayerOfSeeds;
    std::vector<bool> vParlindromeSeed;
    std::vector<bool> vOverlappingSeed;
    size_t uiCurrSocID = 0;
    std::vector<size_t> vSocIds;
    std::vector<geom::Rectangle<nucSeqIndex>> vRectangles;
    std::vector<size_t> vRectangleLayers;
    std::vector<double> vRectangleFillPercentage;
    std::vector<size_t> vRectangleReferenceAmbiguity;
    std::vector<size_t> vRectangleKMerSize;
    std::vector<bool> vRectangleUsedDp;
    std::vector<std::pair<Seed, Seed>> vJumpSeeds;

    HelperRetVal( ) : pSeeds( std::make_shared<Seeds>( ) ), pRemovedSeeds( std::make_shared<Seeds>( ) ){ };

    inline void computeOverlappingSeedsVector( std::shared_ptr<Seeds> pFilteredRet )
    {
        this->vOverlappingSeed.clear( );
        SeedCmp xSort;
        std::sort( pFilteredRet->begin( ), pFilteredRet->end( ), xSort );
        for( auto& rSeed : *this->pSeeds )
            this->vOverlappingSeed.push_back(
                !std::binary_search( pFilteredRet->begin( ), pFilteredRet->end( ), rSeed, xSort ) );
    } // method
}; // class

#if 1
/**
 * @brief Extends seeds at the end to create maximally extened seeds
 * @ingroup module
 */
class SeedExtender : public libMS::Module<Seeds, false, Seeds, NucSeq, Pack>
{
    template <typename Func1_t, typename Func2_t>
    static void extendSeedHelper( Seed& rSeed, Func1_t&& fGetNucQuery, Func2_t&& fGetNucRef, size_t uiQSize,
                                  size_t uiRStart, size_t uiREnd )
    {
        size_t uiForw = 1;
        if( rSeed.bOnForwStrand )
            while( uiForw <= rSeed.start( ) && //
                   uiForw + uiRStart <= rSeed.start_ref( ) && //
                   fGetNucQuery( rSeed.start( ) - uiForw ) == fGetNucRef( rSeed.start_ref( ) - uiForw ) )
                uiForw++;
        else
            while( uiForw <= rSeed.start( ) && //
                   uiForw + rSeed.start_ref( ) < uiREnd && //
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
                   uiBackw + rSeed.end_ref( ) < uiREnd && //
                   fGetNucQuery( rSeed.end( ) + uiBackw ) == fGetNucRef( rSeed.end_ref( ) + uiBackw ) )
                uiBackw++;
        else
            while( uiBackw + rSeed.end( ) < uiQSize && //
                   rSeed.start_ref( ) >= uiRStart + uiBackw + rSeed.size( ) && //
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
        auto uiContigID = pRef->uiSequenceIdForPosition( rSeed.start_ref( ) );
        extendSeedHelper(
            rSeed, [ & ]( size_t uiI ) { return pQuery->pxSequenceRef[ uiI ]; },
            [ & ]( size_t uiI ) { return pRef->vExtract( uiI ); }, pQuery->length( ),
            pRef->startOfSequenceWithId( uiContigID ), pRef->endOfSequenceWithId( uiContigID ) );
    }

    static void extendSeed( Seed& rSeed, std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
    {
        extendSeedHelper(
            rSeed, [ & ]( size_t uiI ) { return pQ1->pxSequenceRef[ uiI ]; },
            [ & ]( size_t uiI ) { return pQ2->pxSequenceRef[ uiI ]; }, pQ1->length( ), 0, pQ2->length( ) );
    }

    static void extendSeed( Seed& rSeed, NucSeq& rQ1, NucSeq& rQ2 )
    {
        extendSeedHelper(
            rSeed, [ & ]( size_t uiI ) { return rQ1.pxSequenceRef[ uiI ]; },
            [ & ]( size_t uiI ) { return rQ2.pxSequenceRef[ uiI ]; }, rQ1.length( ), 0, rQ2.length( ) );
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
            [ & ]( const std::pair<Seed*, int64_t>& rA, const std::pair<Seed*, int64_t>& rB ) //
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
        auto pRet = execute_helper(
            pIn, [ & ]( Seed& rSeed ) { SeedExtender::extendSeed( rSeed, pQuery, pRef ); },
            [ & ]( Seed& rLast, Seed& rSeed ) {
                auto uiContigID = pRef->uiSequenceIdForPosition( rLast.start_ref( ) );

                size_t uiBackw = 0;
                if( rLast.bOnForwStrand )
                    while( rLast.end( ) + uiBackw < rSeed.start( ) &&
                           rLast.end( ) >= uiBackw + pRef->startOfSequenceWithId( uiContigID ) &&
                           pQuery->pxSequenceRef[ rLast.end( ) + uiBackw ] ==
                               pRef->vExtract( rLast.end_ref( ) + uiBackw ) )
                        uiBackw++;
                else
                    while( rLast.end( ) + uiBackw < rSeed.start( ) &&
                           rLast.end( ) + uiBackw < pRef->endOfSequenceWithId( uiContigID ) &&
                           pQuery->pxSequenceRef[ rLast.end( ) + uiBackw ] ==
                               3 - pRef->vExtract( rLast.start_ref( ) - rLast.size( ) - uiBackw ) )
                        uiBackw++;
                rLast.iSize += uiBackw;
            } );
        return pRet;
    } // method

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pIn, NucSeq& rQ1, NucSeq& rQ2 )
    {
        return execute_helper(
            pIn, [ & ]( Seed& rSeed ) { SeedExtender::extendSeed( rSeed, rQ1, rQ2 ); },
            [ & ]( Seed& rLast, Seed& rSeed ) {
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

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute_py( std::shared_ptr<Seeds> pIn, std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2 )
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

/**
 * @brief Filters a set of seeds according to the seeds reference positions
 * @details
 * assumes seeds are not bridging
 * @ingroup module
 */
class FilterContigBorder : public libMS::Module<Seeds, false, Seeds, Pack>
{
  public:
    size_t uiMaxDist = 25000;

    FilterContigBorder( const ParameterSetManager& rParameters )
    {} // default constructor
    FilterContigBorder( )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<Pack> pPack )
    {
        auto pRet = std::make_shared<Seeds>( );
        for( auto& rSeed : *pSeeds )
        {
            auto uiStart = rSeed.bOnForwStrand ? rSeed.start_ref( ) : rSeed.start_ref( ) - rSeed.size( ) - 1;
            auto uiEnd = rSeed.bOnForwStrand ? rSeed.end_ref( ) : rSeed.start_ref( ) - 1;
            auto uiIdx = pPack->uiSequenceIdForPosition( uiStart );
            // sanity check
            if( pPack->uiSequenceIdForPosition( uiEnd ) != uiIdx )
                continue;
            if( pPack->startOfSequenceWithId( uiIdx ) + uiMaxDist >= uiStart )
                continue;
            if( pPack->endOfSequenceWithId( uiIdx ) <= uiEnd + uiMaxDist )
                continue;
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
                                       [ & ]( const Seed& rSeed ) { return rSeed.size( ) < uiMinLen; } ),
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
            xTree.visit_overlapping(
                uiX, uiX, [ & ]( const interval_tree::IntervalTree<nucSeqIndex, Seed*>::interval& rInterval ) {
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
 * @brief Filters out seeds that are overlapping
 * @details
 * shortens / removes seeds until no overlaps on query are remaining
 * shortens / removes both overlapping seeds -> there will be gaps after this
 * @ingroup module
 */
class FilterOverlappingSeeds : public libMS::Module<Seeds, false, Seeds>
{
    nucSeqIndex uiMinNtNonOverlap = 16; // or if there are 16 non overlapping NT

  public:
    FilterOverlappingSeeds( const ParameterSetManager& rParameters )
    {} // default constructor

    // overload
    std::shared_ptr<Seeds> execute_helper( std::shared_ptr<Seeds> pSeeds, HelperRetVal* pOutExtra )
    {
        std::sort( pSeeds->begin( ), pSeeds->end( ), []( const Seed& rA, const Seed& rB ) {
            if( rA.start( ) == rB.start( ) )
                return rA.size( ) > rB.size( );
            return rA.start( ) < rB.start( );
        } );

        auto pRet = std::make_shared<Seeds>( );
        pRet->reserve( pSeeds->size( ) ); // allocate all space at once


        nucSeqIndex uiMax = 0;
        for( size_t uiI = 0; uiI < pSeeds->size( ); uiI++ )
        {
            size_t uiJ = uiI + 1;
            nucSeqIndex uiLocalMax = std::max( uiMax, ( *pSeeds )[ uiI ].start( ) );
            while( uiJ <= pSeeds->size( ) && uiLocalMax < ( *pSeeds )[ uiI ].end( ) )
            {
                auto uiLocalEnd = ( *pSeeds )[ uiI ].end( );
                if( uiJ < pSeeds->size( ) && ( *pSeeds )[ uiJ ].start( ) < uiLocalEnd )
                    uiLocalEnd = ( *pSeeds )[ uiJ ].start( );

                // only keep this non overlapping section of the seed if it is at least uiMinNtNonOverlap long
                // or if it is the entire seed
                if( uiLocalMax + uiMinNtNonOverlap < uiLocalEnd ||
                    ( uiLocalMax == ( *pSeeds )[ uiI ].start( ) && uiLocalEnd == ( *pSeeds )[ uiI ].end( ) ) )
                {
                    auto uiLen = uiLocalEnd - uiLocalMax;
                    auto uiRefPos = ( *pSeeds )[ uiI ].start_ref( );
                    if( ( *pSeeds )[ uiI ].bOnForwStrand )
                        uiRefPos += uiLocalMax - ( *pSeeds )[ uiI ].start( );
                    else
                        uiRefPos -= uiLocalMax - ( *pSeeds )[ uiI ].start( );
                    pRet->emplace_back( uiLocalMax, uiLen, uiRefPos, ( *pSeeds )[ uiI ].bOnForwStrand );
                    pRet->back( ).uiAmbiguity = ( *pSeeds )[ uiI ].uiAmbiguity;
                    pRet->back( ).uiSoCNt = ( *pSeeds )[ uiI ].uiSoCNt;
                    pRet->back( ).uiDelta = ( *pSeeds )[ uiI ].uiDelta;
                }
                if( uiJ < pSeeds->size( ) )
                    uiLocalMax = std::max( uiLocalMax, ( *pSeeds )[ uiJ ].end( ) );
                uiJ++;
            } // while
            // record seed as removed (even if it was just broken into segments)
            // if( pOutExtra != nullptr && !( pRet->back( ) == ( *pSeeds )[ uiI ] ) )
            //    pOutExtra->pRemovedSeeds->push_back( ( *pSeeds )[ uiI ] );

            uiMax = std::max( uiMax, ( *pSeeds )[ uiI ].end( ) );
        } // for

        return pRet;
    } // method
    virtual std::shared_ptr<Seeds> DLL_PORT( MA ) execute( std::shared_ptr<Seeds> pSeeds )
    {
        return execute_helper( pSeeds, nullptr );
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
 * @brief Filters out seeds that are overlapping via SoCs
 * @ingroup module
 */
template <bool WITH_SOC> class FilterOverlappingSoCs : public libMS::Module<Seeds, false, Seeds, NucSeq, Pack>
{
    float fMinNonOverlap = 0.50;
    nucSeqIndex uiMinNonOverlap = 5;
    float fValueFac = 2;
    StripOfConsiderationSeeds xSoCs;
    GetAllFeasibleSoCsAsSet xExtractSoCs;

    static void adjustSeed( nucSeqIndex uiFrom, nucSeqIndex uiTo, Seed& xSeed, std::shared_ptr<NucSeq> pQuerySeq,
                            std::shared_ptr<Pack> pRefSeq, HelperRetVal* pOutExtra )
    {
        if( xSeed.start( ) >= uiFrom )
        {
            if( xSeed.start( ) < uiTo )
            {
                // seed actually overlaps region
                if( xSeed.end( ) <= uiTo )
                {
                    // seed completely enclosed
                    if( pOutExtra != nullptr )
                        pOutExtra->pRemovedSeeds->push_back( xSeed );
                    xSeed.iSize = 0;
                } // if
                else
                {
                    // seed's start overlaps region
                    auto uiShortenBy = uiTo - xSeed.start( );
                    if( pOutExtra != nullptr )
                        pOutExtra->pRemovedSeeds->emplace_back( xSeed.start( ), uiShortenBy, xSeed.uiPosOnReference,
                                                                xSeed.bOnForwStrand );
                    // adjust start pos of seed
                    xSeed.iSize -= uiShortenBy;
                    xSeed.iStart += uiShortenBy;
                    if( xSeed.bOnForwStrand ) // take care of orientation of rev strand seeds
                        xSeed.uiPosOnReference += uiShortenBy;
                    else
                        xSeed.uiPosOnReference -= uiShortenBy;
                    libMA::ExtractSeeds::setDeltaOfSeed( xSeed, pQuerySeq->length( ), *pRefSeq, true );
                } // else
            } // if
            // else do nothing (seed does not overlap region)
        } // if
        else
        {
            if( xSeed.end( ) > uiFrom )
            {
                // seed actually overlaps region
                if( xSeed.end( ) <= uiTo )
                {
                    // seed's end overlaps region
                    auto uiShortenBy = xSeed.end( ) - uiFrom;
                    if( pOutExtra != nullptr )
                    {
                        auto uiR = xSeed.uiPosOnReference;
                        if( xSeed.bOnForwStrand ) // take care of orientation of rev strand seeds
                            uiR += uiFrom - xSeed.start( );
                        else
                            uiR -= uiFrom - xSeed.start( );
                        pOutExtra->pRemovedSeeds->emplace_back( uiFrom, uiShortenBy, uiR, xSeed.bOnForwStrand );
                    }
                    assert( uiShortenBy < xSeed.iSize );
                    xSeed.iSize -= uiShortenBy;
                    libMA::ExtractSeeds::setDeltaOfSeed( xSeed, pQuerySeq->length( ), *pRefSeq, true );
                } // if
                else
                {
                    // region cuts seed in half (seed encloses region)
                    if( pOutExtra != nullptr )
                        pOutExtra->pRemovedSeeds->push_back( xSeed );
                    xSeed.iSize = 0;
                } // else
            } // if
            // else do noting (seed does not overlap region)
        } // else
    } // method

    static void removeSeedsInRange( nucSeqIndex uiFrom, nucSeqIndex uiTo, std::shared_ptr<Seeds> pSeeds,
                                    HelperRetVal* pOutExtra, std::shared_ptr<NucSeq> pQuerySeq,
                                    std::shared_ptr<Pack> pRefSeq )
    {
        // shorten overlapping seeds
        for( auto& rSeed : *pSeeds )
            adjustSeed( uiFrom, uiTo, rSeed, pQuerySeq, pRefSeq, pOutExtra );
        // remove seeds that became size 0
        pSeeds->erase(
            std::remove_if( pSeeds->begin( ), pSeeds->end( ), []( const auto& rS ) { return rS.size( ) == 0; } ),
            pSeeds->end( ) );
    } // method

    static void removeSeeds( std::shared_ptr<Seeds> pSeeds, HelperRetVal* pOutExtra )
    {
        if( pOutExtra != nullptr )
            pOutExtra->pRemovedSeeds->append( pSeeds );
        pSeeds->clear( );
    } // method

    static nucSeqIndex
    valueInRangeSeeds( nucSeqIndex uiFrom, nucSeqIndex uiTo,
                       std::tuple<nucSeqIndex, nucSeqIndex, std::shared_ptr<Seeds>, std::shared_ptr<Seeds>> xTuple )
    {
        nucSeqIndex uiValue = 0;
        for( auto& rS : *std::get<3>( xTuple ) )
            if( rS.end( ) > uiFrom && rS.start( ) < uiTo )
                uiValue += std::min( rS.end( ), uiTo ) - std::max( rS.start( ), uiFrom );
        return uiValue;
    } // method

    static void removeSeedInRange( nucSeqIndex uiFrom, nucSeqIndex uiTo, Seed& xSeed, HelperRetVal* pOutExtra,
                                   std::shared_ptr<NucSeq> pQuerySeq, std::shared_ptr<Pack> pRefSeq )
    {
        adjustSeed( uiFrom, uiTo, xSeed, pQuerySeq, pRefSeq, pOutExtra );
    } // method

    static void removeSeed( Seed& xSeed, HelperRetVal* pOutExtra )
    {
        if( pOutExtra != nullptr )
            pOutExtra->pRemovedSeeds->push_back( xSeed );
        xSeed.iSize = 0;
    } // method

    static nucSeqIndex valueInRangeSeed( nucSeqIndex uiFrom, nucSeqIndex uiTo,
                                         std::tuple<nucSeqIndex, nucSeqIndex, Seed, Seed> xTuple )
    {
        if( std::get<3>( xTuple ).end( ) > uiFrom && std::get<3>( xTuple ).start( ) < uiTo )
            return std::min( std::get<3>( xTuple ).end( ), uiTo ) - std::max( std::get<3>( xTuple ).start( ), uiFrom );
        return 0;
    } // method


  public:
    FilterOverlappingSoCs( const ParameterSetManager& rParameters ) : xSoCs( rParameters ), xExtractSoCs( rParameters )
    {
        xSoCs.uiCurrHarmScoreMin = 1;
        xSoCs.fGiveUp = 0;
        xExtractSoCs.uiMinNt = 1; // get All SoCs
    } // default constructor


    template <typename LineSweepVec, typename RemoveSeedsInRange, typename RemoveSeeds, typename ValueInRange>
    void core( LineSweepVec& vLineSweepVec, RemoveSeedsInRange&& fRemoveSeedsInRange, RemoveSeeds&& fRemoveSeeds,
               HelperRetVal* pOutExtra, ValueInRange&& fValueInRange, std::shared_ptr<NucSeq> pQuerySeq,
               std::shared_ptr<Pack> pRefSeq )
    {
        // sort by query start
        std::sort( vLineSweepVec.begin( ), vLineSweepVec.end( ), []( const auto& rA, const auto& rB ) {
            if( std::get<0>( rA ) == std::get<0>( rB ) )
                return std::get<1>( rA ) > std::get<1>( rB );
            return std::get<0>( rA ) < std::get<0>( rB );
        } );

        // line sweep
        for( size_t uiI = 0; uiI < vLineSweepVec.size( ); uiI++ )
        {
            auto uiIStart = std::get<0>( vLineSweepVec[ uiI ] );
            auto uiIEnd = std::get<1>( vLineSweepVec[ uiI ] );
            auto uiPercentageOfI = std::max( ( nucSeqIndex )( ( uiIEnd - uiIEnd ) * fMinNonOverlap ), uiMinNonOverlap );
            size_t uiJ = uiI + 1;
            while( uiJ < vLineSweepVec.size( ) && uiIEnd > std::get<0>( vLineSweepVec[ uiJ ] ) )
            {
                auto uiJStart = std::get<0>( vLineSweepVec[ uiJ ] );
                auto uiJEnd = std::get<1>( vLineSweepVec[ uiJ ] );
                auto uiPercentageOfJ =
                    std::max( ( nucSeqIndex )( ( uiJEnd - uiJStart ) * fMinNonOverlap ), uiMinNonOverlap );

                bool bStartOfIUncovered = uiIStart + uiPercentageOfI < uiJStart;
                bool bEndOfIUncovered = uiJEnd + uiPercentageOfI < uiIEnd;
                bool bEndOfJUncovered = uiIEnd + uiPercentageOfJ < uiJEnd;

                if( bStartOfIUncovered && bEndOfIUncovered )
                {
                    // I encloses J
                    auto uiIValue = fValueInRange( uiJStart, uiJEnd, vLineSweepVec[ uiI ] );
                    auto uiJValue = fValueInRange( uiJStart, uiJEnd, vLineSweepVec[ uiJ ] );
                    if( uiJValue > uiIValue * fValueFac )
                        // j is enclosed by i (and more valuable than i) -> remove seeds in i that overlap j
                        fRemoveSeedsInRange( uiJStart, uiJEnd, std::get<2>( vLineSweepVec[ uiI ] ), pOutExtra,
                                             pQuerySeq, pRefSeq );
                    else
                        // j is enclosed by i (but less valuable than i) -> remove j completely
                        fRemoveSeeds( std::get<2>( vLineSweepVec[ uiJ ] ), pOutExtra );
                } // if
                else if( bStartOfIUncovered && bEndOfJUncovered )
                {
                    // both SoCs have unique regions
                    // make a cut in the center of the overlapping region
                    // and remove seed from one side of the cut in each SoC
                    auto uiCenter = ( uiIEnd + uiJStart ) / 2;
                    fRemoveSeedsInRange( uiCenter, uiIEnd, std::get<2>( vLineSweepVec[ uiI ] ), pOutExtra, pQuerySeq,
                                         pRefSeq );
                    fRemoveSeedsInRange( uiJStart, uiCenter, std::get<2>( vLineSweepVec[ uiJ ] ), pOutExtra, pQuerySeq,
                                         pRefSeq );
                } // if
                else
                {
                    // SoCs cover roughly the same area
                    auto uiIValue = fValueInRange( std::max( uiIStart, uiJStart ), std::min( uiIEnd, uiJEnd ),
                                                   vLineSweepVec[ uiI ] );
                    auto uiJValue = fValueInRange( std::max( uiIStart, uiJStart ), std::min( uiIEnd, uiJEnd ),
                                                   vLineSweepVec[ uiJ ] );
                    // the SoCs are completely overlapping
                    // eliminiate both or keep one if it's value is larger than the other by fValueFac
                    if( uiIValue <= uiJValue * fValueFac )
                        fRemoveSeeds( std::get<2>( vLineSweepVec[ uiI ] ), pOutExtra );
                    if( uiJValue <= uiIValue * fValueFac )
                        fRemoveSeeds( std::get<2>( vLineSweepVec[ uiJ ] ), pOutExtra );
                } // if
                uiJ += 1;
            } // while
        } // for
    }

    // overload
    inline std::shared_ptr<Seeds> execute_helper( std::shared_ptr<Seeds> pSeeds,
                                                  std::shared_ptr<NucSeq>
                                                      pQuerySeq,
                                                  std::shared_ptr<Pack>
                                                      pRefSeq,
                                                  HelperRetVal* pOutExtra )
    {
        if( WITH_SOC )
        {
            auto pSoCQueue = xSoCs.execute( pSeeds, pQuerySeq, pRefSeq );
            auto pSoCs = xExtractSoCs.execute( pSoCQueue );

            //                     query start, query end  , seeds in soc          , seeds in soc (copy for value calc)
            std::vector<std::tuple<nucSeqIndex, nucSeqIndex, std::shared_ptr<Seeds>, std::shared_ptr<Seeds>>> vSoCs;
            for( auto pSeeds : pSoCs->xContent )
            {
                nucSeqIndex uiMin = std::numeric_limits<nucSeqIndex>::max( );
                nucSeqIndex uiMax = 0;
                for( auto& xSeed : *pSeeds )
                {
                    uiMin = std::min( uiMin, xSeed.start( ) );
                    uiMax = std::max( uiMax, xSeed.end( ) );
                } // for
                vSoCs.emplace_back( uiMin, uiMax, pSeeds, pSeeds );
            } // for

            core( vSoCs, removeSeedsInRange, removeSeeds, pOutExtra, valueInRangeSeeds, pQuerySeq, pRefSeq );

            // combine output
            auto pRet = std::make_shared<Seeds>( );

            // for( size_t uiI = 0; uiI < 5 && !pSoCQueue->empty( ); uiI++ )
            //{
            //    std::cout << "score: " << std::get<2>( pSoCQueue->peek_info( ) ) << std::endl;
            //    pRet->append( pSoCQueue->pop( ) );
            //} // for

            // for( auto pSeeds : pSoCs->xContent )
            //    pRet->append( pSeeds );
            for( auto xTuple : vSoCs )
                pRet->append( std::get<2>( xTuple ) );
            if( pOutExtra != nullptr )
                pOutExtra->computeOverlappingSeedsVector( pRet );
            return pRet;
        } // if
        else
        {
            //                     query start, query end  , seed, seed (copy for value calc)
            std::vector<std::tuple<nucSeqIndex, nucSeqIndex, Seed, Seed>> vSeeds;
            for( auto& xSeed : *pSeeds )
                vSeeds.emplace_back( xSeed.start( ), xSeed.end( ), xSeed, xSeed );

            core( vSeeds, removeSeedInRange, removeSeed, pOutExtra, valueInRangeSeed, pQuerySeq, pRefSeq );

            // combine output
            auto pRet = std::make_shared<Seeds>( );
            for( auto& xSeed : vSeeds )
                if( std::get<2>( xSeed ).size( ) > 0 )
                    pRet->push_back( std::get<2>( xSeed ) );
            if( pOutExtra != nullptr )
                pOutExtra->computeOverlappingSeedsVector( pRet );
            return pRet;
        } // else
    } // method

    virtual std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuerySeq, std::shared_ptr<Pack> pRefSeq )
    {
        return execute_helper( pSeeds, pQuerySeq, pRefSeq, nullptr );
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
        if( pParlindromes == nullptr )
            pParlindromes = std::make_shared<Seeds>( );
        else
            pParlindromes->clear( );
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
                xIntervals[ &*xForwStartIt ] = true; // insert
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
