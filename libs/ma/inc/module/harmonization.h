/**
 * @file harmonization.h
 * @brief Implements the harmonization @ref libMA::Module "module"
 * @author Markus Schmidt
 */
#ifndef LINESWEEP_H
#define LINESWEEP_H

#include "IntervalTree.h"
#include "container/segment.h"
#include "container/soc.h"
#include "module/module.h"
#include <algorithm>
#include <csignal>
#include <set>
#include <unordered_map>

//#define PRINT_BREAK_CRITERIA(x) x
#define PRINT_BREAK_CRITERIA( x )

#define USE_RANSAC ( 1 )
#define PI 3.14159265

namespace libMA
{
/**
 * @brief Implements the Harmonization algorithm.
 * @ingroup module
 * @details
 * Removes all contradicting seeds.
 * This should only be used in combination with the StripOfConsideration module.
 */
class HarmonizationSingle : public Module<Seeds, false, Seeds, NucSeq, FMIndex>
{
  private:
    /**
     * @brief The shadow of a Seed.
     * @details
     * Each perfect match "casts a shadow" at the left and right border of the strip.
     * Each shadow is stored in one of these data structures.
     */
    class ShadowInterval : public Interval<int64_t>
    {
      public:
        Seeds::iterator pSeed;

        /**
         * @brief Creates a new shadow.
         * @details
         * The linesweep algorithm disables seeds.
         * Therefore the iterator is required in order to delete the respective seed from its list.
         */
        ShadowInterval( int64_t iBegin, int64_t iSize, Seeds::iterator pSeed )
            : Interval( iBegin, iSize ), pSeed( pSeed )
        {} // constructor

        /**
         * @brief Copy constructor
         */
        ShadowInterval( const ShadowInterval& rOther ) : Interval( rOther ), pSeed( rOther.pSeed )
        {} // copy constructor

        bool within( const ShadowInterval& rOther )
        {
            return start( ) >= rOther.start( ) && end( ) <= rOther.end( );
        } // function
    }; // class


    /**
     * @brief Implements the linesweep algorithm.
     * @details
     * The algorithm has to be run on left and right shadows,
     * therefore it is provided as individual function.
     */
    std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>> EXPORTED
    linesweep( std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>> pShadows,
               const int64_t uiRStart, const double fAngle );

    inline double deltaDistance( const Seed& rSeed, const double fAngle, const int64_t uiRStart )
    {
        double y = rSeed.start_ref( ) + rSeed.start( ) / std::tan( PI / 2 - fAngle );
        double x = ( y - uiRStart ) * std::sin( fAngle );
        double x_1 = rSeed.start( ) / std::sin( PI / 2 - fAngle );

        return std::abs( x - x_1 );
    } // method

    std::pair<double, double> ransac( std::shared_ptr<libMA::Seeds> pSeedsIn );

    std::shared_ptr<libMA::Seeds> applyLinesweeps( std::shared_ptr<libMA::Seeds> pSeedsIn
#if DEBUG_LEVEL > 0
                                                   ,
                                                   bool bRecord
#endif
    );

    std::shared_ptr<Seeds> applyFilters( std::shared_ptr<Seeds>& pIn ) const;

    std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>> cluster( std::shared_ptr<Seeds> pSeedsIn,
                                                                      nucSeqIndex uiQLen ) const;

  public:
    /**
     * @brief true: estimate all possible position as matches in the gap cost filter
     * @details
     * After the linesweep picks a subset of all seeds so that the overall score
     * becomes optimal.
     * This filter has a linear time complexity.
     * It estimates the penalty for gaps between seeds.
     * We have two options here:
     * True -> the gapcost is estimated optimistically (as small as possible)
     * False -> we assume that the score for matches/missmatches is roughly equal within the gap
     *
     * @note not beeing optimistic here has negative affect on the accuracy
     * but improves runtime significantly
     */
    const bool optimisticGapEstimation;

    /// @brief If the seeds cover less that x percent of the query we use SW,
    /// otherwise we fill in the gaps.
    // const double fMinimalQueryCoverage;

    /// @brief Stop the SoC extraction if the harmonized score drops more than x.
    const double fScoreTolerace;

    /// @brief Extract at least x SoCs.
    const size_t uiMinTries;

    /// @brief Lookahead for the equality break criteria.
    const size_t uiMaxEqualScoreLookahead;

    /// @brief Consider two scores equal if they do not differ by more than x (relative to the total
    /// score).
    const double fScoreDiffTolerance;

    /// @brief switch between the two break criteria based on weather the query len is larger than
    /// x.
    const nucSeqIndex uiSwitchQLen;

    /// @brief minimum accumulated seed length after the harmonization as absolute value.
    const nucSeqIndex uiCurrHarmScoreMin;
    /// @brief minimum accumulated seed length after the harmonization relative to query length.
    const double fCurrHarmScoreMinRel;

    const bool bDoHeuristics;

    const bool bDoGapCostEstimationCutting;

    const double dMaxDeltaDist;

    const nucSeqIndex uiMinDeltaDist;

    // const double dMaxSVRatio;

    // const int64_t iMinSVDistance;

    // const nucSeqIndex uiMaxGapArea;

    const size_t uiSVPenalty;

    const nucSeqIndex uiMaxDeltaDistanceInCLuster;

    HarmonizationSingle( const ParameterSetManager& rParameters )
        : optimisticGapEstimation( rParameters.getSelected( )->xOptimisticGapCostEstimation->get( ) ),
          // fMinimalQueryCoverage( rParameters.getSelected( )->xMinQueryCoverage->get( ) ),
          fScoreTolerace( rParameters.getSelected( )->xSoCScoreDecreaseTolerance->get( ) ),
          uiMinTries( rParameters.getSelected( )->xMinNumSoC->get( ) ),
          uiMaxEqualScoreLookahead( rParameters.getSelected( )->xMaxScoreLookahead->get( ) ),
          fScoreDiffTolerance( rParameters.getSelected( )->xScoreDiffTolerance->get( ) ),
          uiSwitchQLen( rParameters.getSelected( )->xSwitchQlen->get( ) ),
          uiCurrHarmScoreMin( rParameters.getSelected( )->xHarmScoreMin->get( ) ),
          fCurrHarmScoreMinRel( rParameters.getSelected( )->xHarmScoreMinRel->get( ) ),
          bDoHeuristics( !rParameters.getSelected( )->xDisableHeuristics->get( ) ),
          bDoGapCostEstimationCutting( !rParameters.getSelected( )->xDisableGapCostEstimationCutting->get( ) ),
          dMaxDeltaDist( rParameters.getSelected( )->xMaxDeltaDist->get( ) ),
          uiMinDeltaDist( rParameters.getSelected( )->xMinDeltaDist->get( ) ),
          // dMaxSVRatio( rParameters.getSelected( )->xMaxSVRatio->get( ) ),
          // iMinSVDistance( rParameters.getSelected( )->xMinSVDistance->get( ) ),
          // uiMaxGapArea( rParameters.getSelected( )->xMaxGapArea->get( ) ),
          uiSVPenalty( rParameters.getSelected( )->xSVPenalty->get( ) ),
          uiMaxDeltaDistanceInCLuster( rParameters.getSelected( )->xMaxDeltaDistanceInCLuster->get( ) )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED
    execute( std::shared_ptr<Seeds> pPrimaryStrand, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<FMIndex> pFMIndex );

}; // class

/**
 * @brief Extracts SoC from a priority queue and performs harmonization on all extracted SoC;s.
 * @ingroup module
 */
class Harmonization : public Module<ContainerVector<std::shared_ptr<Seeds>>, false, SoCPriorityQueue, NucSeq, FMIndex>
{
  public:
    HarmonizationSingle xSingle;

    /// @brief Extract at most x SoCs.
    const size_t uiMaxTries;

    Harmonization( const ParameterSetManager& rParameters )
        : xSingle( rParameters ), uiMaxTries( rParameters.getSelected( )->xMaxNumSoC->get( ) )
    {} // default constructor

    // overload
    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>>
        EXPORTED execute( std::shared_ptr<SoCPriorityQueue> pSoCSsIn, std::shared_ptr<NucSeq> pQuery,
                          std::shared_ptr<FMIndex> pFMIndex )
    {
        unsigned int uiNumTries = 0;

        auto pSoCs = std::make_shared<ContainerVector<std::shared_ptr<Seeds>>>( );

        for( ; uiNumTries < uiMaxTries && !pSoCSsIn->empty( ); uiNumTries++ )
        {
            auto pSoC = pSoCSsIn->pop( );
            // for(auto seed : *pSoC)
            //     std::cout << seed.start() << ", " << seed.start_ref() << ", " << seed.size() << std::endl;
            DEBUG( pSoC->pSoCIn = pSoCSsIn; ) // DEBUG
            // split seeds into forward and reverse strand
            auto pSecondaryStrand = pSoC->extractStrand( false );
            // deal with seeds on forward strand
            while( !pSoC->empty( ) )
                pSoCs->push_back( xSingle.execute( pSoC, pQuery, pFMIndex ) );
            // deal with seeds on reverse strand
            while( !pSecondaryStrand->empty( ) )
                pSoCs->push_back( xSingle.execute( pSecondaryStrand, pQuery, pFMIndex ) );
        } // for

        PRINT_BREAK_CRITERIA( if( pSoCSsIn->empty( ) ) std::cout << "exhausted all SoCs" << std::endl;
                              else std::cout << "break after " << uiNumTries << " tries." << std::endl;
                              std::cout << "computed " << pSoCs->size( ) << " SoCs." << std::endl; )

        return pSoCs;
    } // method
}; // class


/**
 * @brief Combines overlapping seeds
 * @ingroup module
 * @details
 * Uses a n log(n) algorithm to combine overlapping seeds
 */
class SeedLumping : public Module<Seeds, false, Seeds>
{
  public:
    const bool bReduceToMaxSpan = false;

    SeedLumping( const ParameterSetManager& rParameters )
    {} // default constructor

#if 1
    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pIn )
    {
        if( pIn->size( ) == 0 )
            return std::make_shared<Seeds>( );
        /*
         * currently seeds on the reverse strand are pointing to the top left instead of the top right.
         * The following algorithm requires all seeds to point to the top right
         * therefore we mirror the reverse strand seeds to the reverse complement strand (and back once were done)
         *
         * We do not need to know the actual size of the reference in order to do that
         * we can just use the maximal possible reference size.
         */
        const nucSeqIndex uiMaxRefSize = std::numeric_limits<nucSeqIndex>::max( ) / 2 - 10;
        for( Seed& rSeed : *pIn )
            if( !rSeed.bOnForwStrand )
                rSeed.uiPosOnReference = uiMaxRefSize - rSeed.uiPosOnReference;

        std::sort( //
            pIn->begin( ),
            pIn->end( ),
            []( const Seed& rA, const Seed& rB ) //
            { //
                if( rA.bOnForwStrand != rB.bOnForwStrand )
                    return rA.bOnForwStrand;
                if( rA.start_ref( ) - (int64_t)rA.start( ) != rB.start_ref( ) - (int64_t)rB.start( ) )
                    return rA.start_ref( ) - (int64_t)rA.start( ) < rB.start_ref( ) - (int64_t)rB.start( );
                if( rA.start( ) == rB.start( ) )
                    return rA.size( ) < rB.size( );
                return rA.start( ) < rB.start( );
            } // lambda
        ); // sort call

        auto pRet = std::make_shared<Seeds>( );
        pRet->reserve( pIn->size( ) );

        int64_t iDelta = pIn->front( ).start_ref( ) - (int64_t)pIn->front( ).start( );
        pRet->push_back( pIn->front( ) );

        auto xIt = pIn->begin( ) + 1;
        while( xIt != pIn->end( ) )
        {
            int64_t iNewDelta = xIt->start_ref( ) - (int64_t)xIt->start( );
            if( iDelta == iNewDelta && xIt->start( ) <= pRet->back( ).end( ) )
            {
                if( xIt->end( ) > pRet->back( ).end( ) )
                    pRet->back( ).iSize = xIt->end( ) - pRet->back( ).start( );
            } // if
            else
            {
                pRet->push_back( *xIt );
                iDelta = iNewDelta;
            } // else
            xIt++;
        } // while

        // see note on top (mirror result seeds back)
        for( Seed& rSeed : *pRet )
            if( !rSeed.bOnForwStrand )
                rSeed.uiPosOnReference = uiMaxRefSize - rSeed.uiPosOnReference;
        // see note on top (mirror input seeds back; we don't want to destroy the input maybe someone will keep using
        // it)
        for( Seed& rSeed : *pIn )
            if( !rSeed.bOnForwStrand )
                rSeed.uiPosOnReference = uiMaxRefSize - rSeed.uiPosOnReference;

        if( bReduceToMaxSpan )
        {
            // only keep seeds that are not overlapping with longer seeds on query
            std::sort( //
                pRet->begin( ),
                pRet->end( ),
                []( const Seed& rA, const Seed& rB ) //
                { //
                    if( rA.start( ) == rB.start( ) )
                        return rA.end( ) > rB.end( );
                    return rA.start( ) < rB.start( );
                } // lambda
            ); // sort call
            auto pRet2 = std::make_shared<Seeds>( );
            pRet2->reserve( pRet->size( ) );
            for( Seed& rSeed : *pRet )
                if( pRet2->size( ) == 0 || pRet2->back( ).end( ) <= rSeed.end( ) )
                    pRet2->push_back( rSeed );
            return pRet2;
        } // if
        return pRet;
    } // method
#else
    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pIn )
    {
        if( pIn->size( ) == 0 )
            return std::make_shared<Seeds>( );
        /*
         * currently seeds on the reverse strand are pointing to the top left instead of the top right.
         * The following algorithm requires all seeds to point to the top right
         * therefore we mirror the reverse strand seeds to the reverse complement strand (and back once were done)
         *
         * We do not need to know the actual size of the reference in order to do that
         * we can just use the maximal possible reference size.
         */
        const nucSeqIndex uiMaxRefSize = std::numeric_limits<nucSeqIndex>::max( ) / 2 - 10;
        for( Seed& rSeed : *pIn )
            if( !rSeed.bOnForwStrand )
                rSeed.uiPosOnReference = uiMaxRefSize - rSeed.uiPosOnReference;
        std::unordered_map<int64_t, std::vector<Seed>> xHashTable;
        for( Seed& rSeed : *pIn )
            xHashTable[ rSeed.start_ref( ) - (int64_t)rSeed.start( ) ].push_back( rSeed );

        auto pRet = std::make_shared<Seeds>( );
        pRet->reserve( pIn->size( ) );
        for( std::pair<const int64_t, std::vector<Seed>>& rPair : xHashTable )
        {
            std::vector<Seed>& rSeeds = rPair.second;
            std::sort( rSeeds.begin( ), //
                       rSeeds.end( ),
                       []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );
            for( Seed& rSeed : rSeeds )
            {
                if( rSeed.start( ) <= pRet->back( ).end( ) )
                {
                    if( rSeed.end( ) > pRet->back( ).end( ) )
                        pRet->back( ).iSize = rSeed.end( ) - pRet->back( ).start( );
                    assert( pRet->back( ).end( ) == rSeed.end( ) );
                    assert( pRet->back( ).end_ref( ) == rSeed.end_ref( ) );
                } // if
                else
                    pRet->push_back( rSeed );
            } // for
        } // for

        // see note on top (mirror result seeds back)
        for( Seed& rSeed : *pRet )
            if( !rSeed.bOnForwStrand )
                rSeed.uiPosOnReference = uiMaxRefSize - rSeed.uiPosOnReference;
        // see note on top (mirror input seeds back; we don't want to destroy the input maybe someone will keep using
        // it)
        for( Seed& rSeed : *pIn )
            if( !rSeed.bOnForwStrand )
                rSeed.uiPosOnReference = uiMaxRefSize - rSeed.uiPosOnReference;

        if( bReduceToMaxSpan )
        {
            // only keep seeds that are not overlapping with longer seeds on query
            std::sort( //
                pRet->begin( ),
                pRet->end( ),
                []( const Seed& rA, const Seed& rB ) //
                { //
                    if( rA.start( ) == rB.start( ) )
                        return rA.end( ) > rB.end( );
                    return rA.start( ) < rB.start( );
                } // lambda
            ); // sort call
            auto pRet2 = std::make_shared<Seeds>( );
            pRet2->reserve( pRet->size( ) );
            for( Seed& rSeed : *pRet )
                if( pRet2->size( ) == 0 || pRet2->back( ).end( ) <= rSeed.end( ) )
                    pRet2->push_back( rSeed );
            return pRet2;
        } // if
        return pRet;
    } // method
#endif

    virtual std::vector<std::shared_ptr<libMA::Seeds>> lump( std::vector<std::shared_ptr<libMA::Seeds>> vIn )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        for( auto pSeeds : vIn )
            vRet.push_back( execute( pSeeds ) );
        return vRet;
    } // method
}; // class


#if 1
/**
 * @brief Extends seeds at the end to create maximally extened seeds
 * @ingroup module
 */
class SeedExtender : public Module<Seeds, false, Seeds, NucSeq, Pack>
{

  public:
    SeedExtender( const ParameterSetManager& rParameters )
    {} // default constructor
    SeedExtender( )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                                     std::shared_ptr<Pack> pRef )
    {
        for( auto& rSeed : *pSeeds )
        {
            size_t uiForw = 1;
            if( rSeed.bOnForwStrand )
                while( uiForw <= rSeed.start( ) && //
                       uiForw <= rSeed.start_ref( ) && //
                       pQuery->pxSequenceRef[ rSeed.start( ) - uiForw ] ==
                           pRef->getNucleotideOnPos( rSeed.start_ref( ) - uiForw ) )
                    uiForw++;
            else
                while( uiForw <= rSeed.start( ) && //
                       uiForw + rSeed.start_ref( ) < pRef->uiUnpackedSizeForwardStrand && //
                       pQuery->pxSequenceRef[ rSeed.start( ) - uiForw ] ==
                           3 - pRef->getNucleotideOnPos( rSeed.start_ref( ) + uiForw ) )
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
                while( uiBackw + rSeed.end( ) < pQuery->length( ) && //
                       uiBackw + rSeed.end_ref( ) < pRef->uiUnpackedSizeForwardStrand && //
                       pQuery->pxSequenceRef[ rSeed.end( ) + uiBackw ] ==
                           pRef->getNucleotideOnPos( rSeed.end_ref( ) + uiBackw ) )
                    uiBackw++;
            else
                while( uiBackw + rSeed.end( ) < pQuery->length( ) && //
                       rSeed.start_ref( ) >= uiBackw + rSeed.size( ) && //
                       pQuery->pxSequenceRef[ rSeed.end( ) + uiBackw ] ==
                           3 - pRef->getNucleotideOnPos( rSeed.start_ref( ) - rSeed.size( ) - uiBackw ) )
                    uiBackw++;
            rSeed.iSize += uiBackw;
        } // for
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
 * @brief Filters a set of seeds until there are only unique seeds left
 * @details
 * allows up to n missmatches
 * implemented very inefficiently at the moment.
 * @ingroup module
 */
class FilterToUnique : public Module<Seeds, false, Seeds, NucSeq, NucSeq>
{
  public:
    size_t uiNumMissmatchesAllowed = 3;

    FilterToUnique( const ParameterSetManager& rParameters )
    {} // default constructor
    FilterToUnique( )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                                     std::shared_ptr<NucSeq> pRef )
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
class MaxExtendedToSMEM : public Module<Seeds, false, Seeds>
{
  public:
    MaxExtendedToSMEM( const ParameterSetManager& rParameters )
    {} // constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds )
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
class MinLength : public Module<Seeds, false, Seeds>
{
    size_t uiMinLen;

  public:
    MinLength( const ParameterSetManager& rParameters, size_t uiMinLen ) : uiMinLen( uiMinLen )
    {} // constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds )
    {
        pSeeds->erase( std::remove_if( pSeeds->begin( ), pSeeds->end( ),
                                       [&]( const Seed& rSeed ) { return rSeed.size( ) < uiMinLen; } ),
                       pSeeds->end( ) );
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
class MaxExtendedToMaxSpanning : public Module<Seeds, false, Seeds>
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

  public:
    MaxExtendedToMaxSpanning( const ParameterSetManager& rParameters )
    {} // default constructor

    // overload
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds )
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
 * @brief Filters a set of seeds removing all parlindrome seeds
 * @details
 * parlindrome seeds are two seeds that are crossing and on opposite strands.
 * @ingroup module
 */
class ParlindromeFilter : public Module<Seeds, false, Seeds>
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
            // comparing the pointers since we want to std::find to get the exact seed searched for.
            if( pA->start_ref( ) == pB->start_ref( ) )
            {
                if( pA->start( ) == pB->start( ) )
                    return pA->size( ) < pB->size( );
                return pA->start( ) < pB->start( );
            } // if
            return pA->start_ref( ) < pB->start_ref( );
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
    virtual std::shared_ptr<Seeds> EXPORTED execute( std::shared_ptr<Seeds> pSeeds )
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
#if DEBUG_CODE
                std::cout << "AAA " << ( xForwStartIt->bOnForwStrand ? "f " : "r " ) << xForwStartIt->start( ) << ","
                          << xForwStartIt->start_ref( ) << "," << xForwStartIt->size( ) << std::endl;
#endif
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

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Harmonization @ref libMA::Module "module" to boost python.
 * @ingroup export
 */
void exportHarmonization( py::module& rxPyModuleId );
#endif

#endif