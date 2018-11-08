/**
 * @file soc.h
 * @brief Some Helper classes for the SoC computation
 * @author Markus Schmidt
 */
#ifndef SOC_H
#define SOC_H

#include "container/seed.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <algorithm>
/// @endcond

namespace libMA
{
/**
 * @brief Used to determine more complex score-orders of SoCs
 * @details
 * In this case we determine the SoC scores as follows:
 * We count the accumulative seed length and the accumulative ambiguity.
 * SoCs are sorted according to their accumulative seed length.
 * If two SoCs have the same accumulative seed length we break the tie by sorting the SoC with
 * less accumulative ambiguity first.
 */
class SoCOrder
{
  public:
    /// @brief Stores the accumulative seed length.
    nucSeqIndex uiAccumulativeLength = 0;
    /// @brief Stores the accumulative seed ambiguity.
    unsigned int uiSeedAmbiguity = 0;
    /// @brief Stores the amount of seeds in the SoC. (Currently not used for the Score.)
    unsigned int uiSeedAmount = 0;

    /**
     * @brief Add another seed to the SoC score.
     * @details
     * Adjusts all three stored values accordingly.
     */
    inline void operator+=( const Seed& rS )
    {
        uiSeedAmbiguity += rS.uiAmbiguity;
        uiSeedAmount++;
        uiAccumulativeLength += rS.getValue( );
    } // operator

    /**
     * @brief Remove a seed from the SoC score.
     * @details
     * Adjusts all three stored values accordingly.
     * The user needs to ensure that the seed that is removed
     * was added to this SoC score before.
     * Otherwise this may result in undefined behaviour while sorting SoCs.
     */
    inline void operator-=( const Seed& rS )
    {
        assert( uiSeedAmbiguity >= rS.uiAmbiguity );
        uiSeedAmbiguity -= rS.uiAmbiguity;
        assert( uiAccumulativeLength >= rS.getValue( ) );
        uiAccumulativeLength -= rS.getValue( );
        uiSeedAmount--;
    } // operator

    /**
     * @brief Compare two scores for smaller.
     * @details
     * Compares the accumulative seed length and
     * uses the accumulative ambiguity as a tie breaker.
     */
    inline bool operator<( const SoCOrder& rOther ) const
    {
        if( uiAccumulativeLength == rOther.uiAccumulativeLength )
            return uiSeedAmbiguity > rOther.uiSeedAmbiguity;
        return uiAccumulativeLength < rOther.uiAccumulativeLength;
    } // operator

    /**
     * @brief Compare two scores for equality.
     * @details
     * Compares the accumulative seed length and
     * uses the accumulative ambiguity as a tie breaker.
     */
    inline void operator=( const SoCOrder& rOther )
    {
        uiAccumulativeLength = rOther.uiAccumulativeLength;
        uiSeedAmbiguity = rOther.uiSeedAmbiguity;
    } // operator
}; // class


/**
 * @brief Acts as stack during collection of SoCs then turns into a max-heap for SoC extraction.
 */
class SoCPriorityQueue : public Container
{
  public:
#if DEBUG_LEVEL >= 1
    /// Confirms the respective functions are always called in the correct mode
    bool bInPriorityMode = false;
    std::vector<std::pair<nucSeqIndex, nucSeqIndex>> vScores;
    class blub
    {
      public:
        nucSeqIndex first = 0, second = 0, qCoverage = 0, rStart = 0, rEnd = 0, rStartSoC = 0, rEndSoC = 0;
        inline void operator=( const blub& rOther )
        {
            first = rOther.first;
            second = rOther.second;
            qCoverage = rOther.qCoverage;
            rStart = rOther.rStart;
            rEnd = rOther.rEnd;
            rStartSoC = rOther.rStartSoC;
            rEndSoC = rOther.rEndSoC;
        } // operator
        inline bool operator==( const blub& rOther ) const
        {
            return first == rOther.first && second == rOther.second && qCoverage == rOther.qCoverage &&
                   rStart == rOther.rStart && rEnd == rOther.rEnd && rStartSoC == rOther.rStartSoC &&
                   rEndSoC == rOther.rEndSoC;
        } // operator
    };
    std::vector<blub> vExtractOrder;
    std::vector<std::shared_ptr<Seeds>> vSoCs;
    std::vector<std::shared_ptr<Seeds>> vHarmSoCs;
    std::vector<double> vSlopes;
    std::vector<double> vIntercepts;
    std::vector<std::shared_ptr<Seeds>> vIngroup;
    nucSeqIndex uiQLen;
#endif
    /// @brief The index of the next SoC during the extraction process.
    unsigned int uiSoCIndex = 0;
    /// @brief The complete seed set.
    std::shared_ptr<Seeds> pSeeds;
    /// @brief Contains the SoCs in for of tuples (score, start, end).
    std::vector<std::tuple<SoCOrder, std::vector<Seed>::iterator, std::vector<Seed>::iterator>> vMaxima;
    /// @brief End position of the last SoC during collection; required to determine overlaps.
    nucSeqIndex uiLastEnd = 0;
    /// @brief Function used to for the make_max_heap call.
    static bool heapOrder( const std::tuple<SoCOrder, std::vector<Seed>::iterator, std::vector<Seed>::iterator>& rA,
                           const std::tuple<SoCOrder, std::vector<Seed>::iterator, std::vector<Seed>::iterator>& rB )
    {
        return std::get<0>( rA ) < std::get<0>( rB );
    } // method

    /**
     * @brief Create a new SoC priority queue with pSeeds as the seed set.
     * @details
     * The queue has two states and push_back_no_overlap can only be called in the first state,
     * while pop can only be called in the second state.
     * Initially the queue is in the first state.
     * Call make_heap to switch to the second state.
     * If compiled in debug mode this class checks for correct usage.
     * In release mode incorrect usage results in undefined behaviour.
     */
    SoCPriorityQueue( std::shared_ptr<Seeds> pSeeds ) : pSeeds( pSeeds ), vMaxima( )
    {} // constructor

    /**
     * @brief Create a new SoC priority queue without a seed set.
     */
    SoCPriorityQueue( ) : pSeeds( nullptr ), vMaxima( )
    {} // constructor

    /// @brief Returns the size (usable in either state).
    inline size_t size( ) const
    {
        return vMaxima.size( );
    }

    // overload
    inline bool canCast( const std::shared_ptr<Container>& c ) const
    {
        return std::dynamic_pointer_cast<SoCPriorityQueue>( c ) != nullptr;
    } // method

    // overload
    inline std::string getTypeName( ) const
    {
        return "SoCPriorityQueue";
    } // method

    // overload
    inline std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new SoCPriorityQueue( ) );
    } // method

    /// @brief Returns weather the queue is empty (usable in either state).
    inline bool empty( ) const
    {
        return vMaxima.empty( );
    } // method

    inline void rectangularSoC( )
    {
        std::vector<std::pair<nucSeqIndex, nucSeqIndex>> vRefPosMaxima;
        // populate vRefPosMaxima with the reference start and end positions of all SoCs.
        for( auto& tSoC : vMaxima )
        {
            auto xItStart = std::get<1>( tSoC );
            auto xItEnd = std::get<2>( tSoC );
            vRefPosMaxima.emplace_back( xItStart->start_ref( ), xItStart->start_ref( ) );
            while( xItStart != xItEnd )
            {
                vRefPosMaxima.back( ).first = std::min( vRefPosMaxima.back( ).first, xItStart->start_ref( ) );
                vRefPosMaxima.back( ).second = std::max( vRefPosMaxima.back( ).second, xItStart->start_ref( ) );
                xItStart++;
            } // while
        } // for
        // sort seeds after the start and end positions.
        std::sort( pSeeds->begin( ), pSeeds->end( ),
                   []( const Seed& rA, const Seed& rB ) { return rA.start_ref( ) < rB.start_ref( ); } );
        // repopulate vMaxima
        vMaxima.clear( );
        for( auto& rP : vRefPosMaxima )
        {
            vMaxima.emplace_back( );
            std::get<1>( vMaxima.back( ) ) =
                std::lower_bound( pSeeds->begin( ), pSeeds->end( ), rP.first,
                                  []( const Seed& rS, nucSeqIndex uiPos ) { return rS.start_ref( ) < uiPos; } );
            auto xIt = std::get<1>( vMaxima.back( ) );
            while( xIt != pSeeds->end( ) && xIt->start_ref( ) <= rP.second )
            {
                std::get<0>( vMaxima.back( ) ) += *xIt;
                xIt++;
            } // while
            std::get<2>( vMaxima.back( ) ) = xIt;
        } // for
    } // method

    /**
     * @brief Returns the first SoC (usable in second state only).
     * @details
     * After the queue has been turned into a max heap you may extract the SoCs
     * in order of their scores.
     * Pop removes and returns the best SoC.
     */
    inline std::shared_ptr<Seeds> pop( )
    {
        DEBUG( assert( !empty( ) ); assert( bInPriorityMode ); ) // DEBUG
        // the strip that shall be collected
        std::shared_ptr<Seeds> pRet( new Seeds( ) );
        // get the expected amount of seeds in the SoC from the order class and reserve memory
        pRet->reserve( std::get<0>( vMaxima.front( ) ).uiSeedAmount );
        // iterator walking till the end of the strip that shall be collected
        auto xCollect2 = std::get<1>( vMaxima.front( ) );
        assert( xCollect2 != pSeeds->end( ) );
        // save SoC index
        pRet->xStats = pSeeds->xStats;
        pRet->xStats.index_of_strip = uiSoCIndex++;
        // all these things are not used at the moment...
        pRet->xStats.uiInitialQueryBegin = xCollect2->start( );
        pRet->xStats.uiInitialRefBegin = xCollect2->start_ref( );
        pRet->xStats.uiInitialQueryEnd = xCollect2->end( );
        pRet->xStats.uiInitialRefEnd = xCollect2->end_ref( );
        auto xCollectEnd = std::get<2>( vMaxima.front( ) );
        while( xCollect2 != pSeeds->end( ) && xCollect2 != xCollectEnd )
        {
            // save the beginning and end of the SoC
            // all these things are not used at the moment...
            if( xCollect2->start( ) < pRet->xStats.uiInitialQueryBegin )
                pRet->xStats.uiInitialQueryBegin = xCollect2->start( );
            if( xCollect2->start_ref( ) < pRet->xStats.uiInitialRefBegin )
                pRet->xStats.uiInitialRefBegin = xCollect2->start_ref( );
            if( xCollect2->end( ) > pRet->xStats.uiInitialQueryEnd )
                pRet->xStats.uiInitialQueryEnd = xCollect2->end( );
            if( xCollect2->end_ref( ) > pRet->xStats.uiInitialRefEnd )
                pRet->xStats.uiInitialRefEnd = xCollect2->end_ref( );
            pRet->xStats.num_seeds_in_strip++;
            assert( xCollect2->start( ) <= xCollect2->end( ) );
            // if the iterator is still within the strip add the seed and increment the iterator
            pRet->push_back( *( xCollect2++ ) );
        } // while


        // move to the next strip
        std::pop_heap( vMaxima.begin( ), vMaxima.end( ), heapOrder );
        vMaxima.pop_back( );

        return pRet;
    } // method

    /**
     * @brief Returns information about the first SoC (usable in second state only).
     * @details
     * After the queue has been turned into a max heap you may extract the SoCs
     * in order of their scores.
     * Pop removes and the best SoC and returns information about it's start end and score.
     */
    inline std::tuple<nucSeqIndex, nucSeqIndex, uint32_t, double> pop_info( )
    {
        DEBUG( assert( !empty( ) ); assert( bInPriorityMode ); ) // DEBUG

        // the information that shall be returned
        auto xCollect = std::get<1>( vMaxima.front( ) );
        auto xCollectEnd = std::get<2>( vMaxima.front( ) );
        nucSeqIndex uiStart = xCollect->start_ref( );
        nucSeqIndex uiEnd = (xCollectEnd-1)->end_ref( );
        uint32_t uiScore = std::get<0>( vMaxima.front( ) ).uiAccumulativeLength;
        size_t uiForwSeeds = 0;
        size_t uiSeedsTotal = 0;
        while( xCollect != pSeeds->end( ) && xCollect != xCollectEnd )
        {
            if( xCollect->bOnForwStrand )
                uiForwSeeds++;
            uiSeedsTotal++;
            xCollect++;
        } // while

        // move to the next strip
        std::pop_heap( vMaxima.begin( ), vMaxima.end( ), heapOrder );
        vMaxima.pop_back( );

        return std::make_tuple( uiStart, uiEnd, uiScore, uiForwSeeds / (double)uiSeedsTotal );
    } // method

    /**
     * @brief Add a new SoC (usable in first state only).
     * @details
     * While collecting SoCs this queue functions as a Stack.
     * If the new and the last SoC are overlapping we only want to keep the higher scored one.
     * This functions either replaces the last SoC if the new one has a higher score or
     * discards the new SoC in case of overlaps.
     * Otherwise the new SoC is merely pushed onto the stack.
     * @note itStripEnd points to one element past end of the SoC.
     */
    inline void push_back_no_overlap( const SoCOrder& rCurrScore, const std::vector<Seed>::iterator itStrip,
                                      const std::vector<Seed>::iterator itStripEnd, const nucSeqIndex uiCurrStart,
                                      const nucSeqIndex uiCurrEnd )
    {
        DEBUG( assert( !bInPriorityMode ); )
        if( vMaxima.empty( ) || uiLastEnd < uiCurrStart || std::get<0>( vMaxima.back( ) ) < rCurrScore )
        {
            // if we reach this point we want to save the current SoC
            if( !vMaxima.empty( ) && uiLastEnd >= uiCurrStart )
                // the new and the last SoC overlap and the new one has a higher score
                // so we want to replace the last SoC
                vMaxima.pop_back( );
            // else the new and the last SoC do not overlapp
            // so we want to save the new SOC
            vMaxima.push_back( std::make_tuple( rCurrScore, itStrip, itStripEnd ) );
            // we need to remember the end position of the new SoC
            uiLastEnd = uiCurrEnd;
        } // if
        // else new and last SoC overlap and the new one has a lower score => ignore the new one
    } // method

    /**
     * @brief switch to second state (usable in first state only).
     * @details
     * Use make max heap to turn this stack into a priority queue.
     * Once this function is called push_back_no_overlap results in undefined behaviour
     * however pop can now be used safely.
     */
    inline void make_heap( )
    {
        DEBUG( assert( !bInPriorityMode ); bInPriorityMode = true; )
        // make a max heap from the SOC starting points according to the scores,
        // so that we can extract the best SOC first
        std::make_heap( vMaxima.begin( ), vMaxima.end( ), heapOrder );
    } // method
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportSoC( );
#else
void exportSoC( py::module& rxPyModuleId );
#endif
#endif

#endif