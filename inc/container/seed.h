/** 
 * @file seed.h
 * @brief Implements Seed.
 * @author Markus Schmidt
 */
#ifndef SEED_H
#define SEED_H

#include "container/container.h"
#include "container/interval.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <list>
/// @endcond

#define DELTA_CACHE ( 1 )
#define CONTIG_ID_CACHE ( 0 )

namespace libMA
{
    ///@brief any index on the query or reference nucleotide sequence is given in this datatype
    typedef uint64_t nucSeqIndex;

    /**
     * @brief A seed.
     * @details
     * A extracted seed, that comprises two intervals, one on the query one on the reference.
     * Both intervals are equal in size.
     * @note the overloaded functions of Interval refer to the Interval on the query.
     * @ingroup container
     */
    class Seed: public Container, public Interval<nucSeqIndex>
    {
    public:
        ///@brief the beginning of the match on the reference
        nucSeqIndex uiPosOnReference;
        unsigned int uiAmbiguity;
#if DELTA_CACHE == ( 1 )
        nucSeqIndex uiDelta = 0;
#endif
#if CONTIG_ID_CACHE == ( 1 )
        size_t uiContigId = 0;
#endif

        /**
         * @brief Creates a new Seed.
         */
        Seed(
                const nucSeqIndex uiPosOnQuery, 
                const nucSeqIndex uiLength, 
                const nucSeqIndex uiPosOnReference
            )
                :
            Interval(uiPosOnQuery, uiLength),
            uiPosOnReference(uiPosOnReference),
            uiAmbiguity(0)
        {}//constructor

        /**
         * @brief Creates a new Seed.
         */
        Seed(
                const nucSeqIndex uiPosOnQuery, 
                const nucSeqIndex uiLength, 
                const nucSeqIndex uiPosOnReference,
                const unsigned int uiAmbiguity
            )
                :
            Interval(uiPosOnQuery, uiLength),
            uiPosOnReference(uiPosOnReference),
            uiAmbiguity(uiAmbiguity)
        {}//constructor

        /**
         * @brief Copys from a Seed.
         */
        Seed(const Seed& rOther)
                :
            Interval(rOther),
            uiPosOnReference(rOther.uiPosOnReference),
            uiAmbiguity(rOther.uiAmbiguity)
#if DELTA_CACHE == ( 1 )
            , uiDelta(rOther.uiDelta)
#endif
#if CONTIG_ID_CACHE == ( 1 )
            , uiContigId(rOther.uiContigId)
#endif
        {}//copy constructor

        /**
         * @brief Default Constructor.
         */
        Seed()
                :
            Interval()
        {}//default constructor
        
        /**
         * @brief Returns the beginning of the seed on the reference.
         */
        inline nucSeqIndex start_ref() const
        {
            return uiPosOnReference;
        }//function
        
        /**
         * @brief Returns the end of the seed on the reference.
         */
        inline nucSeqIndex end_ref() const
        {
            return uiPosOnReference + size();
        }//function
        
        /**
         * @brief Returns the value of the seed.
         * @details
         * A seeds value corresponds to its size.
         */
        inline nucSeqIndex getValue() const
        {
            return size();
        }//function

        /**
         * @brief Copys from another Seed.
         */
        inline Seed& operator=(const Seed& rxOther)
        {
            Interval::operator=(rxOther);
            uiPosOnReference = rxOther.uiPosOnReference;
            uiAmbiguity = rxOther.uiAmbiguity;
#if DELTA_CACHE == ( 1 )
            uiDelta = rxOther.uiDelta;
#endif
#if CONTIG_ID_CACHE == ( 1 )
            uiContigId = rxOther.uiContigId;
#endif
            return *this;
        }// operator
        
        /*
        * @brief compares two Seeds.
        * @returns true if start and size are equal, false otherwise.
        */
        inline bool operator==(const Seed& rxOther)
        {
            return Interval::operator==(rxOther) && uiPosOnReference == rxOther.uiPosOnReference;
        }// operator

        //overload
        inline bool canCast(const std::shared_ptr<Container>& c) const
        {
            return std::dynamic_pointer_cast<Seed>(c) != nullptr;
        }//function

        //overload
        inline std::string getTypeName() const
        {
            return "Seed";
        }//function

        //overload
        inline std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new Seed());
        }//function

    }; //class

    class Alignment;
    /**
     * @brief Used to store some statistics to each alignment
     * @details
     * Intended for figuring out optimal thresholds.
     */
    class AlignmentStatistics
    {
    public:
        unsigned int index_of_strip;
        unsigned int num_seeds_in_strip;
        unsigned int anchor_size;
        unsigned int anchor_ambiguity;
        std::weak_ptr<Alignment> pOther;
        bool bFirst;
        std::string sName;
        nucSeqIndex uiInitialQueryBegin;
        nucSeqIndex uiInitialRefBegin;
        nucSeqIndex uiInitialQueryEnd;
        nucSeqIndex uiInitialRefEnd;

        AlignmentStatistics()
                :
            index_of_strip(0),
            num_seeds_in_strip(0),
            anchor_size(0),
            anchor_ambiguity(0),
            pOther(),
            bFirst(false),
            sName("unknown"),
            uiInitialQueryBegin(0),
            uiInitialRefBegin(0),
            uiInitialQueryEnd(0),
            uiInitialRefEnd(0)
        {}

        void operator=(const AlignmentStatistics &rOther)
        {
            index_of_strip = rOther.index_of_strip;
            num_seeds_in_strip = rOther.num_seeds_in_strip;
            anchor_size = rOther.anchor_size;
            anchor_ambiguity = rOther.anchor_ambiguity;
            pOther = rOther.pOther;
            bFirst = rOther.bFirst;
            sName = rOther.sName;
            uiInitialQueryBegin = rOther.uiInitialQueryBegin;
            uiInitialRefBegin = rOther.uiInitialRefBegin;
            uiInitialQueryEnd = rOther.uiInitialQueryEnd;
            uiInitialRefEnd = rOther.uiInitialRefEnd;
        }//function
    };//class

    class SVInfo
    {
    public:
        /**
         * @details
         * Contains indices of seeds, where a structural variant (deletion or insertion)
         * follows the seed.
         */
        std::vector<size_t> vSeedIndicesOfSVIndels;
    };// class

    /**
     * @brief A list where one element is a Seed.
     * @details
     * Also holds the summed up score of the seeds within the list.
     * @ingroup Container
     */
    class Seeds
        :
            public std::vector<Seed>,// @todo :@
            public Container
    {
    public:
        nucSeqIndex mem_score = 0;
        //some statistics
        AlignmentStatistics xStats;
        //SVInfo xSvInfo;
        bool bConsistent;

        Seeds(std::shared_ptr<Seeds> pOther)
                :
            vector(),
            Container(),
            xStats(pOther->xStats),
            bConsistent(pOther->bConsistent)
        {
            append(pOther);
        }//copy constructor

        Seeds()
                :
            vector(),
            Container(),
            xStats(),
            bConsistent(false)
        {}//default constructor

        template< class InputIt >
        Seeds(InputIt xBegin, InputIt xEnd)
                :
            vector(xBegin, xEnd)
        {}//iterator constructor

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<Seeds>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "Seeds";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new Seeds());
        }//function

        /*returns the sum off all scores within the list*/
        nucSeqIndex EXPORTED getScore() const;

        /*append a copy of another list*/
        void append(std::shared_ptr<Seeds> pOther)
        {
            for(Seed& rS : *pOther)
                push_back(rS);
        }//function

        bool larger(const std::shared_ptr<Container> pOther) const
        {
            const std::shared_ptr<Seeds> pSeeds = std::dynamic_pointer_cast<Seeds>(pOther);
            if(pSeeds == nullptr)
                return true;
            return getScore() > pSeeds->getScore();
        }// operator

        // inline bool hasSV() const
        // {
        //     return ! xSvInfo.vSeedIndicesOfSVIndels.empty();
        // }// method
        // 
        // inline std::shared_ptr<Seeds> partitionOnSV()
        // {
        //     assert(hasSV());
        // 
        //     auto pRet = std::make_shared<Seeds>();
        //     for(size_t uiPos = xSvInfo.vSeedIndicesOfSVIndels.back(); uiPos<this->size();uiPos++)
        //         if( (*this)[uiPos].size() != 0 )
        //             pRet->emplace_back( (*this)[uiPos] );
        //     this->resize(xSvInfo.vSeedIndicesOfSVIndels.back());
        //     xSvInfo.vSeedIndicesOfSVIndels.pop_back();
        //     assert( ! this->empty());
        //     assert( ! pRet->empty());
        //     return pRet;
        // }// method

    };//class
}//namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief exports the Seed and Seedlist classes to python.
 * @ingroup export
 */
void exportSeed();
#endif

#endif