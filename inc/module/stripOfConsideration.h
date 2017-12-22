/** 
 * @file stripOfConsideration.h
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "container/segment.h"
#include "module/module.h"

namespace libMABS
{
    /**
     * @brief Used to quickly find areas with high density of @ref Seed "seeds".
     * @ingroup module
     */
    class StripOfConsideration: public Module
    {
    public:
        /// @brief The strip of consideration size.
        nucSeqIndex uiStripSize = 10000;
        /// @brief Maximum ambiguity for a seed to be considered.
        unsigned int uiMaxAmbiguity = 500;
        /**
        * @brief skip seeds with too much ambiguity
        * @details
        * True: skip all seeds with to much ambiguity
        * False: use max_hits instances of the seeds with more ambiguity
        */
        bool bSkipLongBWTIntervals = true;
        
    private:
        inline nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS) const 
        { 
            return xS.start_ref() + (uiQueryLength - xS.start()); 
        }//function

        void forEachNonBridgingSeed(
                std::shared_ptr<SegmentVector> pVector,
                std::shared_ptr<FMIndex> pxFM_index,std::shared_ptr<Pack> pxRefSequence,
                std::shared_ptr<NucSeq> pxQuerySeq,
                std::function<void(Seed)> fDo,
                nucSeqIndex addSize// = 0 (default)
            );

    public:

        StripOfConsideration(){}//default constructor

        StripOfConsideration(nucSeqIndex uiStripSize, unsigned int uiMaxAmbiguity)
                :
            uiStripSize(uiStripSize),
            uiMaxAmbiguity(uiMaxAmbiguity)
        {}//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - SegmentVector
         * - Seeds
         * - NucSeq
         * - Pack
         * - FMIndex
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Seeds)
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "StripOfConsideration";
        }
    };//class
}//namspace libMABS

/**
 * @brief export the bucketing @ref Module "module" to python.
 * @ingroup export
 */
void exportStripOfConsideration();

#endif