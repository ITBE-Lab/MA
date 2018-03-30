/** 
 * @file linesweep.h
 * @brief Implements the linesweep @ref Module "module"
 * @author Markus Schmidt
 */
#ifndef LINESWEEP
#define LINESWEEP

#include "container/segment.h"
#include "module/module.h"

namespace libMA
{
    /**
     * @brief Implements the LinearLineSweep algorithm.
     * @ingroup module
     * @details
     * Removes all contradicting seeds.
     * This should only be used in combination with the StripOfConsideration module.
     */
    class LinearLineSweep: public Module
    {
    private:
        /**
         * @brief The shadow of a Seed.
         * @details
         * Each perfect match "casts a shadow" at the left and right border of the strip.
         * Each shadow is stored in one of these data structures.
         */
        class ShadowInterval: public Interval<int64_t>{
        public:
            Seeds::iterator pSeed;

            /**
             * @brief Creates a new shadow.
             * @details
             * The linesweep algorithm disables seeds.
             * Therefore the iterator is required in order to delete the respective seed from its list.
             */
            ShadowInterval(
                    int64_t iBegin, 
                    int64_t iSize, 
                    Seeds::iterator pSeed
                )
                    :
                Interval(iBegin, iSize),
                pSeed(pSeed)
            {}//constructor

            /**
             * @brief Copy constructor
             */
            ShadowInterval( const ShadowInterval& rOther )
                    :
                Interval(rOther),
                pSeed(rOther.pSeed)
            {}//copy constructor

            bool within(const ShadowInterval& rOther)
            {
                return start() >= rOther.start() && end() <= rOther.end();
            }//function
        };//class


        /**
        * @brief Implements the linesweep algorithm.
        * @details
        * The algorithm has to be run on left and right shadows,
        * therefore it is provided as individual function.
        */
        std::shared_ptr<std::vector<std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>>>
        EXPORTED linesweep(
            std::shared_ptr<std::vector<
                std::tuple<Seeds::iterator, nucSeqIndex, nucSeqIndex>
            >> pShadows
        );

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
         * @note not beeing optimistic here has no negative affect on the accuracy
         * but improves runtime significantly
         */
        bool optimisticGapEstimation = false;

        LinearLineSweep() {}//default constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> pInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Seeds
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Seeds
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "LineSweep2";
        }

        std::string getFullDesc() const
        {
            return std::string("LineSweep2(") + 
                std::to_string(optimisticGapEstimation) + ")"
                ;
        }//function
    };//class
}//namespace libMA

/**
 * @brief Exposes the LineSweep @ref Module "module" to boost python.
 * @ingroup export
 */
void exportLinesweep();

#endif