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

            /**
             * @brief Removes this or all interfering seeds.
             * @details
             * Checks weather this seed is more valuable than the interfering ones.
             */
            void remove(
                    std::shared_ptr<Seeds> pSeeds
                )
            {
                //reached an already invalidated iterator
                if(pSeed == pSeeds->end())
                {
                    DEBUG(
                        std::cout << "\treached invalid iterator" << std::endl;
                    )
                    return;
                }
                pSeeds->erase(pSeed);
                pSeed = pSeeds->end();
            }//function

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
        void EXPORTED linesweep(
                std::vector<ShadowInterval>& vShadows, 
                std::shared_ptr<Seeds> pSeeds
            );

        /**
        * @brief Returns the left shadow of a seed.
        * @details
        * "Casts" the left shadows.
        */
        ShadowInterval EXPORTED getLeftShadow(Seeds::iterator pSeed) const;

        /**
        * @brief Returns the right shadow of a seed.
        * @details
        * "Casts" the right shadows.
        */
        ShadowInterval EXPORTED getRightShadow(Seeds::iterator pSeed) const;
    public:
        bool optimisticGapEstimation = false;

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
    };//class
}//namespace libMA

/**
 * @brief Exposes the LineSweep @ref Module "module" to boost python.
 * @ingroup export
 */
void exportLinesweep();

#endif