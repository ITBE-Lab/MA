#ifndef MAPPING_QUALITY_H
#define MAPPING_QUALITY_H

#include "module/module.h"
#include "container/alignment.h"

namespace libMABS
{

    /**
     * @brief Sets the mapping quality on alignment
     * @ingroup module
     * @details
     * Given a vector of alignments this module computes the mapping quality for the
     * first alignment on the basis of the second
     */
    class MappingQuality: public Module
    {
    public:

        MappingQuality()
        {}//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - NucSeq
         * - ContainerVector(Alignment)
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Nil
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "MappingQuality";
        }
    };//class
}//namspace libMABS

/**
 * @brief export the MappingQuality @ref Module "module" to python.
 * @ingroup export
 */
void exportMappingQuality();


#endif
