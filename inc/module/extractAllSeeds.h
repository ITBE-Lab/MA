/**
 * @file extractAllSeeds.h
 * @brief Implements ExtractAllSeeds.
 * @author Markus Schmidt
 */
#ifndef EXTRACT_ALL_SEEDS_H
#define EXTRACT_ALL_SEEDS_H

#include "container/segment.h"
#include "module/module.h"

namespace libMA
{
/**
 * @brief Extracts all Seeds from a SegmentList.
 * @ingroup module
 */
class ExtractAllSeeds : public Module
{
   public:
    /**
     * @brief the maximal ambiguity a extracted seed shall have.
     * @details
     * Will skip any segments that would lead to more than maxAmbiguity seeds.
     */
    unsigned int maxAmbiguity = defaults::uiMaxAmbiguity;
    /**
     * @brief the minimal seed length required
     * @details
     * Will skip any segments that are shorter than uiMinLen.
     */
    unsigned int uiMinLen = defaults::uiMinLen;

    ExtractAllSeeds( )
    {
    } // default constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - SegmentVector
     * - FMIndex
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - Seeds
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "ExtractAllSeeds";
    }

    std::string getFullDesc( ) const
    {
        return std::string( "ExtractAllSeeds(" ) + std::to_string( maxAmbiguity ) + ")";
    } // function
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the ExtractAllSeeds @ref Module "module" to python.
 * @ingroup export
 */
void exportExtractAllSeeds( );
#endif

#endif