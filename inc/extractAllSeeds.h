/** 
 * @file extractAllSeeds.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef EXTRACT_ALL_SEEDS_H
#define EXTRACT_ALL_SEEDS_H

#include "intervalTree.h"
#include "cppModule.h"

/**
 * @brief Used to quickly find areas with high density of @ref Seed "seeds".
 * @ingroup module
 */
class ExtractAllSeeds: public CppModule
{


public:
	unsigned int maxNum = 10;

	ExtractAllSeeds(){}//constructor

	std::shared_ptr<Container> execute(ContainerVector vpInput);

    ContainerVector getInputType() const;

	std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "ExtractAllSeeds";
    }
};//class

/**
 * @brief export the ExtractAllSeeds @ref CppModule "module" to python.
 * @ingroup export
 */
void exportExtractAllSeeds();

#endif