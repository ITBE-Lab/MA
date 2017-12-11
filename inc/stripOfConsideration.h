/** 
 * @file bucketing.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "segmentList.h"
#include "cppModule.h"

/**
 * @brief Used to quickly find areas with high density of @ref Seed "seeds".
 * @ingroup module
 */
class StripOfConsideration: public CppModule
{
private:

public:
	/// @brief The strip of consideration size.
	nucSeqIndex uiStripSize = 10000;
	/// @brief Maximum ambiguity for a seed to be considered.
	unsigned int uiMaxHitsPerInterval = 500;
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
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::function<void(Seed)> fDo,
			nucSeqIndex addSize// = 0 (default)
		);

public:

	StripOfConsideration(){}//constructor

	std::shared_ptr<Container> execute(ContainerVector vpInput);

    ContainerVector getInputType() const;

	std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "StripOfConsideration";
    }
};//class

/**
 * @brief export the bucketing @ref CppModule "module" to python.
 * @ingroup export
 */
void exportStripOfConsideration();

#endif