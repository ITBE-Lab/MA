/** 
 * @file bucketing.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef BUCKETING_H
#define BUCKETING_H

#include "graphicalMethod.h"
#include "module.h"

class Bucketing: public Module
{
private:

public:
	unsigned int uiNumThreads = 8;
	nucSeqIndex uiStripSize = 1000;
	unsigned int uiMaxHitsPerInterval = 1000;
	bool bSkipLongBWTIntervals = true;
	
private:
	nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength, const Seed xS) const 
	{ 
		return xS.start_ref() + (uiQueryLength - xS.start()); 
	}//function

	void addSeed(
			nucSeqIndex uiQueryLength, 
			const Seed xNew, 
			std::vector<SeedBucket>& raxSeedBuckets
		)
	{
		raxSeedBuckets[getPositionForBucketing(uiQueryLength, xNew) / uiStripSize].addSeed(xNew);
	}//function

	void forEachNonBridgingSeed(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			bool bAnchorOnly,
			std::shared_ptr<FM_Index> pxFM_index,
			std::shared_ptr<FM_Index> pxRev_FM_Index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::function<void(Seed)> fDo
		);
	
	void saveSeeds(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pxFM_index,
			std::shared_ptr<FM_Index> pxRev_FM_Index,
			std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::vector<SeedBucket>& raxSeedBuckets
		);

public:

	Bucketing(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

	std::vector<ContainerType> getOutputType();
};//class

void exportBucketing();

#endif