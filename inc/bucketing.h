/** 
 * @file bucketing.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef BUCKETING_H
#define BUCKETING_H

#include "intervalTree.h"
#include "module.h"

/**
 * @brief a bucket for the seeds
 * @note TODO: remove me!
 */
 class SeedBucket
 {
 private:
	 nucSeqIndex uiTotalScore;
 
	 std::list<Seed> lxContent;
 
	 std::mutex xMutex;
 
 public:
	SeedBucket()
		:
		uiTotalScore(0),
		lxContent(),
		xMutex()
		//,bUsed(false)
	{}//constructor

	////disable copying of buckets
	SeedBucket(const SeedBucket&) = delete;

	void addSeed(const Seed xNew)
	{
		std::lock_guard<std::mutex> xGuard(xMutex);
		lxContent.push_back(xNew);
		uiTotalScore += xNew.size();
		//end of scope xGuard
	}//function

	nucSeqIndex getValue()
	{
		return uiTotalScore;
	}//function

	void forall(std::function<void(const Seed&)> fDo)
	{
		for (auto xSeed : lxContent)
		{
			fDo(xSeed);
		}//for
	}//function

	const std::list<Seed>& seeds() const
	{
		return lxContent;
	}//function
 };//class

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
		assert(getPositionForBucketing(uiQueryLength, xNew) / uiStripSize <raxSeedBuckets.size());
		raxSeedBuckets[getPositionForBucketing(uiQueryLength, xNew) / uiStripSize].addSeed(xNew);
	}//function

	void forEachNonBridgingSeed(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pxFM_index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::function<void(Seed)> fDo
		);
	
	void saveSeeds(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pxFM_index,
			std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::vector<SeedBucket>& raxSeedBuckets
		);

public:

	Bucketing(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

	ContainerType getOutputType();
};//class

void exportBucketing();

#endif