/** 
 * @file bucketing.h
 * @brief Implements a Bucketing.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "intervalTree.h"
#include "cppModule.h"

/**
 * @brief a bucket for the seeds
 * @ingroup container
 */
 class SeedBucket
 {
 private:
	nucSeqIndex uiTotalScore;
 
	std::vector<Seed> lxContent;
 
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
		SYNC(std::lock_guard<std::mutex> xGuard(xMutex);)
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

	const std::vector<Seed>& seeds() const
	{
		return lxContent;
	}//function
 };//class

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

	void addSeed(
			nucSeqIndex uiQueryLength, 
			const Seed xNew, 
			std::vector<SeedBucket>& raxSeedBuckets
		)
	{
		assert(getPositionForBucketing(uiQueryLength, xNew) / uiStripSize < raxSeedBuckets.size());
		raxSeedBuckets[getPositionForBucketing(uiQueryLength, xNew) / uiStripSize].addSeed(xNew);
	}//function

	void forEachNonBridgingSeed(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pxFM_index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::function<void(Seed)> fDo,
			nucSeqIndex addSize// = 0 (default)
		);
	
	void saveSeeds(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pxFM_index,
			std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::vector<SeedBucket>& raxSeedBuckets
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