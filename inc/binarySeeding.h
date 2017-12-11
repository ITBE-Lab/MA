/** 
 * @file binarySeeding.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef BINARY_SEEDING_H
#define BINARY_SEEDING_H

#define DEBUG_ENABLED

#include "system.h"
#include "cppModule.h"
#include "segmentList.h"
#include "threadPool.h"

class PerfectMatch;

/**
 * @brief Computes all non-enclosed seeds.
 * @details
 * This is the the BWA seeding strategy, parallelized using binary seeding.
 * Gives you all non-enclosed seeds.
 * Therefore there are a lot of bad quality seeds and strict filters need to be applied.
 * @ingroup module
 */
class BinarySeeding : public CppModule{
private:
	bool bLrExtension;

	static SAInterval extend_backward(
			const SAInterval &ik, 
			const uint8_t c, 
			std::shared_ptr<FMIndex> pFM_index
		);

	/**
	 * TODO:
	 */
	Interval<nucSeqIndex> lrExtension(
			nucSeqIndex center,
			std::shared_ptr<FMIndex> pFM_index,
			std::shared_ptr<NucSeq> pQuerySeq,
			std::shared_ptr<SegmentVector> pSegmentVector
		);

	/**
	 * TODO:
	 */
	Interval<nucSeqIndex> nonEnclosedExtension(
			nucSeqIndex center,
			std::shared_ptr<FMIndex> pFM_index,
			std::shared_ptr<NucSeq> pQuerySeq,
			std::shared_ptr<SegmentVector> pSegmentVector
		);

	/*
	*	does nothing if the given interval can be found entirely on the genome.
	*	if the interval cannot be found this method splits the interval in half and repeats the step with the first half,
	*	while queuing the second half as a task in the thread pool.
	*/
	void procesInterval(
			Interval<nucSeqIndex> xAreaToCover,
			std::shared_ptr<SegmentVector> pSegmentVector,
			std::shared_ptr<FMIndex> pFM_index,
			std::shared_ptr<NucSeq> pQuerySeq,
			ThreadPoolAllowingRecursiveEnqueues* pxPool
		);
	
	/*
	 * functions need to be static in order to enqueue them into a threadpool
	 * we need to enqueue procesInterval with an object associated
	 * workaround: give the object as second (since pool will give tId as first) parameter 
	 */
	static void procesIntervalStatic(
			size_t uiThreadId,
			BinarySeeding *obj,
			Interval<nucSeqIndex> xAreaToCover,
			std::shared_ptr<SegmentVector> pSegmentVector,
			std::shared_ptr<FMIndex> pFM_index,
			std::shared_ptr<NucSeq> pQuerySeq,
			ThreadPoolAllowingRecursiveEnqueues* pxPool
		)
	{
		obj->procesInterval(
				xAreaToCover,
				pSegmentVector,
				pFM_index,
				pQuerySeq,
				pxPool
			);
	}//function

public:
	BinarySeeding(bool bLrExtension = true)
			:
		bLrExtension(bLrExtension)
	{}//constructor
	
	std::shared_ptr<Container> execute(ContainerVector vpInput);

	/**
	 * @brief Used to check the input of execute.
	 * @details
	 * Returns:
	 * - FMIndex
	 * - NucSeq
	 */
	ContainerVector getInputType() const;

	/**
	 * @brief Used to check the output of execute.
	 * @details
	 * Returns:
	 * - SegmentVector
	 */
    std::shared_ptr<Container> getOutputType() const;

    std::string getName() const
    {
        return "BinarySeeding";
    }
};//class


/**
 * @brief exports the Segmentation @ref CppModule "module" to python.
 * @ingroup export
 */
void exportBinarySeeding();

#endif