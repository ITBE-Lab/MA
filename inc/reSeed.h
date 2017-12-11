/** 
 * @file longestNonEnclosedSegments.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */

#ifndef RE_SEED_H
#define RE_SEED_H

//#define DEBUG_ENABLED

#include "system.h"
#include "cppModule.h"
#include "segmentList.h"
#include "threadPool.h"


/**
 * @brief Computes a set of maximal non-enclosed seeds.
 * @ingroup module
 */
class ReSeed : public CppModule{
private:

	static SAInterval extend_backward(
			const SAInterval &ik, 
			const uint8_t c, 
			std::shared_ptr<FMIndex> pFM_index
		);

	/*
	 * bwa style extension
	 */
	static void extend(
			std::shared_ptr<SegmentVector> pxVector,
			nucSeqIndex min,
			nucSeqIndex max,
			std::shared_ptr<FMIndex> pFM_index,
			std::shared_ptr<NucSeq> pQuerySeq
		);


public:
	nucSeqIndex minSplitLen = 16;
	
	std::shared_ptr<Container> execute(ContainerVector vpInput);

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - FMIndex
     * - SegmentVector
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
};//class


/**
 * @brief exports the Segmentation @ref CppModule "module" to python.
 * @ingroup export
 */
void exportReSeed();

#endif