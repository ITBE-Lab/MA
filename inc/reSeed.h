/** 
 * @file longestNonEnclosedSegments.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */

#ifndef LONGEST_NON_ENCLOSED_SEGMENTS_H
#define LONGEST_NON_ENCLOSED_SEGMENTS_H

#define DEBUG_ENABLED

#include "system.h"
#include "cppModule.h"
#include "intervalTree.h"
#include "threadPool.h"


/**
 * @brief Computes a set of maximal non-enclosed seeds.
 * @ingroup module
 */
class ReSeed : public CppModule{
private:

	static SA_IndexInterval extend_backward(
			const SA_IndexInterval &ik, 
			const uint8_t c, 
			std::shared_ptr<FM_Index> pFM_index
		);

	/*
	 * bwa style extension
	 */
	static SaSegment extend(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pFM_index,
			std::shared_ptr<NucleotideSequence> pQuerySeq
		);


public:
	nucSeqIndex minSplitLen = 16;
	
	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

	std::vector<std::shared_ptr<Container>> getInputType() const;
    std::shared_ptr<Container> getOutputType() const;
};//class


/**
 * @brief exports the Segmentation @ref CppModule "module" to python.
 * @ingroup export
 */
void exportReSeed();

#endif