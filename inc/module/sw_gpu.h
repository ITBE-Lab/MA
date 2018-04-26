/**
 * @file sw_gpu.h
 * @author Markus Schmidt
 * @brief execute a given module on each element of a given libMA::ContainerVector.
 */

#ifndef SW_GPU_H
#define SW_GPU_H

#include <vector>
#include <cstddef>

class GPUReturn
{
public:
    int iMaxScore;
    std::vector<size_t> vMaxPos;
    GPUReturn(int iMaxScore, std::vector<size_t> vMaxPos)
            :
        iMaxScore(iMaxScore),
        vMaxPos(vMaxPos)
    {}// default constructor
    GPUReturn(){}

    bool operator==(const GPUReturn& rOther)
    {
        return iMaxScore == rOther.iMaxScore;
    }// operator
};// class

std::vector<GPUReturn> cudaAlign
(
    std::vector<char> &rvRefSeq, // reference sequence
	std::vector<std::vector<char>> &rvQuerySeqs // vector of query sequences
);

/**
 * @brief Exposes the SweepAllReturnBest @ref Module "module" to boost python.
 * @ingroup export
 */
void exportSW_GPU();

#endif
