/** 
 * @file alignment.h
 * @brief Implements a Container that holds a finished alignment.
 * @author Markus Schmidt
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "intervalTree.h"

/**
 * @brief Holds a finished alignment.
 * @details
 * Contains a list of MatchTypes (match, missmatch, insertion, deletion).
 * @ingroup container
 */
class Alignment : public Container
{
public:
    /**
     * @brief Describes the type of match at one specific position of the alignment.
     * @details
     * @li @c match: query and reference have the same nucleotide.
     * @li @c seed: query and reference have the same nucleotide 
     * (the match was found as part of a seed).
     * @li @c missmatch: query and reference have different nucleotides, 
     *      but they are aligned to the same position nonetheless.
     * @li @c insertion: a nucleotide is present on the query that has no counterpart on the reference.
     * @li @c deletion: a nucleotide is present on the reference that has no counterpart on the query.
     */
    enum MatchType{
        seed,
        match,
        missmatch,
        insertion,
        deletion
    };//enum
private:
    /// The sparse list of MatchTypes that describe the alignment.
    std::vector<std::tuple<MatchType, nucSeqIndex>> data;
    /// The length of the alignment.
    nucSeqIndex uiLength;
    /// Where the actual alignment starts within data.
    nucSeqIndex uiDataStart;
    /// The start of the alignment on the reference sequence.
    nucSeqIndex uiBeginOnRef;
    /// The end of the alignment on the reference sequence.
    nucSeqIndex uiEndOnRef;

public:
    /**
     * @brief Creates an empty alignment.
     */
    Alignment()
            :
        data(),
        uiLength(0),
        uiDataStart(0),
        uiBeginOnRef(0),
        uiEndOnRef(0)
    {}//constructor
    /**
     * @brief Creates an empty aligment, 
     * where the interval of the reference that is used is already known.
     */
    Alignment(nucSeqIndex uiBeginOnRef, nucSeqIndex uiEndOnRef)
            :
        data(),
        uiLength(0),
        uiDataStart(0),
        uiBeginOnRef(uiBeginOnRef),
        uiEndOnRef(uiEndOnRef)
    {}//constructor

    //overload
    ContainerType getType(){return ContainerType::alignment;}

    /**
     * @returns the type of math for the given position i.
     * @brief Type of math at i.
     */
    MatchType at(nucSeqIndex i)
    {
        //everything after the query is a deletion
        if(i >= uiLength)
            return MatchType::deletion;

        //the MatchType match type is stored in a compressed format -> extract it
        nucSeqIndex j = 0;
        unsigned int k = uiDataStart;
        while(k < data.size() && (j += std::get<1>(data[k++])) <= i);

        return std::get<0>(data[k-1]);
    }//function

    /**
     * @returns the type of math for the given position i.
     * @brief Type of math at i.
     */
    MatchType operator[](nucSeqIndex i)
    {
        return at(i);
    }//operator

    /**
     * @brief appends multiple matchTypes to the aligment
     * @details
     * This is used for appending seeds,
     * where simply size of the seed matches need to be appended.
     */
    void append(MatchType type, nucSeqIndex size)
    {
        if(data.size() != 0 && std::get<0>(data.back()) == type)
            std::get<1>(data.back()) += size;
        else
            data.push_back(std::make_tuple(type, size));
        uiLength += size;
    }//function

    /**
     * @brief appends a matchType to the aligment
     */
    void append(MatchType type)
    {
        append(type, 1);
    }//function

    ///@brief wrapper for boost-python
    void append_boost1(MatchType type, nucSeqIndex size)
    {
        append(type, size);
    }//function

    ///@brief wrapper for boost-python
    void append_boost2(MatchType type)
    {
        append(type, 1);
    }//function

    ///@returns the length of the alignment
    ///@brief Length of the alignment
    nucSeqIndex length()
    {
        return uiLength;
    }//function

    ///@brief Start of the alignment.
    ///@returns the start of the alignment on the reference sequence.
    nucSeqIndex beginOnRef()
    {
        return uiBeginOnRef;
    }//function

    ///@brief End of the alignment.
    ///@returns the end of the alignment on the reference sequence.
    nucSeqIndex endOnRef()
    {
        return uiEndOnRef;
    }//function

    /**
     * @brief Remove dangeling deletions.
     * @details
     * Removes parts of the reference at the front and back that overhang the aligned query.
     */
    void removeDangelingDeletions()
    {
        if(data.size() <= 2)
            return;
        if(std::get<0>(data[uiDataStart]) == MatchType::deletion)
        {
            uiBeginOnRef += std::get<1>(data[uiDataStart]);
            uiLength -= std::get<1>(data[uiDataStart]);
            uiDataStart++;
        }//if
        if(std::get<0>(data.back()) == MatchType::deletion)
        {
            uiEndOnRef -= std::get<1>(data.back());
            uiLength -= std::get<1>(data.back());
            data.pop_back();
        }//if
    }//function
};//class

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportAlignment();

#endif