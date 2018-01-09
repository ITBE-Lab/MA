/** 
 * @file alignment.h
 * @brief Implements a Container that holds a finished alignment.
 * @author Markus Schmidt
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "container/segment.h"
#include "container/seed.h"


namespace libMABS
{
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

    /**
     * @brief Holds a finished alignment.
     * @details
     * Contains a list of MatchTypes (match, missmatch, insertion, deletion).
     * @ingroup container
     */
    class Alignment : public Container
    {
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
        /// The start of the alignment on the query sequence.
        nucSeqIndex uiBeginOnQuery;
        /// The end of the alignment on the query sequence.
        nucSeqIndex uiEndOnQuery;

        int iScore = 0;

    public:
        //some statistics
        AlignmentStatistics xStats;

        /**
         * @brief Creates an empty alignment.
         */
        Alignment()
                :
            data(),
            uiLength(0),
            uiDataStart(0),
            uiBeginOnRef(0),
            uiEndOnRef(0),
            uiBeginOnQuery(0),
            uiEndOnQuery(0),
            xStats()
        {}//constructor
        /**
         * @brief Creates an empty alignment, 
         * where the interval of the reference that is used is already known.
         */
        Alignment(
                nucSeqIndex uiBeginOnRef, 
                nucSeqIndex uiEndOnRef,
                nucSeqIndex uiBeginOnQuery,
                nucSeqIndex uiEndOnQuery
            )
                :
            data(),
            uiLength(0),
            uiDataStart(0),
            uiBeginOnRef(uiBeginOnRef),
            uiEndOnRef(uiEndOnRef),
            uiBeginOnQuery(uiBeginOnQuery),
            uiEndOnQuery(uiEndOnQuery),
            xStats()
        {}//constructor

        /**
         * @brief Creates an empty alignment, 
         * where the interval of the reference that is used is already known.
         */
        Alignment(
                nucSeqIndex uiBeginOnRef
            )
                :
            data(),
            uiLength(0),
            uiDataStart(0),
            uiBeginOnRef(uiBeginOnRef),
            uiEndOnRef(0),
            xStats()
        {}//constructor

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<Alignment>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "Alignment";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new Alignment());
        }//function

        /**
         * @returns the type of math for the given position i.
         * @brief Type of math at i.
         */
        MatchType at(nucSeqIndex i) const
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
        MatchType operator[](nucSeqIndex i) const
        {
            return at(i);
        }//operator

        /**
         * @brief appends multiple matchTypes to the alignment
         * @details
         * This is used for appending seeds,
         * where simply size of the seed matches need to be appended.
         */
        void append(MatchType type, nucSeqIndex size);

        /**
         * @brief appends a matchType to the alignment
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
        nucSeqIndex length() const
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
         * @brief the NMW score for this alignment
         */
        int score() const
        {
            return iScore;
        }

        /**
         * @brief returns how many nucleotides within this alignment are determined by seeds 
         * as a percentage
         */
        float seedCoverage() const
        {
            unsigned int iCount = 0;
            for(unsigned int i = 0; i < length(); i++)
                if(at(i) == MatchType::seed)
                    iCount++;
            return ((float)iCount)/ (float)length();
        }

        /*
        * @brief returns how many nucleotides within this alignment are determined by seeds
        */
        unsigned int numBySeeds() const
        {
            unsigned int iCount = 0;
            for(unsigned int i = 0; i < length(); i++)
                if(at(i) == MatchType::seed)
                    iCount++;
            return iCount;
        }

        /**
         * @brief Remove dangeling deletions.
         * @details
         * Removes parts of the reference at the front and back that overhang the aligned query.
         */
        void removeDangelingDeletions();

        /**
         * @brief for sorting alignment by their score
         * @details
         * When multiple alignments are created we use this function to sort them 
         * overload from Module
         */
        bool smaller(const std::shared_ptr<Container> pOther) const
        {
            const std::shared_ptr<Alignment> pAlign = std::dynamic_pointer_cast<Alignment>(pOther);
            if(pAlign == nullptr)
                return false;
            return score() < pAlign->score();
        }//function
    };//class
}//namespace libMABS

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportAlignment();

#endif