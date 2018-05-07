/** 
 * @file alignment.h
 * @brief Implements a Container that holds a finished alignment.
 * @author Markus Schmidt
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "container/segment.h"


namespace libMA
{
    /**
     * @brief Describes the type of match at one specific position of the alignment.
     * @details
     * @li @c match: query and reference have the same nucleotide.
     * @li @c seed: query and reference have the same nucleotide 
     * (the match was found as part  of a seed).
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
    public:
#if DEBUG_LEVEL >= 1
            std::vector<std::pair<nucSeqIndex, nucSeqIndex>> vGapsScatter;
#endif
        /// The sparse list of MatchTypes that describe the alignment.
        std::vector<std::tuple<MatchType, nucSeqIndex>> data;
        /// The length of the alignment.
        nucSeqIndex uiLength;
        /// The start of the alignment on the reference sequence.
        nucSeqIndex uiBeginOnRef;
        /// The end of the alignment on the reference sequence.
        nucSeqIndex uiEndOnRef;
        /// The start of the alignment on the query sequence.
        nucSeqIndex uiBeginOnQuery;
        /// The end of the alignment on the query sequence.
        nucSeqIndex uiEndOnQuery;

        int iScore = 0;
        double fMappingQuality = NAN;

        //some statistics
        AlignmentStatistics xStats;
        bool bSecondary;

        /**
         * @brief Creates an empty alignment.
         */
        Alignment()
                :
            data(),
            uiLength(0),
            uiBeginOnRef(0),
            uiEndOnRef(0),
            uiBeginOnQuery(0),
            uiEndOnQuery(0),
            xStats(),
            bSecondary(false)
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
            uiBeginOnRef(uiBeginOnRef),
            uiEndOnRef(uiBeginOnRef),
            uiBeginOnQuery(0),
            uiEndOnQuery(0),
            xStats(),
            bSecondary(false)
        {}//constructor

        /**
         * @brief Creates an empty alignment, 
         * where the interval of the reference that is used is already known.
         */
        Alignment(
                nucSeqIndex uiBeginOnRef,
                nucSeqIndex uiBeginOnQuery
            )
                :
            data(),
            uiLength(0),
            uiBeginOnRef(uiBeginOnRef),
            uiEndOnRef(uiBeginOnRef),
            uiBeginOnQuery(uiBeginOnQuery),
            uiEndOnQuery(uiBeginOnQuery),
            xStats(),
            bSecondary(false)
        {}//constructor

        /**
         * @brief Creates an empty alignment, 
         * where the interval of the reference that is used is already known.
         */
        Alignment(
                nucSeqIndex uiBeginOnRef,
                nucSeqIndex uiBeginOnQuery,
                nucSeqIndex uiEndOnRef,
                nucSeqIndex uiEndOnQuery
            )
                :
            data(),
            uiLength(0),
            uiBeginOnRef(uiBeginOnRef),
            uiEndOnRef(uiEndOnRef),
            uiBeginOnQuery(uiBeginOnQuery),
            uiEndOnQuery(uiEndOnQuery),
            xStats(),
            bSecondary(false)
        {}//constructor

        inline std::string toString() const
        {
            std::string sRet = "Alignment Dump: ";
            const char vTranslate[5] = {'S', '=', 'X', 'I', 'D'};
            for(auto& tuple : data)
                sRet += vTranslate[ (unsigned int)std::get<0>(tuple)] + std::to_string(std::get<1>(tuple)) + " ";
            sRet += std::to_string(iScore);
            return sRet; 
        }

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

        int EXPORTED reCalcScore() const;

        /**
         * @returns the type of math for the given position i.
         * @brief Type of math at i.
         */
        MatchType at(nucSeqIndex i) const
        {
            //everything after and before the query is a deletion
            if(i >= uiLength || i < 0)
                return MatchType::deletion;

            //the MatchType match type is stored in a compressed format -> extract it
            nucSeqIndex j = 0;
            unsigned int k = 0;
            while(k < data.size())
            {
                j += std::get<1>(data[k]);
                if(j > i)
                    break;
                k++;
            }// while

            return std::get<0>(data[k]);
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
         * @brief extract the alignment as vector
         */
        std::vector<MatchType> extract() const
        {
            std::vector<MatchType> aRet;
            for(std::tuple<MatchType, nucSeqIndex> xElement : data)
                for(unsigned int i = 0; i < std::get<1>(xElement); i++)
                    aRet.push_back(std::get<0>(xElement));
            return aRet;
        }//function

        /**
         * @brief appends multiple matchTypes to the alignment
         * @details
         * This is used for appending seeds,
         * where simply size of the seed matches need to be appended.
         */
        void EXPORTED append(MatchType type, nucSeqIndex size);

        /**
         * @brief appends a matchType to the alignment
         */
        void append(MatchType type)
        {
            append(type, 1);
        }//function

        /**
         * @brief appends another alignment
         */
        void append(const Alignment& rOther)
        {
            for(auto xTuple : rOther.data)
                append(std::get<0>(xTuple), std::get<1>(xTuple));
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
            DEBUG(
                nucSeqIndex uiCheck = 0;
                for(auto xTup : data)
                    uiCheck += std::get<1>(xTup);
                if(uiCheck != uiLength)
                {
                    std::cout << "Alignment length check failed: " << uiCheck << " != " << uiLength 
                        << std::endl;
                    assert(false);
                }// if
            )//DEBUG
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
            // the data.size() == 0 is to allow SW to set the score directly
            assert(data.size() == 0 || reCalcScore() == iScore);
            return iScore;
        }

        /**
         * @brief the NMW score for this alignment
         */
        unsigned int EXPORTED localscore() const;

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
         * @brief for sorting alignment by their score
         * @details
         * When multiple alignments are created we use this function to sort them.
         * Overload from Module.
         */
        bool larger(const std::shared_ptr<Container> pOther) const
        {
            const std::shared_ptr<Alignment>& pAlign = std::dynamic_pointer_cast<Alignment>(pOther);
            if(pAlign == nullptr)
                return false;
            auto uiS1 = score();
            auto uiS2 = pAlign->score();
            if(uiS1 == uiS2)
                // if both alignments have the same score output the one with 
                // the higher SoC score first (this is determined by the lower SoC index)
                return xStats.index_of_strip < pAlign->xStats.index_of_strip;
            return uiS1 > uiS2;
        }//function

        /**
         * @brief transform any alignment into a local one
         * @details
         * When an alignment is computed on the foundation of seeds it might not be local.
         * This function has a linear complexity with regard to the compressed alignment length.
         */
        void EXPORTED makeLocal();

        /**
         * @brief removes dangeling Deletions
         * @details
         * When the alignment is created there might be some dangeling deletions at 
         * the beginning or end. This function removes them
         */
        void EXPORTED removeDangeling();

        void operator=(const std::shared_ptr<Alignment> pOther)
        {
            data = pOther->data;
            uiLength = pOther->uiLength;
            uiBeginOnRef = pOther->uiBeginOnRef;
            uiEndOnRef = pOther->uiEndOnRef;
            uiBeginOnQuery = pOther->uiBeginOnQuery;
            iScore = pOther->iScore;
            fMappingQuality = pOther->fMappingQuality;
            bSecondary = pOther->bSecondary;
            xStats = pOther->xStats;
        }//function

        inline void shiftOnRef(nucSeqIndex uiBy)
        {
            uiBeginOnRef += uiBy;
            uiEndOnRef += uiBy;
        }// method

        inline void shiftOnQuery(nucSeqIndex uiBy)
        {
            uiBeginOnQuery += uiBy;
            uiEndOnQuery += uiBy;
        }// method
    };//class
}//namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportAlignment();
#endif

#endif