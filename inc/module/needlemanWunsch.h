/** 
 * @file needlemanWunsch.h
 * @brief Implements NMW
 * @author Markus Schmidt
 */
#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include "container/alignment.h"
#include "module/module.h"
// The NW library:
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#include "parasail/matrix_lookup.h"

namespace libMA
{
    /**
     * @brief implements NMW
     * @details
     * Returns a finished alignment if given a sound selection of seeds.
     * @ingroup module
     */
    class NeedlemanWunsch : public Module
    {
    public:
        bool bLocal;
        /// @brief The realtive padding before the first and after the last seed.
        double fRelativePadding = 1.1;
        /// @brief If the seeds cover less that x percent of the query we use SW, 
        /// otherwise we fill in the gaps.
        double fMinimalQueryCoverage = .25;

        //the match missmatch matrix
        parasail_matrix_t matrix;
        std::vector<int> vMatrixContent;

        NeedlemanWunsch(bool bLocal);

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief the NW dynamic programming algorithm
         * @details 
         * Uses parasail to have an efficient vectorized implementation.
         * For the moment: is either bNoGapAtBeginning or bNoGapAtEnd it jumps to 
         * the naive implementation.
         * 
         * @TODO: at the moment ugly C code is used here (free functions).. 
         * find a way to replace that?
         */
        nucSeqIndex needlemanWunsch(
            std::shared_ptr<NucSeq> pQuery, 
            std::shared_ptr<NucSeq> pRef,
            nucSeqIndex fromQuery,
            nucSeqIndex toQuery,
            nucSeqIndex fromRef,
            nucSeqIndex toRef,
            std::shared_ptr<Alignment> pAlignment,
            bool bNoGapAtBeginning,
            bool bNoGapAtEnd
            DEBUG_PARAM(bool bPrintMatrix = false)
        );

        /**
         * @brief the SW dynamic programming algorithm
         * @details 
         * Uses parasail to have an efficient vectorized implementation.
         * We use SW if the query coverage of the seeds is so little
         * that there is no point in filling gaps.
         * 
         * 
         * @TODO: at the moment ugly C code is used here (free functions).. 
         * find a way to replace that?
         */
        std::shared_ptr<Alignment> smithWaterman(
            std::shared_ptr<NucSeq> pQuery, 
            std::shared_ptr<NucSeq> pRef,
            nucSeqIndex uiOffsetRef
        );

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Seeds
         * - NucSeq
         * - Pack
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Alignment
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "NeedlemanWunsch";
        }

        std::string getFullDesc() const;
    };//class

    
    /**
     * @brief Converts a local to an global alignment
     * @details
     * Checks the mapping quality of the best alignment.
     * If the mapping quality is low we convert into a global alignment 
     * hopefully improving our decision.
     * @ingroup module
     */
    class LocalToGlobal : public Module
    {
    public:
        /// @brief The realtive padding before the first and after the last seed.
        const double fMappingQualMin;
        const unsigned int uiReturnBestN;
        /// @brief The realtive padding before the first and after the last seed.
        double fRelativePadding = 1.1;

        LocalToGlobal(const double fMappingQualMin, const unsigned int uiReturnBestN)
                :
            fMappingQualMin(fMappingQualMin),
            uiReturnBestN(uiReturnBestN)
        {}//constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Seeds
         * - NucSeq
         * - Pack
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - Alignment
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "LocalToGlobal";
        }

        std::string getFullDesc() const;
    };//class
}//namespace libMA

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportNeedlemanWunsch();

#endif