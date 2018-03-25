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
        /// @brief the realtive padding before the first and after the last seed in global alignment
        double fRelativePadding = 1.1;

        //the match missmatch matrix
        const parasail_matrix_t *matrix;

        NeedlemanWunsch(bool bLocal);

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

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
            DEBUG_PARAM(bool bPrintMatrix)
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
}//namespace libMA

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportNeedlemanWunsch();

#endif