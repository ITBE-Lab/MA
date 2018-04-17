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
#include "gaba.h"

namespace libMA
{
    /**
     * @todo this class should not be in the .h
     * @brief wrapper for parsail results.
     * @details
     * This will automatically call free the result in the deconstructor and therefore make everything
     * exception save.
     */
    class ParsailResultWrapper
    {
        parasail_result_t* pContent;
    public:

        ParsailResultWrapper(parasail_result_t* pContent)
                :
            pContent(pContent)
        {}//constructor

        ~ParsailResultWrapper()
        {
            parasail_result_free(pContent);
        }//deconstructor

        const parasail_result_t& operator*() const
        {
            return *pContent;
        }//operator

        const parasail_result_t* get() const
        {
            return pContent;
        }//operator

        parasail_result_t* get()
        {
            return pContent;
        }//operator

        const parasail_result_t* operator->() const
        {
            return pContent;
        }//operator
    };//class

    /**
     * @todo this class should not be in the .h
     * @brief wrapper for parsail cigars.
     * @details
     * This will automatically call free the result in the deconstructor and therefore make everything
     * exception save.
     */
    class ParsailCigarWrapper
    {
        parasail_cigar_t* pContent;
    public:

        ParsailCigarWrapper(parasail_cigar_t* pContent)
                :
            pContent(pContent)
        {}//constructor

        ~ParsailCigarWrapper()
        {
            parasail_cigar_free(pContent);
        }//deconstructor

        const parasail_cigar_t& operator*() const
        {
            return *pContent;
        }//operator

        const parasail_cigar_t* get() const
        {
            return pContent;
        }//operator

        parasail_cigar_t* get()
        {
            return pContent;
        }//operator

        const parasail_cigar_t* operator->() const
        {
            return pContent;
        }//operator
    };//class

    /**
     * @todo this class should not be in the .h
     */
    class Gaba_tWrapper
    {
    public:
        gaba_params_s xParams;
        gaba_t* pContext;

        Gaba_tWrapper(const gaba_params_s &rxParams)
                :
            xParams(rxParams)
        {
            pContext = gaba_init(&this->xParams);
        }//constructor

        Gaba_tWrapper()
                :
            pContext(nullptr)
        {}//constructor

        ~Gaba_tWrapper()
        {
            if(pContext != nullptr)
                gaba_clean(pContext);
        }//deconstructor

        void operator=(const Gaba_tWrapper &rOther) = delete;
        Gaba_tWrapper(const Gaba_tWrapper& rOther) = delete;
    };//class

    class Gaba_dp_tWrapper
    {
    public:
        gaba_dp_t* pDp;
        struct gaba_alignment_s * pR;

        Gaba_dp_tWrapper(gaba_dp_t* pDp)
                :
            pDp(pDp),
            pR(nullptr)
        {}//constructor

        ~Gaba_dp_tWrapper()
        {
            if(pR != nullptr)
                gaba_dp_res_free(pDp, pR);
            gaba_dp_clean(pDp);
        }//deconstructor
    };//class

    //@todo: this should not be in a .h file
    void EXPORTED dynPrg(
        const std::shared_ptr<NucSeq> pQuery, 
        const std::shared_ptr<NucSeq> pRef,
        const nucSeqIndex fromQuery, const nucSeqIndex toQuery,
        const nucSeqIndex fromRef, const nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment, // in & output
        const bool bLocalBeginning,
        const bool bLocalEnd
        );


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
        /// @brief If the seeds cover less that x percent of the query we use SW, 
        /// otherwise we fill in the gaps.
        double fMinimalQueryCoverage = .25;

        NeedlemanWunsch(bool bLocal);

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

#if 0
        /**
         * @brief the NW dynamic programming algorithm
         * @details 
         * Uses parasail to have an efficient vectorized implementation.
         * For the moment: is either bNoGapAtBeginning or bNoGapAtEnd it jumps to 
         * the naive implementation.
         * 
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
#endif

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

        LocalToGlobal(double fMappingQualMin)
                :
            fMappingQualMin(fMappingQualMin)
        {}//constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
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

    /**
     * @brief @todo
     * @ingroup module
     */
    class CombatRepetitively : public Module
    {
    public:
        const double fMappingQualMax;
        const nucSeqIndex uiRegionLength;

        CombatRepetitively(double fMappingQualMax, nucSeqIndex uiRegionLength)
                :
            fMappingQualMax(fMappingQualMax),
            uiRegionLength(uiRegionLength)
        {}//constructor

        //overload
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
         * - NucSeq
         * - Pack
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "CombatRepetitively";
        }

        std::string getFullDesc() const
        {
            return std::string("CombatRepetitively(") + 
                std::to_string(fMappingQualMax) + "," +
                std::to_string(uiRegionLength) + ")"
                ;
        }//function
    };//class
}//namespace libMA

/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportNeedlemanWunsch();

#endif