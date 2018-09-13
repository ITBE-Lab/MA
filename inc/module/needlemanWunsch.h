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
#include "gaba.h"

#define ALLOCATE_ONCE ( 0 ) // naive approach is disabled
#define NAIVE_MAX_SIZE 5

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
        /**
         * @brief wrapper that take care of deallocation for the Gaba library
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

        std::vector<std::vector<std::vector<int>>> s;
        std::vector<std::vector<std::vector<char>>> dir;

        void naiveNeedlemanWunsch(
                std::shared_ptr<NucSeq> pQuery, 
                std::shared_ptr<NucSeq> pRef,
                nucSeqIndex fromQuery,
                nucSeqIndex toQuery,
                nucSeqIndex fromRef,
                nucSeqIndex toRef,
                bool bNoGapAtBeginning,
                bool bNoGapAtEnd,
                std::shared_ptr<Alignment> pAlignment
            );

        void EXPORTED dynPrg(
            const std::shared_ptr<NucSeq> pQuery, 
            const std::shared_ptr<NucSeq> pRef,
            const nucSeqIndex fromQuery, const nucSeqIndex toQuery,
            const nucSeqIndex fromRef, const nucSeqIndex toRef,
            std::shared_ptr<Alignment> pAlignment, // in & output
            const bool bLocalBeginning,
            const bool bLocalEnd
            );

        /*
        * @todo: fix this
        * 
        * arghh this is really ugly...
        * At the moment there is only one sequence that sets all NW and SW parameters correctly:
        * 1) set the parameters.
        * 2) create a new NW module.
        * 3) Then create all other modules that you want to use.....
        */
        //the match missmatch matrix
        std::shared_ptr<Gaba_tWrapper> pGabaScoring;

        nucSeqIndex uiMaxGapArea = defaults::uiMaxGapArea;
    public:
        bool bLocal = false;

        NeedlemanWunsch();

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
            return "NeedlemanWunsch";
        }// method

        std::string getFullDesc() const;
    };//class

}//namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportNeedlemanWunsch();
#endif

#endif