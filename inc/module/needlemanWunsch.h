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
    std::vector<std::vector<std::vector<int>>> s;
    std::vector<std::vector<std::vector<char>>> dir;

    void naiveNeedlemanWunsch( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef,
                               nucSeqIndex fromQuery, nucSeqIndex toQuery, nucSeqIndex fromRef,
                               nucSeqIndex toRef, bool bNoGapAtBeginning, bool bNoGapAtEnd,
                               std::shared_ptr<Alignment> pAlignment );

    void EXPORTED dynPrg( const std::shared_ptr<NucSeq> pQuery, const std::shared_ptr<NucSeq> pRef,
                          const nucSeqIndex fromQuery, const nucSeqIndex toQuery,
                          const nucSeqIndex fromRef, const nucSeqIndex toRef,
                          std::shared_ptr<Alignment> pAlignment, // in & output
                          const bool bLocalBeginning, const bool bLocalEnd );

    nucSeqIndex uiMaxGapArea = defaults::uiMaxGapArea;

   public:
    bool bLocal = false;

    NeedlemanWunsch( );

    // overload
    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );


    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Seeds
     * - NucSeq
     * - Pack
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - Alignment
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    std::string getName( ) const
    {
        return "NeedlemanWunsch";
    } // method

    std::string getFullDesc( ) const;
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
void exportNeedlemanWunsch( );
#endif

#endif