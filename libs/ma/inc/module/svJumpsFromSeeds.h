/**
 * @file svJumpsFromSeeds.h
 * @brief Implements a way to compute SV-jumps from seeds.
 * @author Markus Schmidt
 */
#pragma once

#include "container/segment.h"
#include "container/svJump.h"
#include "module/module.h"

namespace libMA
{
class PerfectMatch;

#define complement( x ) ( uint8_t ) NucSeq::nucleotideComplement( x )


/**
 * @brief Computes Sv-Jumps from a given seed set
 */
class SvJumpsFromSeeds : public Module<ContainerVector<SvJump>, false, SegmentVector, Pack, FMIndex, NucSeq>
{
  public:
    /**
     * @brief Initialize a SvJumpsFromSeeds Module
     */
    SvJumpsFromSeeds( const ParameterSetManager& rParameters )
    {} // constructor

    virtual std::shared_ptr<ContainerVector<SvJump>> EXPORTED execute( std::shared_ptr<SegmentVector> pSegments,
                                                                       std::shared_ptr<Pack> pRefSeq,
                                                                       std::shared_ptr<FMIndex> pFM_index,
                                                                       std::shared_ptr<NucSeq> pQuery );
}; // class

}; // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief exports the SvJumpsFromSeeds @ref Module "module" to python.
 * @ingroup export
 */
void exportSvJumpsFromSeeds( py::module& rxPyModuleId );
#endif