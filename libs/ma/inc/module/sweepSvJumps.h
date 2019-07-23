/**
 * @file sweepSvJumps.h
 * @author Markus Schmidt
 */
#pragma once

#include "container/squeezedVector.h"
#include "module/module.h"

namespace libMA
{

class CompleteBipartiteSubgraphCluster: public Container
{

}; // class

/**
 * @brief
 * @details
 */
class SweepSvJumps : public Module<CompleteBipartiteSubgraphCluster, true, TP_CONTAINER>
{
  public:
    /**
     * @brief
     * @details
     */
    SweepSvJumps( const ParameterSetManager& rParameters )
    {} // constructor

    virtual std::shared_ptr<TP_CONTAINER> EXPORTED execute( std::shared_ptr<TP_CONTAINER> pInput )
    {
        return pInput;
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportSweepSvJump( py::module& rxPyModuleId );
#endif