/**
 * @file mappingQuality.h
 * @brief Computes the mapping quality of Alignments.
 * @author Markus Schmidt
 */
#ifndef MAPPING_QUALITY_H
#define MAPPING_QUALITY_H

#include "container/alignment.h"
#include "module/module.h"

namespace libMA
{
/**
 * @brief Sets the mapping quality on alignment
 * @ingroup module
 * @details
 * Given a vector of alignments this module computes the mapping quality for the
 * first alignment on the basis of the second
 * @note the name quality is missleading it rather is a mapping confidence
 */
class MappingQuality
    : public Module<ContainerVector<std::shared_ptr<Alignment>>, false, NucSeq, ContainerVector<std::shared_ptr<Alignment>>>
{
  public:
    size_t uiReportNBest = defaults::uiReportN;
    size_t uiMinAlignmentScore = defaults::uiMinAlignmentScore;
    double dMaxOverlapSupplementary = defaults::dMaxOverlapSupplementary;
    size_t uiMaxSupplementaryPerPrim = defaults::uiMaxSupplementaryPerPrim;

    MappingQuality( )
    {} // constructor

    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
        EXPORTED execute( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments );

}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the MappingQuality @ref Module "module" to python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportMappingQuality( );
#else
void exportMappingQuality( py::module& rxPyModuleId );
#endif
#endif


#endif
