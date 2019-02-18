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
class MappingQuality : public Module<ContainerVector<std::shared_ptr<Alignment>>, false, NucSeq,
                                     ContainerVector<std::shared_ptr<Alignment>>>
{
  public:
    const size_t uiReportNBest;
    const size_t uiMinAlignmentScore;
    const double dMaxOverlapSupplementary;
    const size_t uiMaxSupplementaryPerPrim;

    MappingQuality( const ParameterSetManager& rParameters )
        : uiReportNBest( rParameters.getSelected( )->xReportN.get( ) ),
          uiMinAlignmentScore( rParameters.getSelected( )->xMinAlignmentScore.get( ) ),
          dMaxOverlapSupplementary( rParameters.getSelected( )->xMaxOverlapSupplementary.get( ) ),
          uiMaxSupplementaryPerPrim( rParameters.getSelected( )->xMaxSupplementaryPerPrim.get( ) )
    {} // constructor

    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> EXPORTED
    execute( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments );

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
