/**
 * @file fileWriter.cpp
 * @author Markus Schmidt
 */
#include "module/dbWriter.h"

#ifdef WITH_POSTGRES
using namespace libMA;


std::shared_ptr<Container> DbWriter::execute( std::shared_ptr<NucSeq> pQuery,
                                              std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                  pAlignments,
                                              std::shared_ptr<Pack>
                                                  pPack )
{

    for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
    {
        if( pAlignment->length( ) == 0 )
            continue;
        std::string sCigar = pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sSegment = pAlignment->getQuerySequence( *pQuery, *pPack );
        std::string sTlen = std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery );

        std::string sRefName = pAlignment->getContig( *pPack );
        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition( *pPack );

        double fMappingQual = pAlignment->fMappingQuality;


        std::string sSQL = "";

        sSQL += "INSERT INTO alignment (cigar, position, mapping_quality, query_id, run_id, sam_flags, "
                "contig, query_sequence, contig_other, position_other, length) VALUES ( \'";

        sSQL += sCigar + "\', ";
        sSQL += std::to_string( uiRefPos ) + ", ";
        sSQL += std::to_string( fMappingQual ) + ", \'";
        sSQL += pAlignment->xStats.sName + "\', ";
        sSQL += std::to_string( iRunId ) + ", ";
        sSQL += std::to_string( flag ) + ", \'";
        sSQL += sRefName + "\', \'";
        sSQL += sSegment + "\', \'";
        sSQL += sContigOther + "\', ";
        sSQL += sPosOther + ", ";
        sSQL += sTlen + " )";

        // std::cerr << sSQL << std::endl;

        xConnection.exec( sSQL );
    } // for

    return std::shared_ptr<Container>( new Container( ) ); // return empty container
} // function


std::shared_ptr<Container> PairedDbWriter::execute( std::shared_ptr<NucSeq> pQuery1,
                                                    std::shared_ptr<NucSeq>
                                                        pQuery2,
                                                    std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                        pAlignments,
                                                    std::shared_ptr<Pack>
                                                        pPack )
{

    for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
    {
        if( pAlignment->length( ) == 0 )
            continue;
        std::string sCigar = pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sSegment = "";
        std::string sTlen = std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery );
        // paired
        if( pAlignment->xStats.pOther.lock( ) != nullptr )
        {
            // flag |= pAlignment->xStats.bFirst ? FIRST_IN_TEMPLATE // flag not actually required
            //                                  : LAST_IN_TEMPLATE;
            flag |= MULTIPLE_SEGMENTS_IN_TEMPLATE | SEGMENT_PROPERLY_ALIGNED;
            if( pPack->bPositionIsOnReversStrand( pAlignment->xStats.pOther.lock( )->uiBeginOnRef ) )
                flag |= NEXT_REVERSE_COMPLEMENTED;

            sContigOther = pAlignment->xStats.pOther.lock( )->getContig( *pPack );
            sPosOther = std::to_string( pAlignment->xStats.pOther.lock( )->getSamPosition( *pPack ) );

            if( !pAlignment->xStats.bFirst )
            {
                sSegment = pAlignment->getQuerySequence( *pQuery2, *pPack );
                sTlen = "-" + sTlen;
            } // if
            else
                sSegment = pAlignment->getQuerySequence( *pQuery1, *pPack );
            DEBUG( if( pQuery1->uiFromLine != pQuery2->uiFromLine ) {
                std::cerr << "outputting paired alignment for reads from different lines: " << pQuery1->uiFromLine
                          << " and " << pQuery2->uiFromLine << "; query names are: " << pQuery1->sName << " and "
                          << pQuery2->sName << std::endl;
            } // if
                   ) // DEBUG
        } // if
        else
            sSegment = pAlignment->getQuerySequence( *pQuery1, *pPack );

        std::string sRefName = pAlignment->getContig( *pPack );
        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition( *pPack );

        double fMappingQual = pAlignment->fMappingQuality;


        std::string sSQL = "";

        sSQL += "INSERT INTO alignment (cigar, position, mapping_quality, query_id, run_id, sam_flags, "
                "contig, query_sequence, contig_other, position_other, length) VALUES ( \'";

        sSQL += sCigar + "\', ";
        sSQL += std::to_string( uiRefPos ) + ", ";
        sSQL += std::to_string( fMappingQual ) + ", \'";
        sSQL += pAlignment->xStats.sName + "\', ";
        sSQL += std::to_string( iRunId ) + ", ";
        sSQL += std::to_string( flag ) + ", \'";
        sSQL += sRefName + "\', \'";
        sSQL += sSegment + "\', \'";
        sSQL += sContigOther + "\', ";
        sSQL += sPosOther + ", ";
        sSQL += sTlen + " )";

        // std::cerr << sSQL << std::endl;

        xConnection.exec( sSQL );
    } // for

    return std::shared_ptr<Container>( new Container( ) ); // return empty container
} // function

#ifdef WITH_PYTHON
void exportDBWriter( )
{
    // export the DbWriter class
    exportModule<DbWriter, std::string, uint32_t>( "DbWriter" );
    // export the PairedDbWriter class
    exportModule<PairedDbWriter, std::string, uint32_t>( "PairedDbWriter" );

} // function

#endif // WITH_PYTHON

#endif // WITH_POSTGRES