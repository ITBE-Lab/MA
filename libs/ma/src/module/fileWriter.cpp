/**
 * @file fileWriter.cpp
 * @author Markus Schmidt
 * @todo this file contains duplicate code...
 */
#include "module/fileWriter.h"

using namespace libMA;

std::shared_ptr<libMS::Container> FileWriter::execute( std::shared_ptr<NucSeq> pQuery,
                                                std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                    pAlignments,
                                                std::shared_ptr<Pack>
                                                    pPack )
{
    std::string sCombined = "";
    for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
    {
        if( pAlignment->length( ) == 0 )
            continue;
        if( bNoSecondary && pAlignment->bSecondary )
            continue;
        if( bNoSupplementary && pAlignment->bSupplementary )
            continue;

        if( bForcedConsistentConsequtiveInsertionDeletionOrder &&
            pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) )
        {
            pAlignment->invertSuccessiveInserionAndDeletion( );
        } // if

        std::string sCigar;
        if( bCGTag && pAlignment->data.size( ) >= uiMaxCigarLen )
            sCigar = std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery ).append( "S" );
        else
            sCigar = bOutputMInsteadOfXAndEqual ? pAlignment->cigarStringWithMInsteadOfXandEqual( *pPack )
                                                : pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sName = pQuery->sName;
        std::string sSegment = pAlignment->getQuerySequence( *pQuery, *pPack );
        std::string sQual = pAlignment->getQueryQuality( *pQuery );

        std::string sRefName = pAlignment->getContig( *pPack );
        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition( *pPack );

#if DEBUG_LEVEL > 0
        bool bWrong = false;
        if( pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) )
        {
            //@todo frill in this self check...
        } // if
        else
        {
            if( pAlignment->uiBeginOnRef != pPack->startOfSequenceWithName( sRefName ) + uiRefPos - 1 )
                bWrong = true;
        } // else

        if( bWrong )
        {
            std::cerr << "Error: Tried to write wrong index to file" << std::endl;
            std::cerr << "Have: " << sRefName << " (= " << pPack->startOfSequenceWithName( sRefName ) << ") "
                      << uiRefPos << std::endl;
            std::cerr << "Wanted: " << pAlignment->uiBeginOnRef << " " << uiRefPos << std::endl;
            if( pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) )
                std::cerr << "Begin is reverse: True" << std::endl;
            else
                std::cerr << "Begin is reverse: False" << std::endl;
            if( pPack->bPositionIsOnReversStrand( pAlignment->uiEndOnRef ) )
                std::cerr << "End is reverse: True" << std::endl;
            else
                std::cerr << "End is reverse: False" << std::endl;
            throw std::runtime_error( "Error: Tried to write wrong index to file" );
        } // if
#endif

        std::string sTag = this->computeTag( pQuery, pAlignment, pPack, pAlignments );

        std::string sMapQual;
        if( std::isnan( pAlignment->fMappingQuality ) )
            sMapQual = "255";
        else
            sMapQual = std::to_string( static_cast<int>( std::ceil( pAlignment->fMappingQuality * 254 ) ) );

        assert(sTag.empty() || sTag[0] == '\t');
        sCombined +=
            // query name
            sName + "\t" +
            // alignment flag
            std::to_string( flag ) + "\t" +
            // reference name
            sRefName + "\t" +
            // pos
            std::to_string( uiRefPos ) + "\t" +
            // mapping quality
            sMapQual + "\t" +
            // cigar
            sCigar + "\t" +
            // Ref. name of the mate/next read
            sContigOther + "\t" +
            // Position of the mate/next read
            sPosOther + "\t" +
            // observed Template length
            "0\t" +
            // segment sequence
            sSegment + "\t" +
            // ASCII of Phred-scaled base Quality+33
            sQual + // \t is in sTag
            // Tag
            sTag + "\n";
    } // for
    // if we have not computed any alignment then we should still output the query as unaligned:
    if( pAlignments->size( ) == 0 )
    {
        sCombined +=
            // query name
            pQuery->sName + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED ) + "\t*\t0\t255\t*\t*\t0\t0\t" +
            // segment sequence
            pQuery->toString( ) + "\t" + pQuery->toQualString() + "\n";
    } // if
    if( sCombined.size( ) == 0 )
    {
        sCombined +=
            // query name
            pQuery->sName + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED ) +
            "\t*\t0\t0\t*\t*\t0\t0\t" +
            // segment sequence
            pQuery->toString( ) + "\t" + pQuery->toQualString() + "\n";
    } // if

    { // scope xGuard
        // synchronize file output
        std::lock_guard<std::mutex> xGuard( *pLock );

        // print alignment
        // flushing will be done in the outstream class
        *pOut << sCombined;
    } // scope xGuard
    return std::shared_ptr<libMS::Container>( new libMS::Container( ) );
} // function

std::shared_ptr<libMS::Container> PairedFileWriter::execute( std::shared_ptr<NucSeq> pQuery1,
                                                      std::shared_ptr<NucSeq>
                                                          pQuery2,
                                                      std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                          pAlignments,
                                                      std::shared_ptr<Pack>
                                                          pPack )
{
    std::string sCombined = "";
    bool bFirstQueryHasAlignment = false;
    bool bSecondQueryHasAlignment = false;
    for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
    {
        if( pAlignment->length( ) == 0 )
            continue;
        if( bNoSecondary && pAlignment->bSecondary )
            continue;
        if( bNoSupplementary && pAlignment->bSupplementary )
            continue;

        if( bForcedConsistentConsequtiveInsertionDeletionOrder &&
            pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) )
        {
            pAlignment->invertSuccessiveInserionAndDeletion( );
        } // if

        if( pAlignment->xStats.bFirst )
            bFirstQueryHasAlignment = true;
        if( !pAlignment->xStats.bFirst )
            bSecondQueryHasAlignment = true;

        std::string sCigar;
        if( bCGTag && pAlignment->data.size( ) >= uiMaxCigarLen )
            sCigar = std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery ).append( "S" );
        else
            sCigar = bOutputMInsteadOfXAndEqual ? pAlignment->cigarStringWithMInsteadOfXandEqual( *pPack )
                                                : pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sName = pAlignment->xStats.bFirst ? pQuery1->sName : pQuery2->sName;
        // DEBUG( std::cout << "Aligned: " << sName << std::endl; )
        std::string sSegment = pAlignment->getQuerySequence( pAlignment->xStats.bFirst ? *pQuery1 : *pQuery2, *pPack );
        std::string sQual = pAlignment->getQueryQuality( pAlignment->xStats.bFirst ? *pQuery1 : *pQuery2 );
        // paired
        flag |= MULTIPLE_SEGMENTS_IN_TEMPLATE | SEGMENT_PROPERLY_ALIGNED;
        flag |= pAlignment->xStats.bFirst ? FIRST_IN_TEMPLATE : LAST_IN_TEMPLATE;
        std::string sTlen = "0";
        std::string sRefName = pAlignment->getContig( *pPack );
        if( pAlignment->xStats.pOther.lock( ) != nullptr )
        {

            nucSeqIndex uiP1 = pAlignment->beginOnRef( );
            // for illumina the reads are always on opposite strands
            nucSeqIndex uiP2 = pPack->uiPositionToReverseStrand( pAlignment->xStats.pOther.lock( )->beginOnRef( ) );
            // get the distance of the alignments on the reference
            nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;
            sTlen = ( pAlignment->xStats.bFirst ? "" : "-" ) + std::to_string( std::llabs( d ) );

            // assert( pQuery2 != nullptr );
            if( pPack->bPositionIsOnReversStrand( pAlignment->xStats.pOther.lock( )->uiBeginOnRef ) )
                flag |= NEXT_REVERSE_COMPLEMENTED;

            sContigOther = pAlignment->xStats.pOther.lock( )->getContig( *pPack );
            // SAM file specification:
            // This field is set as ‘*’ when the information is unavailable, and set as ‘=’ if RNEXT is identical RNAME
            if( sContigOther == sRefName )
                sContigOther = "=";
            sPosOther = std::to_string( pAlignment->xStats.pOther.lock( )->getSamPosition( *pPack ) );

#if DEBUG_LEVEL > 0
            if( pQuery1->uiFromLine != pQuery2->uiFromLine )
            {
                std::cerr << "outputting paired alignment for reads from different lines: " << pQuery1->uiFromLine
                          << " and " << pQuery2->uiFromLine << "; query names are: " << pQuery1->sName << " and "
                          << pQuery2->sName << std::endl;
            } // if
#endif
        } // if
#if DEBUG_LEVEL > 0
        else if( !pAlignment->bSecondary && !pAlignment->bSupplementary )
            std::cerr << "[Info] read " << pQuery1->sName << " is unpaired." << std::endl;
#endif

        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition( *pPack );

#if DEBUG_LEVEL > 0
        bool bWrong = false;
        if( pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) )
        {
            //@todo frill in this self check...
        } // if
        else
        {
            if( pAlignment->uiBeginOnRef != pPack->startOfSequenceWithName( sRefName ) + uiRefPos - 1 )
                bWrong = true;
        } // else

        if( bWrong )
        {
            std::cerr << "Error: Tried to write wrong index to file" << std::endl;
            std::cerr << "Have: " << sRefName << " (= " << pPack->startOfSequenceWithName( sRefName ) << ") "
                      << uiRefPos << std::endl;
            std::cerr << "Wanted: " << pAlignment->uiBeginOnRef << " " << uiRefPos << std::endl;
            if( pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) )
                std::cerr << "Begin is reverse: True" << std::endl;
            else
                std::cerr << "Begin is reverse: False" << std::endl;
            if( pPack->bPositionIsOnReversStrand( pAlignment->uiEndOnRef ) )
                std::cerr << "End is reverse: True" << std::endl;
            else
                std::cerr << "End is reverse: False" << std::endl;
            throw std::runtime_error( "Error: Tried to write wrong index to file" );
        } // if
#endif
        std::string sTag =
            this->computeTag( pAlignment->xStats.bFirst ? pQuery1 : pQuery2, pAlignment, pPack, pAlignments );

        std::string sMapQual;
        if( std::isnan( pAlignment->fMappingQuality ) )
            sMapQual = "255";
        else
        {
            assert( pAlignment->fMappingQuality >= 0 );
            sMapQual =
                std::to_string( std::min( static_cast<int>( std::ceil( pAlignment->fMappingQuality * 254 ) ), 255 ) );
        }

        assert(sTag.empty() || sTag[0] == '\t');
        sCombined +=
            // query name
            sName + "\t" +
            // alignment flag
            std::to_string( flag ) + "\t" +
            // reference name
            sRefName + "\t" +
            // pos
            std::to_string( uiRefPos ) + "\t" +
            // mapping quality
            sMapQual + "\t" +
            // cigar
            sCigar + "\t" +
            // Ref. name of the mate/next read
            sContigOther + "\t" +
            // Position of the mate/next read
            sPosOther + "\t" +
            // observed Template length
            /*sTlen*/ "0" + "\t" + // output information unavailable for now...
            // segment sequence
            sSegment + "\t" +
            // ASCII of Phred-scaled base Quality+33
            sQual + // \t is in sTag
            // Tag
            sTag + "\n";
    } // for
    // if we have not computed any alignment then we should still output the query as unaligned:
    if( !bFirstQueryHasAlignment && !bSecondQueryHasAlignment )
    {
        sCombined +=
            // query name
            pQuery1->sName + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED | MULTIPLE_SEGMENTS_IN_TEMPLATE | FIRST_IN_TEMPLATE |
                            NEXT_SEGMENT_UNMAPPED ) +
            "\t*\t0\t0\t*\t*\t0\t0\t" +
            // segment sequence
            pQuery1->toString( ) + "\t" + pQuery1->toQualString() + "\n";
        sCombined +=
            // query name
            pQuery2->sName + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED | MULTIPLE_SEGMENTS_IN_TEMPLATE | LAST_IN_TEMPLATE |
                            NEXT_SEGMENT_UNMAPPED ) +
            "\t*\t0\t0\t*\t*\t0\t0\t" +
            // segment sequence
            pQuery2->toString( ) + "\t" + pQuery2->toQualString() + "\n";
    } // if
    // if we have not computed an alignment for only one of the queries output the other one:
    else if( !bFirstQueryHasAlignment || !bSecondQueryHasAlignment )
    {
        // both asserts are guaranteed due to bFirstQueryHasAlignment and bSecondQueryHasAlignment
        assert( pAlignments->size( ) > 0 );
        assert( ( *pAlignments )[ 0 ]->xStats.bFirst == ( !bFirstQueryHasAlignment ? false : true ) );

        std::string sPosOther = std::to_string( ( *pAlignments )[ 0 ]->getSamPosition( *pPack ) );
        std::string sContigOther = ( *pAlignments )[ 0 ]->getContig( *pPack );
        sCombined +=
            // query name
            ( !bFirstQueryHasAlignment ? pQuery1->sName : pQuery2->sName ) + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED | MULTIPLE_SEGMENTS_IN_TEMPLATE |
                            ( !bFirstQueryHasAlignment ? FIRST_IN_TEMPLATE : LAST_IN_TEMPLATE ) ) +
            "\t" + sContigOther + "\t" + sPosOther +
            "\t0" // outputing mapq0 because bwa mem does that as well...
            "\t*\t=\t" +
            sPosOther + "\t0\t" +
            // segment sequence
            ( !bFirstQueryHasAlignment ? pQuery1->toString( ) : pQuery2->toString( ) ) + "\t*\n";
    } // if

    if( sCombined.size( ) > 0 )
    { // scope xGuard
        // synchronize file output
        std::lock_guard<std::mutex> xGuard( *pLock );

        // print alignment
        // flushing will be done in the the outstream class
        *pOut << sCombined;
    } // if & scope xGuard
    return std::shared_ptr<libMS::Container>( new libMS::Container( ) );
} // function

#ifdef WITH_PYTHON
void exportFileWriter( py::module& rxPyModuleId )
{
    // export the FileWriter class
    exportModule<FileWriter, std::string, std::shared_ptr<Pack>>( rxPyModuleId, "FileWriter" );
    exportModuleAlternateConstructor<FileWriter, std::shared_ptr<FileWriter>>( rxPyModuleId, "SyncFileWriter" );

    // export the PairedFileWriter class
    exportModule<PairedFileWriter, std::string, std::shared_ptr<Pack>>( rxPyModuleId, "PairedFileWriter" );
    exportModuleAlternateConstructor<PairedFileWriter, std::shared_ptr<FileWriter>>( rxPyModuleId,
                                                                                     "SyncPairedFileWriter" );
    exportModuleAlternateConstructor<PairedFileWriter, std::shared_ptr<PairedFileWriter>>(
        rxPyModuleId, "PairedSyncPairedFileWriter" );
} // function
#endif