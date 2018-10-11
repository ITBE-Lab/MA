/**
 * @file fileWriter.cpp
 * @author Markus Schmidt
 */
#include "module/fileWriter.h"

using namespace libMA;


std::shared_ptr<Container> FileWriter::execute( std::shared_ptr<NucSeq> pQuery,
                                                std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                    pAlignments,
                                                std::shared_ptr<Pack>
                                                    pPack )
{
    std::string sCombined = "";
    for( std::shared_ptr<Alignment> pA : *pAlignments )
    {
        std::shared_ptr<Alignment> pAlignment = std::dynamic_pointer_cast<Alignment>( pA ); // dc
        if( pAlignment->length( ) == 0 )
            continue;
        std::string sCigar = pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sName = pQuery->sName;
        std::string sSegment = pAlignment->getQuerySequence( *pQuery, *pPack );
        std::string sTlen = std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery );

        std::string sRefName = pAlignment->getContig( *pPack );
        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition( *pPack ) + 1;

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
            throw AnnotatedException( "Error: Tried to write wrong index to file" );
        } // if
#endif

        std::string sMapQual;
        if( std::isnan( pAlignment->fMappingQuality ) )
            sMapQual = "255";
        else
            sMapQual = std::to_string( static_cast<int>( std::ceil( pAlignment->fMappingQuality * 254 ) ) );

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
            sTlen + "\t" +
            // segment sequence
            sSegment +
            "\t"
            // ASCII of Phred-scaled base Quality+33
            + "*\n";
    } // for

    if( sCombined.size( ) > 0 )
    { // scope xGuard
        // synchronize file output
        std::lock_guard<std::mutex> xGuard( *pLock );

        // print alignment
        // flushing will be done in the deconstructor
        *pOut << sCombined;
    } // if & scope xGuard
    return std::shared_ptr<Container>( new Container( ) );
} // function

std::shared_ptr<Container> PairedFileWriter::execute( std::shared_ptr<NucSeq> pQuery1,
                                                      std::shared_ptr<NucSeq>
                                                          pQuery2,
                                                      std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                          pAlignments,
                                                      std::shared_ptr<Pack>
                                                          pPack )
{
    std::string sCombined = "";
    for( std::shared_ptr<Alignment> pA : *pAlignments )
    {
        std::shared_ptr<Alignment> pAlignment = std::dynamic_pointer_cast<Alignment>( pA ); // dc
        if( pAlignment->length( ) == 0 )
            continue;
        std::string sCigar = pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sName = pQuery1->sName;
        // DEBUG( std::cout << "Aligned: " << sName << std::endl; )
        std::string sSegment = pAlignment->getQuerySequence( *pQuery1, *pPack );
        std::string sTlen = std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery );
        // paired
        if( pAlignment->xStats.pOther.lock( ) != nullptr )
        {
            // assert( pQuery2 != nullptr );
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
                sName = pQuery2->sName;
                sTlen = "-" + sTlen;
            } // if

#if DEBUG_LEVEL > 0
            if( pQuery1->uiFromLine != pQuery2->uiFromLine )
            {
                std::cerr << "outputting paired alignment for reads from different lines: " << pQuery1->uiFromLine
                          << " and " << pQuery2->uiFromLine << "; query names are: " << pQuery1->sName << " and "
                          << pQuery2->sName << std::endl;
            } // if
#endif
        } // if

        std::string sRefName = pAlignment->getContig( *pPack );
        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition( *pPack ) + 1;

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
            throw AnnotatedException( "Error: Tried to write wrong index to file" );
        } // if
#endif

        std::string sMapQual;
        if( std::isnan( pAlignment->fMappingQuality ) )
            sMapQual = "255";
        else
            sMapQual = std::to_string( static_cast<int>( std::ceil( pAlignment->fMappingQuality * 254 ) ) );

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
            sTlen + "\t" +
            // segment sequence
            sSegment + "\t" +
            // ASCII of Phred-scaled base Quality+33
            "*\n";
    } // for

    if( sCombined.size( ) > 0 )
    { // scope xGuard
        // synchronize file output
        std::lock_guard<std::mutex> xGuard( *pLock );

        // print alignment
        // flushing will be done in the deconstructor
        *pOut << sCombined;
    } // if & scope xGuard
    return std::shared_ptr<Container>( new Container( ) );
} // function

#ifdef WITH_PYTHON
void exportFileWriter( )
{
    // export the FileWriter class
    exportModule<FileWriter, std::string, std::shared_ptr<Pack>>( "FileWriter" );

    // export the PairedFileWriter class
    exportModule<PairedFileWriter, std::string, std::shared_ptr<Pack>>( "PairedFileWriter" );
} // function
#endif