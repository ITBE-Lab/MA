/**
 * @file fileWriter.cpp
 * @author Markus Schmidt
 * @todo this file contains duplicate code...
 */
#include "module/fileWriter.h"

using namespace libMA;

std::string computeTag( const std::shared_ptr<NucSeq> pQuery,
                        const std::shared_ptr<Alignment>
                            pAlignment,
                        const std::shared_ptr<Pack>
                            pPack,
                        const bool bMDTag,
                        const bool bSVTag )
{
    std::string sTag = "";
    if( bMDTag )
    {
        sTag.append( "\tMD:Z:" );
        nucSeqIndex uiRPos = pAlignment->uiBeginOnRef;
        nucSeqIndex uiNumMatchesAndSeeds = 0;
        bool bLastWasDeletion = false;
        for( std::pair<MatchType, nucSeqIndex> xTup : pAlignment->data )
        {
            /*
             * If seeds are followed my matches or vice versa, we have to print the accumulated number.
             * I.e a 5nt seed followed by 2 matches is NOT: 52 but rather 7....
             * Therefore we accumulate the number of matches and seed in uiNumMatchesAndSeeds and print it
             * once the next element is neither match nor seed.
             */
            if( ( xTup.first == MatchType::missmatch || xTup.first == MatchType::deletion ) &&
                uiNumMatchesAndSeeds > 0 )
            {
                // append the number of seeds and matches
                sTag.append( std::to_string( uiNumMatchesAndSeeds ) );
                uiNumMatchesAndSeeds = 0;
            } // if
            bool bFirst = !bLastWasDeletion;
            bLastWasDeletion = false;
            switch( xTup.first )
            {
                case MatchType::match:
                case MatchType::seed:
                    // see comment above
                    uiNumMatchesAndSeeds += xTup.second;
                    // update reference position
                    uiRPos += xTup.second;
                    break;
                case MatchType::insertion:
                    // this is simply ignored since it does not represent a loss of information...
                    break;
                case MatchType::missmatch:
                    // append all missmatched nucleotides with zeros in between them.
                    for(const char& rC : pPack->vExtract( uiRPos, uiRPos + xTup.second )->toString( ) )
                    {
                        if(bFirst == true)
                            bFirst = false;
                        else
                            sTag.push_back('0');
                        sTag.push_back(rC);
                    }// for
                    // update reference position
                    uiRPos += xTup.second;
                    break;
                case MatchType::deletion:
                    // prefix for deletion
                    sTag.append( "^" );
                    // deleted sequence
                    sTag.append( pPack->vExtract( uiRPos, uiRPos + xTup.second )->toString( ) );
                    // update reference position
                    uiRPos += xTup.second;
                    bLastWasDeletion = true;
                    break;
                default:
                    throw std::runtime_error( "Invalid symbol in cigar!" );
                    break;
            } // switch
        } // for
        if( uiNumMatchesAndSeeds > 0 )
            sTag.append( std::to_string( uiNumMatchesAndSeeds ) );
    } // if
    if( bSVTag )
    {
        /*
         * @see https://github.com/fritzsedlazeck/Sniffles/issues/51#issuecomment-377471553
         * "
         * Hi Heng,
         *
         * The SV:i tag in the NGMLR output contains two bit flags:
         *
         * `0x1` indicates whether the reference sequence consists of mostly Ns up and/or downstream of the read
         * alignment. Sniffles uses soft-clip information to improve detection of (large) insertions. As a consequence
         * longer regions of Ns in the reference would cause false positives signals as the reads up/down stream of the
         * Ns will all be soft-clipped at the same position of the reference. `0x1` helped to avoid that. `0x2` is set
         * if > 95 % of the read base pairs were aligned in the particular alignment. Sniffles could compute that from
         * the CIGAR string but since we had the SV tag already at that point we just added it. Don't know if Fritz is
         * still using it.
         *
         * Hope that helps,
         * Philipp
         * "
         * For 0x1 we check if either 80% of the 100 nt in front or behind the alignment are covered by a hole.
         */
        size_t uiTag = 0;
        if( pPack->amountOfRegionCoveredByHole( pAlignment->uiBeginOnRef - 100, pAlignment->uiBeginOnRef ) > .8 ||
            pPack->amountOfRegionCoveredByHole( pAlignment->uiEndOnRef, pAlignment->uiEndOnRef + 100 ) > .8 )
            uiTag += 1;
        if( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery >= pQuery->length( ) * 0.95 )
            uiTag += 2;
        sTag.append( "\tSV:i:" ).append( std::to_string( uiTag ) );
    } // if
    return sTag;
} // function

std::shared_ptr<Container> FileWriter::execute( std::shared_ptr<NucSeq> pQuery,
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
        std::string sCigar = pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sName = pQuery->sName;
        std::string sSegment = pAlignment->getQuerySequence( *pQuery, *pPack );

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
            throw AnnotatedException( "Error: Tried to write wrong index to file" );
        } // if
#endif

        std::string sTag = computeTag( pQuery, pAlignment, pPack, bMDTag, bSVTag );

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
            "0\t" +
            // segment sequence
            sSegment + "\t" +
            // ASCII of Phred-scaled base Quality+33
            "*"
            // Tag
            + sTag + "\n";
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
            pQuery->toString( ) + "\t*\n";
    } // if

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
        if( pAlignment->xStats.bFirst )
            bFirstQueryHasAlignment = true;
        if( !pAlignment->xStats.bFirst )
            bSecondQueryHasAlignment = true;
        std::string sCigar = pAlignment->cigarString( *pPack );

        uint32_t flag = pAlignment->getSamFlag( *pPack );

        std::string sContigOther = "*";
        std::string sPosOther = "0";
        std::string sName = pAlignment->xStats.bFirst ? pQuery1->sName : pQuery2->sName;
        // DEBUG( std::cout << "Aligned: " << sName << std::endl; )
        std::string sSegment = pAlignment->getQuerySequence( pAlignment->xStats.bFirst ? *pQuery1 : *pQuery2, *pPack );
        // paired
        flag |= MULTIPLE_SEGMENTS_IN_TEMPLATE | SEGMENT_PROPERLY_ALIGNED;
        std::string sTlen = "0";
        if( pAlignment->xStats.pOther.lock( ) != nullptr )
        {

            nucSeqIndex uiP1 = pAlignment->beginOnRef( );
            // for illumina the reads are always on opposite strands
            nucSeqIndex uiP2 = pPack->uiPositionToReverseStrand( pAlignment->xStats.pOther.lock( )->beginOnRef( ) );
            // get the distance of the alignments on the reference
            nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;
            sTlen = ( pAlignment->xStats.bFirst ? "" : "-" ) + std::to_string( std::llabs( d ) );

            // assert( pQuery2 != nullptr );
            flag |= pAlignment->xStats.bFirst ? FIRST_IN_TEMPLATE : LAST_IN_TEMPLATE;
            if( pPack->bPositionIsOnReversStrand( pAlignment->xStats.pOther.lock( )->uiBeginOnRef ) )
                flag |= NEXT_REVERSE_COMPLEMENTED;

            sContigOther = pAlignment->xStats.pOther.lock( )->getContig( *pPack );
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
            throw AnnotatedException( "Error: Tried to write wrong index to file" );
        } // if
#endif
        std::string sTag =
            computeTag( pAlignment->xStats.bFirst ? pQuery1 : pQuery2, pAlignment, pPack, bMDTag, bSVTag );

        std::string sMapQual;
        if( std::isnan( pAlignment->fMappingQuality ) )
            sMapQual = "255";
        else
        {
            assert( pAlignment->fMappingQuality >= 0 );
            sMapQual =
                std::to_string( std::min( static_cast<int>( std::ceil( pAlignment->fMappingQuality * 254 ) ), 255 ) );
        }

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
            "*"
            // Tag
            + sTag + "\n";
    } // for
    // if we have not computed any alignment then we should still output the query as unaligned:
    if( !bFirstQueryHasAlignment )
    {
        sCombined +=
            // query name
            pQuery1->sName + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED | MULTIPLE_SEGMENTS_IN_TEMPLATE ) + "\t*\t0\t255\t*\t*\t0\t0\t" +
            // segment sequence
            pQuery1->toString( ) + "\t*\n";
    } // if
    // if we have not computed any alignment then we should still output the query as unaligned:
    if( !bSecondQueryHasAlignment )
    {
        sCombined +=
            // query name
            pQuery2->sName + "\t" +
            // alignment flag
            std::to_string( SEGMENT_UNMAPPED | MULTIPLE_SEGMENTS_IN_TEMPLATE ) + "\t*\t0\t255\t*\t*\t0\t0\t" +
            // segment sequence
            pQuery2->toString( ) + "\t*\n";
    } // if

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

#ifdef BOOST_PYTHON
void exportFileWriter( )
{
    // export the FileWriter class
    exportModule<FileWriter, std::string, std::shared_ptr<Pack>>( "FileWriter" );

    // export the PairedFileWriter class
    exportModule<PairedFileWriter, std::string, std::shared_ptr<Pack>>( "PairedFileWriter" );
} // function
#else
void exportFileWriter( py::module& rxPyModuleId )
{
    // export the FileWriter class
    exportModule<FileWriter, std::string, std::shared_ptr<Pack>>( rxPyModuleId, "FileWriter" );

    // export the PairedFileWriter class
    exportModule<PairedFileWriter, std::string, std::shared_ptr<Pack>>( rxPyModuleId, "PairedFileWriter" );
} // function
#endif
#endif