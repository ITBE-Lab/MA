/**
 * @file fileWriter.h
 * @brief Writes alignments to a file.
 * @author Markus Schmidt
 */
#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include "container/alignment.h"
#include "module/module.h"

namespace libMA
{

/**
 * @brief wrapper for various out streams.
 * @details
 * this exists, so that alignments can be written to a GUI element if so desired.
 */
class OutStream
{
  public:
    virtual OutStream& operator<<( std::string )
    {
        return *this;
    };
}; // class

/**
 * @brief wraps the std outstream.
 */
class StdOutStream : public OutStream
{
  public:
    StdOutStream& operator<<( std::string s )
    {
        std::cout << s << std::flush;
        return *this;
    } // function
}; // class

/**
 * @brief wraps a file outstream.
 * @details
 * truncates the file if it already exists
 */
class FileOutStream : public OutStream
{
  public:
    std::ofstream file;

    FileOutStream( std::string sFileName ) : file( sFileName, std::ofstream::out | std::ofstream::trunc )
    {
        if( !file.good( ) )
        {
            throw AnnotatedException( "Unable to open file" + sFileName );
        } // if
    } // constructor

    ~FileOutStream( )
    {
        file.close( );
    } // deconstructor

    FileOutStream& operator<<( std::string s )
    {
        file << s << std::flush;
        return *this;
    } // function
}; // class

class TagGenerator
{
  public:
    const bool bMDTag; // NGMLR SAM emulation
    const bool bSVTag; // NGMLR SAM emulation
    const bool bNMTag; // NGMLR SAM emulation
    // according to NGMLR's documentation this is it's behaviour; according to code though it's not...
    const bool bNMTagDoNOTCountIndels = false;
    const bool bASTag;
    const bool bOutputMInsteadOfXAndEqual;
    const bool bXITag; // NGMLR SAM emulation
    const bool bXETag; // NGMLR SAM emulation
    const bool bXRTag; // NGMLR SAM emulation
    const bool bQS_QETag; // NGMLR SAM emulation
    const bool bCVTag; // NGMLR SAM emulation
    const bool bSATag; // NGMLR SAM emulation
    const bool bCGTag; // NGMLR SAM emulation
    const size_t uiMaxCigarLen = 0x10000; // NGMLR SAM emulation
    const bool bForcedConsistentConsequtiveInsertionDeletionOrder; // NGMLR SAM emulation

    TagGenerator( const ParameterSetManager& rParameters )
        : bMDTag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bSVTag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bNMTag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bASTag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bOutputMInsteadOfXAndEqual( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ||
                                      rParameters.getSelected( )->xOutputMCigar->get( ) ),
          bXITag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bXETag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bXRTag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bQS_QETag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bCVTag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bSATag( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) ),
          bCGTag( rParameters.getSelected( )->xCGTag->get( ) ),
          bForcedConsistentConsequtiveInsertionDeletionOrder( rParameters.getSelected( )->xEmulateNgmlrTags->get( ) )
    {} // constructor

    std::string computeTag( const std::shared_ptr<NucSeq> pQuery,
                            const std::shared_ptr<Alignment>
                                pAlignment,
                            const std::shared_ptr<Pack>
                                pPack,
                            const std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                pvAllAlignments ) const
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
                        for( const char& rC : pPack->vExtract( uiRPos, uiRPos + xTup.second )->toString( ) )
                        {
                            if( bFirst == true )
                                bFirst = false;
                            else
                                sTag.push_back( '0' );
                            sTag.push_back( rC );
                        } // for
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
             * alignment. Sniffles uses soft-clip information to improve detection of (large) insertions. As a
             * consequence longer regions of Ns in the reference would cause false positives signals as the reads
             * up/down stream of the Ns will all be soft-clipped at the same position of the reference. `0x1` helped to
             * avoid that. `0x2` is set if > 95 % of the read base pairs were aligned in the particular alignment.
             * Sniffles could compute that from the CIGAR string but since we had the SV tag already at that point we
             * just added it. Don't know if Fritz is still using it.
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
        if( bASTag )
        {
            // the alignment score
            sTag.append( "\tAS:i:" ).append( std::to_string( pAlignment->score( ) ) );
        } // if
        if( bNMTag )
        {
            // see function pAlignment->getNumDifferences
            sTag.append( "\tNM:i:" )
                .append( std::to_string( pAlignment->getNumDifferences( pPack, bNMTagDoNOTCountIndels ) ) );
        } // if
        if( bXITag )
        {
            /*
             * XI = Identity of the alignment: https://github.com/Cibiv/NextGenMap/wiki/Documentation#custom-sam-tags
             *
             * Explanation of identity
             * https://www.researchgate.net/post/Homology_similarity_and_identity-can_anyone_help_with_these_terms
             * A: AAGGCTT; B: AAGGC -> Here identity(A,B)=100% (5 identical nucleotides / min(length(A),length(B))).
             */
            float fIdentity =
                pAlignment->getNumMatches( ) / (float)std::min( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery,
                                                                pAlignment->uiEndOnRef - pAlignment->uiBeginOnRef );
            sTag.append( "\tXI:f:" ).append( std::to_string( fIdentity ) );
        } // if
        if( bXETag )
        {
            /*
             * XS = Number of supported seeds according to here:
             * https://github.com/Cibiv/NextGenMap/wiki/Documentation#custom-sam-tags
             *
             * @note: actually NGMLR outputs the alignment score as the XE tag
             */
            sTag.append( "\tXE:i:" ).append( std::to_string( pAlignment->score( ) ) );
        } // if
        if( bXRTag )
        {
            /*
             * XR is length of alignment on query
             */
            sTag.append( "\tXR:i:" ).append( std::to_string( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery ) );
        } // if
        if( bCVTag )
        {
            /*
             * XR is length of alignment on query
             */
            float fCoverage =
                100.0f * ( pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery ) / (float)pQuery->length( );
            sTag.append( "\tCV:f:" ).append( std::to_string( fCoverage ) );
        } // if
        if( bSATag && pvAllAlignments->size( ) > 1 )
        {
            /*
             * For all chimeric alignments (supplementary alignments) separated additional output
             * (http://samtools.github.io/hts-specs/SAMtags.pdf)
             *
             * Conventionally, at a supplementary line, the first element points to the primary line:
             * We hold this convention implicitly, since our alignments are sorted by score...
             */
            bool bFoundSupplementarySisters = false;
            std::string sSATag;

            for( auto pOtherAlignment : *pvAllAlignments )
            {
                if( pAlignment == pOtherAlignment )
                    continue;
                if( pOtherAlignment->bSecondary )
                    continue;
                if( pOtherAlignment->xStats.bFirst != pAlignment->xStats.bFirst )
                    continue;
                bFoundSupplementarySisters = true;
                sSATag.append( pOtherAlignment->getContig( *pPack ) ).append( "," );
                sSATag.append( std::to_string( pOtherAlignment->getSamPosition( *pPack ) ) ).append( "," );
                sSATag.append( ( pPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef ) ? "-," : "+," ) );
                sSATag
                    .append( bOutputMInsteadOfXAndEqual ? pOtherAlignment->cigarStringWithMInsteadOfXandEqual( *pPack )
                                                        : pOtherAlignment->cigarString( *pPack ) )
                    .append( "," );

                std::string sMapQual;
                if( std::isnan( pOtherAlignment->fMappingQuality ) )
                    sMapQual = "255";
                else
                    sMapQual =
                        std::to_string( static_cast<int>( std::ceil( pOtherAlignment->fMappingQuality * 254 ) ) );
                sSATag.append( sMapQual ).append( "," );
                sSATag.append( std::to_string( pAlignment->getNumDifferences( pPack, bNMTagDoNOTCountIndels ) ) );

                sSATag.append( ";" );
            } // for

            if( bFoundSupplementarySisters )
                sTag.append( "\tSA:Z:" ).append( sSATag );
        } // if
        if( bQS_QETag )
        {
            /*
             * Output start and end of alignment on read
             */
            sTag.append( "\tQS:i:" )
                .append( std::to_string( pAlignment->uiBeginOnQuery ) )
                .append( "\tQE:i:" )
                .append( std::to_string( pAlignment->uiEndOnQuery ) );
        } // if
        // check if the total number of CIGAR operations is too large
        if( bCGTag && pAlignment->data.size( ) >= uiMaxCigarLen )
        {
            sTag.append( "\tCG:B:I" );
            for( std::pair<MatchType, nucSeqIndex>& rPair : pAlignment->data )
            {
                uint32_t uiOperation = 0;
                switch( rPair.first )
                {
                    case MatchType::seed:
                    case MatchType::match:
                        uiOperation = 7;
                        break;
                    case MatchType::missmatch:
                        uiOperation = 8;
                        break;
                    case MatchType::insertion:
                        uiOperation = 1;
                        break;
                    case MatchType::deletion:
                        uiOperation = 2;
                        break;
                    default:
                        std::cerr << "WARNING invalid cigar symbol" << std::endl;
                        break;
                } // switch

                uint32_t uiOut = ( uint32_t )( rPair.second << 4 ) | uiOperation;
                sTag.append( "," ).append( std::to_string( uiOut ) );
            } // for
        } // if
        return sTag;
    } // method
}; // class

/**
 * @brief Writes SAM output.
 * @note flushing of the outstream; this must be done in the deconstructor of OutStream
 *
 */
class FileWriter : public Module<Container, false, NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>,
                   private TagGenerator
{
  public:
    // holds a file ourstream if necessary
    std::shared_ptr<OutStream> pOut;
    std::shared_ptr<std::mutex> pLock;
    const bool bNoSecondary;
    const bool bNoSupplementary;

    /**
     * @brief creates a new FileWriter.
     * @details
     * if sFileName is "stdout" the writer will output to stdout instead of the file.
     * Otherwise sFileName is used as the filename to write to.
     * The file will be truncated is it already exists.
     */
    FileWriter( const ParameterSetManager& rParameters, std::string sFileName, std::shared_ptr<Pack> pPackContainer )
        : TagGenerator( rParameters ),
          pLock( new std::mutex ),
          bNoSecondary( rParameters.getSelected( )->xNoSecondary->get( ) ),
          bNoSupplementary( rParameters.getSelected( )->xNoSupplementary->get( ) )
    {
        if( sFileName != "stdout" )
            pOut = std::shared_ptr<OutStream>( new FileOutStream( sFileName ) );
        else
            pOut = std::shared_ptr<OutStream>( new StdOutStream( ) );
        //*pOut << "@HD VN:1.5 SO:unknown\n";
        for( auto& rSeqInPack : pPackContainer->xVectorOfSequenceDescriptors )
        {
            *pOut << "@SQ\tSN:" << rSeqInPack.sName << "\tLN:" << std::to_string( rSeqInPack.uiLengthUnpacked ) << "\n";
        } // for
        *pOut << "@PG\tID:ma\tPN:ma\tVN:0.1.0\tCL:na\n";
    } // constructor

    /**
     * @brief creates a new FileWriter.
     * @details
     * Allows more control of the output by using a given OutStream.
     */
    FileWriter( const ParameterSetManager& rParameters, std::shared_ptr<OutStream> pOut,
                std::shared_ptr<Pack> pPackContainer )
        : TagGenerator( rParameters ),
          pOut( pOut ),
          pLock( new std::mutex ),
          bNoSecondary( rParameters.getSelected( )->xNoSecondary->get( ) ),
          bNoSupplementary( rParameters.getSelected( )->xNoSupplementary->get( ) )
    {
        //*pOut << "@HD VN:1.5 SO:unknown\n";
        for( auto& rSeqInPack : pPackContainer->xVectorOfSequenceDescriptors )
        {
            *pOut << "@SQ\tSN:" << rSeqInPack.sName << " LN:" << std::to_string( rSeqInPack.uiLengthUnpacked ) << "\n";
        } // for
        *pOut << "@PG\tID:ma\tPN:ma\tVN:0.1.0\tCL:na\n";
    } // constructor

    virtual std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<NucSeq> pQuery,
                                                         std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                             pAlignments,
                                                         std::shared_ptr<Pack>
                                                             pPack );

}; // class

/**
 * @brief Writes SAM output.
 * @note flushing of the outstream; this must be done in the deconstructor of OutStream
 *
 */
class PairedFileWriter
    : public Module<Container, false, NucSeq, NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>,
      private TagGenerator
{
  public:
    // holds a file ourstream if necessary
    std::shared_ptr<OutStream> pOut;
    std::shared_ptr<std::mutex> pLock;
    const bool bNoSecondary;
    const bool bNoSupplementary;

    /**
     * @brief creates a new FileWriter.
     * @details
     * if sFileName is "stdout" the writer will output to stdout instead of the file.
     * Otherwise sFileName is used as the filename to write to.
     * The file will be truncated is it already exists.
     */
    PairedFileWriter( const ParameterSetManager& rParameters, std::string sFileName,
                      std::shared_ptr<Pack> pPackContainer )
        : TagGenerator( rParameters ),
          pLock( new std::mutex ),
          bNoSecondary( rParameters.getSelected( )->xNoSecondary->get( ) ),
          bNoSupplementary( rParameters.getSelected( )->xNoSupplementary->get( ) )
    {
        if( sFileName != "stdout" )
            pOut = std::shared_ptr<OutStream>( new FileOutStream( sFileName ) );
        else
            pOut = std::shared_ptr<OutStream>( new StdOutStream( ) );
        //*pOut << "@HD VN:1.5 SO:unknown\n";
        for( auto& rSeqInPack : pPackContainer->xVectorOfSequenceDescriptors )
        {
            *pOut << "@SQ\tSN:" << rSeqInPack.sName << "\tLN:" << std::to_string( rSeqInPack.uiLengthUnpacked ) << "\n";
        } // for
        *pOut << "@PG\tID:ma\tPN:ma\tVN:0.1.0\tCL:na\n";
    } // constructor

    /**
     * @brief creates a new FileWriter.
     * @details
     * Allows more control of the output by using a given OutStream.
     */
    PairedFileWriter( const ParameterSetManager& rParameters, std::shared_ptr<OutStream> pOut,
                      std::shared_ptr<Pack> pPackContainer )
        : TagGenerator( rParameters ),
          pOut( pOut ),
          pLock( new std::mutex ),
          bNoSecondary( rParameters.getSelected( )->xNoSecondary->get( ) ),
          bNoSupplementary( rParameters.getSelected( )->xNoSupplementary->get( ) )
    {
        //*pOut << "@HD VN:1.5 SO:unknown\n";
        for( auto& rSeqInPack : pPackContainer->xVectorOfSequenceDescriptors )
        {
            *pOut << "@SQ\tSN:" << rSeqInPack.sName << " LN:" << std::to_string( rSeqInPack.uiLengthUnpacked ) << "\n";
        } // for
        *pOut << "@PG\tID:ma\tPN:ma\tVN:0.1.0\tCL:na\n";
    } // constructor

    virtual std::shared_ptr<Container>
        EXPORTED execute( std::shared_ptr<NucSeq> pQuery1, std::shared_ptr<NucSeq> pQuery2,
                          std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments,
                          std::shared_ptr<Pack> pPack );

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportFileWriter( );
#else
void exportFileWriter( py::module& rxPyModuleId );
#endif
#endif

#endif