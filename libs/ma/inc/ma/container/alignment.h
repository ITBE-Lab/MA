/**
 * @file alignment.h
 * @brief Implements a Container that holds a finished alignment.
 * @author Markus Schmidt
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "ma/container/segment.h"
#include <sstream>

#define MULTIPLE_SEGMENTS_IN_TEMPLATE 0x001
#define SEGMENT_PROPERLY_ALIGNED 0x002
#define SEGMENT_UNMAPPED 0x004
#define NEXT_SEGMENT_UNMAPPED 0x008
#define REVERSE_COMPLEMENTED 0x010
#define NEXT_REVERSE_COMPLEMENTED 0x020
#define FIRST_IN_TEMPLATE 0x040
#define LAST_IN_TEMPLATE 0x080
#define SECONDARY_ALIGNMENT 0x100
#define NOT_PASSING_FILTERS 0x200
#define PCR_OR_DUPLICATE 0x400
#define SUPPLEMENTARY_ALIGNMENT 0x800

namespace libMA
{
/**
 * @brief Describes the type of match at one specific position of the alignment.
 * @details
 * @li @c match: query and reference have the same nucleotide.
 * @li @c seed: query and reference have the same nucleotide
 * (the match was found as part  of a seed).
 * @li @c missmatch: query and reference have different nucleotides,
 *      but they are aligned to the same position nonetheless.
 * @li @c insertion: a nucleotide is present on the query that has no counterpart on the reference.
 * @li @c deletion: a nucleotide is present on the reference that has no counterpart on the query.
 */
enum MatchType
{
    seed,
    match,
    missmatch,
    insertion,
    deletion
}; // enum

/**
 * @brief Holds a finished alignment.
 * @details
 * Contains a list of MatchTypes (match, missmatch, insertion, deletion).
 * @ingroup container
 */
class Alignment : public libMS::Container
{
  public:
#if DEBUG_LEVEL >= 1
    std::vector<std::pair<nucSeqIndex, nucSeqIndex>> vGapsScatter;
#endif
    /// The sparse list of MatchTypes that describe the alignment.
    std::vector<std::pair<MatchType, nucSeqIndex>> data;
    /// The length of the alignment.
    nucSeqIndex uiLength;
    /// The start of the alignment on the reference sequence.
    nucSeqIndex uiBeginOnRef;
    /// The end of the alignment on the reference sequence.
    nucSeqIndex uiEndOnRef;
    /// The start of the alignment on the query sequence.
    nucSeqIndex uiBeginOnQuery;
    /// The end of the alignment on the query sequence.
    nucSeqIndex uiEndOnQuery;

    int64_t iScore = 0;
    double fMappingQuality = NAN;

    // some statistics
    AlignmentStatistics xStats;
    bool bSecondary;
    bool bSupplementary;

    /**
     * @brief Creates an empty alignment.
     */
    Alignment( )
        : data( ),
          uiLength( 0 ),
          uiBeginOnRef( 0 ),
          uiEndOnRef( 0 ),
          uiBeginOnQuery( 0 ),
          uiEndOnQuery( 0 ),
          xStats( ),
          bSecondary( false ),
          bSupplementary( false )
    {} // constructor

    /**
     * @brief Creates an empty alignment,
     * where the interval of the reference that is used is already known.
     */
    Alignment( nucSeqIndex uiBeginOnRef )
        : data( ),
          uiLength( 0 ),
          uiBeginOnRef( uiBeginOnRef ),
          uiEndOnRef( uiBeginOnRef ),
          uiBeginOnQuery( 0 ),
          uiEndOnQuery( 0 ),
          xStats( ),
          bSecondary( false ),
          bSupplementary( false )
    {} // constructor

    /**
     * @brief Creates an empty alignment,
     * where the interval of the reference that is used is already known.
     */
    Alignment( nucSeqIndex uiBeginOnRef, nucSeqIndex uiBeginOnQuery )
        : data( ),
          uiLength( 0 ),
          uiBeginOnRef( uiBeginOnRef ),
          uiEndOnRef( uiBeginOnRef ),
          uiBeginOnQuery( uiBeginOnQuery ),
          uiEndOnQuery( uiBeginOnQuery ),
          xStats( ),
          bSecondary( false ),
          bSupplementary( false )
    {} // constructor

    /**
     * @brief Creates an empty alignment,
     * where the interval of the reference that is used is already known.
     */
    Alignment( nucSeqIndex uiBeginOnRef, nucSeqIndex uiBeginOnQuery, nucSeqIndex uiEndOnRef, nucSeqIndex uiEndOnQuery )
        : data( ),
          uiLength( 0 ),
          uiBeginOnRef( uiBeginOnRef ),
          uiEndOnRef( uiEndOnRef ),
          uiBeginOnQuery( uiBeginOnQuery ),
          uiEndOnQuery( uiEndOnQuery ),
          xStats( ),
          bSecondary( false ),
          bSupplementary( false )
    {} // constructor

    Alignment( const Alignment& rOther ) = delete;

    inline std::string toString( ) const
    {
        std::string sRet = "Alignment Dump: ";
        constexpr char vTranslate[ 5 ] = {'S', '=', 'X', 'I', 'D'};
        for( auto& tuple : data )
            sRet += vTranslate[ (unsigned int)tuple.first ] + std::to_string( tuple.second ) + " ";
        sRet += " | Score: ";
        sRet += std::to_string( iScore );
        return sRet;
    }

    inline std::shared_ptr<Seeds> toSeeds( ) const
    {
        auto pRet = std::make_shared<Seeds>( );
        nucSeqIndex uiQ = uiBeginOnQuery;
        nucSeqIndex uiR = uiBeginOnRef;
        for( auto& rPair : data )
        {
            switch( rPair.first )
            {
                case MatchType::seed:
                case MatchType::match:
                    // @todo make this work for rev strand alignments...
                    pRet->emplace_back( uiQ, rPair.second, uiR, true );
                case MatchType::missmatch: // seed and match case continue till the break
                    uiQ += rPair.second;
                    uiR += rPair.second;
                    break;
                case MatchType::insertion:
                    uiQ += rPair.second;
                    break;
                case MatchType::deletion:
                    uiR += rPair.second;
                    break;
                default:
                    throw std::runtime_error( "Should never reach default case in toSeeds switch case." );
                    break;
            } // switch
        } // for
        return pRet;
    } // method

    // overload
    bool canCast( std::shared_ptr<libMS::Container> c ) const
    {
        return std::dynamic_pointer_cast<Alignment>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "Alignment";
    } // function

    // overload
    std::shared_ptr<libMS::Container> getType( ) const
    {
        return std::shared_ptr<libMS::Container>( new Alignment( ) );
    } // function

    int64_t DLL_PORT( MA ) reCalcScore( ) const;

    /**
     * @returns the type of math for the given position i.
     * @brief Type of math at i.
     */
    MatchType at( nucSeqIndex i ) const
    {
        // everything after and before the query is a deletion
        if( i >= uiLength || i < 0 )
            return MatchType::deletion;

        // the MatchType match type is stored in a compressed format -> extract it
        nucSeqIndex j = 0;
        unsigned int k = 0;
        while( k < data.size( ) )
        {
            j += data[ k ].second;
            if( j > i )
                break;
            k++;
        } // while

        return data[ k ].first;
    } // function

    /*
     * @brief returns how many seeds are within the alignment
     */
    size_t getNumSeeds( ) const
    {
        size_t uiRet = 0;
        for( std::pair<MatchType, nucSeqIndex> section : data )
            if( section.first == MatchType::seed )
                uiRet++;
        return uiRet;
    } // function

    /*
     * @brief returns how many nucleotides within this alignment are determined by seeds
     */
    size_t getNumBySeeds( ) const
    {
        size_t uiRet = 0;
        for( std::pair<MatchType, nucSeqIndex> section : data )
            if( section.first == MatchType::seed )
                uiRet += section.second;
        return uiRet;
    } // function

    /*
     * @brief returns how many nucleotides within this alignment matches (or seeds)
     */
    size_t getNumMatches( ) const
    {
        size_t uiRet = 0;
        for( std::pair<MatchType, nucSeqIndex> section : data )
            if( section.first == MatchType::seed || section.first == MatchType::match )
                uiRet += section.second;
        return uiRet;
    } // function

    /*
     * Essentially computes the NM tag of the SAM format
     * NM:i:count Number of differences (mismatches plus inserted and deleted bases) between the sequence and
     * reference, counting only (case-insensitive) A, C, G and T bases in sequence and reference as potential
     * matches, with everything else being a mismatch. Note this means that ambiguity codes in both
     * sequence and reference that match each other, such as ‘N’ in both, or compatible codes such as ‘A’ and
     * ‘R’, are still counted as mismatches. The special sequence base ‘=’ will always be considered to be a
     * match, even if the reference is ambiguous at that point. Alignment reference skips, padding, soft and
     * hard clipping (‘N’, ‘P’, ‘S’ and ‘H’ CIGAR operations) do not count as mismatches, but insertions and
     * deletions count as one mismatch per base.
     * Note that historically this has been ill-defined and both data and tools exist that disagree with this
     * definition.
     *
     * @note: According to NGMLR documentation, they do not count insertions and deletions
     */
    size_t getNumDifferences( const std::shared_ptr<Pack> pPack, const bool bNMTagDoNOTCountIndels ) const
    {
        NucSeq xRefSeqWithN;
        pPack->vExtractSubsectionN( this->uiBeginOnRef, this->uiEndOnRef, xRefSeqWithN );
        size_t uiNumDiff = 0;
        size_t uiRPos = 0;
        for( const std::pair<MatchType, nucSeqIndex>& rPair : this->data )
        {
            switch( rPair.first )
            {
                case MatchType::seed:
                case MatchType::match:
                    for( size_t uiI = 0; uiI < rPair.second; uiI++ )
                        if( xRefSeqWithN[ uiI + uiRPos ] >= 4 ) // there is a N at this position
                            uiNumDiff++;
                    break;
                case MatchType::insertion:
                case MatchType::deletion:
                    // if we do not want to count indels, bNMTagDoNOTCountIndels must be set to true.
                    if( bNMTagDoNOTCountIndels )
                        break;
                case MatchType::missmatch:
                    uiNumDiff += rPair.second;
                    break;
                default:
                    throw std::runtime_error( "Should never reach default case in computeTag switch case." );
                    break;
            } // switch
            if( rPair.first != MatchType::insertion )
                uiRPos += rPair.second;
        } // for
        return uiNumDiff;
    } // method


    /*
     * sniffles requirement:
     * if we fill a dual extend gap by brute force (due to z drop) the used insertion and deletion must always be in
     * consistent order (among several reads).
     * Therefore we need to check the strand of the alignment and adjust accordingly.
     */
    void invertSuccessiveInserionAndDeletion( )
    {
        // using -1 in loop therefore start at 1 instead of 0
        for( size_t uiI = 1; uiI < this->data.size( ); uiI++ )
        {
            if( ( this->data[ uiI - 1 ].first == MatchType::insertion &&
                  this->data[ uiI ].first == MatchType::deletion ) ||
                ( this->data[ uiI - 1 ].first == MatchType::deletion &&
                  this->data[ uiI ].first == MatchType::insertion ) )
            {
                auto xTmp = this->data[ uiI ];
                this->data[ uiI ] = this->data[ uiI - 1 ];
                this->data[ uiI - 1 ] = xTmp;
                uiI++; // increment uiI twice to skip one (so that we do not swap back and forth)
            } // if
        } // for
    } // method

    /**
     * @returns the type of math for the given position i.
     * @brief Type of math at i.
     */
    MatchType operator[]( nucSeqIndex i ) const
    {
        return at( i );
    } // operator

    /**
     * @brief extract the alignment as vector
     */
    std::vector<MatchType> extract( ) const
    {
        std::vector<MatchType> aRet;
        for( std::pair<MatchType, nucSeqIndex> xElement : data )
            for( unsigned int i = 0; i < xElement.second; i++ )
                aRet.push_back( xElement.first );
        return aRet;
    } // method

    std::string cigarString( Pack& rPack )
    {
        std::string sCigar = "";
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            std::reverse( data.begin( ), data.end( ) );

        for( std::pair<MatchType, nucSeqIndex> section : data )
        {
            sCigar.append( std::to_string( section.second ) );
            switch( section.first )
            {
                case MatchType::seed:
                case MatchType::match:
                    sCigar.append( "=" );
                    break;
                case MatchType::missmatch:
                    sCigar.append( "X" );
                    break;
                case MatchType::insertion:
                    sCigar.append( "I" );
                    break;
                case MatchType::deletion:
                    sCigar.append( "D" );
                    break;
                default:
                    std::cerr << "WARNING invalid cigar symbol" << std::endl;
                    break;
            } // switch
        } // for
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            std::reverse( data.begin( ), data.end( ) );
        return sCigar;
    } // method

    std::string cigarStringWithMInsteadOfXandEqual( Pack& rPack )
    {
        std::string sCigar = "";
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            std::reverse( data.begin( ), data.end( ) );

        size_t uiSequentialM = 0;
        for( std::pair<MatchType, nucSeqIndex> section : data )
        {
            switch( section.first )
            {
                case MatchType::seed:
                case MatchType::match:
                case MatchType::missmatch:
                    // add up sequential M's
                    uiSequentialM += section.second;
                    break;
                case MatchType::insertion:
                case MatchType::deletion:
                    // output added up sequential M's
                    if( uiSequentialM > 0 )
                    {
                        sCigar.append( std::to_string( uiSequentialM ) ).append( "M" );
                        uiSequentialM = 0;
                    } // if
                    // output indel
                    sCigar.append( std::to_string( section.second ) )
                        .append( ( section.first == MatchType::insertion ) ? "I" : "D" );
                    break;
                default:
                    std::cerr << "WARNING invalid cigar symbol" << std::endl;
                    break;
            } // switch
        } // for

        // don't forget potential last M!
        if( uiSequentialM > 0 )
            sCigar.append( std::to_string( uiSequentialM ) ).append( "M" );

        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            std::reverse( data.begin( ), data.end( ) );
        return sCigar;
    } // method

    void appendCigarString( std::string sCigar, std::shared_ptr<NucSeq> pQuery, Pack& rPack )
    {
        std::stringstream xStream( sCigar );
        while( !xStream.eof( ) )
        {
            size_t uiLen;
            char cOp;
            xStream >> uiLen;
            xStream >> cOp;
            switch( cOp )
            {
                case 'M':
                    if( pQuery->pxSequenceRef[ uiEndOnQuery ] == rPack.vExtract( uiEndOnRef ) )
                        append( MatchType::match, uiLen );
                    else
                        append( MatchType::missmatch, uiLen );
                    break;
                case 'I':
                    append( MatchType::insertion, uiLen );
                    break;
                case 'D':
                    append( MatchType::deletion, uiLen );
                    break;
                case 'N':
                    append( MatchType::deletion, uiLen );
                    break;
                case 'S':
                    append( MatchType::insertion, uiLen );
                    break;
                case 'H':
                case 'P':
                    /*clipped sequences NOT present in SEQ*/
                    break;
                case '=':
                    append( MatchType::match, uiLen );
                    break;
                case 'X':
                    append( MatchType::missmatch, uiLen );
                    break;

                default:
                    throw std::runtime_error( "Unrecognized SAM CIGAR operation" );
                    break;
            } // switch
        } // while
    } // method

    uint32_t getSamFlag( Pack& rPack ) const
    {
        uint32_t uiRet = 0;
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            uiRet |= REVERSE_COMPLEMENTED;
        if( bSecondary )
            uiRet |= SECONDARY_ALIGNMENT;
        if( bSupplementary )
            uiRet |= SUPPLEMENTARY_ALIGNMENT;
        return uiRet;
    } // method

    std::string getContig( Pack& rPack ) const
    {
        return rPack.nameOfSequenceForPosition( uiBeginOnRef );
    } // method

    nucSeqIndex getSamPosition( Pack& rPack ) const
    {
        std::string sRefName = getContig( rPack );
        auto uiRet = rPack.posInSequence( uiBeginOnRef, uiEndOnRef );
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            uiRet += 1;
        return uiRet + 1;
    } // method

    std::string getQuerySequence( NucSeq& rQuery, Pack& rPack )
    {
        std::string sRet;
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            sRet = rQuery.fromToComplement( uiBeginOnQuery, uiEndOnQuery );
        else
            sRet = rQuery.fromTo( uiBeginOnQuery, uiEndOnQuery );
        int64_t iOff = sRet.length( ) - static_cast<int64_t>( this->uiEndOnQuery - this->uiBeginOnQuery );
        if( iOff != 0 )
            throw std::runtime_error(
                "Query length is off by " + std::to_string( iOff ) + ". Query is on " +
                ( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) ? "refstrand" : "forwstrand" ) +
                ". Alignment name:" + xStats.sName + " " + this->toString( ) + ". Expected length: " +
                std::to_string( static_cast<int64_t>( this->uiEndOnQuery - this->uiBeginOnQuery ) ) +
                ". Querysegment: " + sRet + ". Query size: " + std::to_string( rQuery.uiSize ) +
                ". Alignment pos on query: " + std::to_string( this->uiBeginOnQuery ) + " to " +
                std::to_string( this->uiEndOnQuery ) );
        return sRet;
    } // method

    std::string getQueryQuality( NucSeq& rQuery )
    {
        return rQuery.fromToQual( uiBeginOnQuery, uiEndOnQuery );
    } // method

    std::string getRevCompQuerySequence( NucSeq& rQuery, Pack& rPack )
    {
        if( rPack.bPositionIsOnReversStrand( uiBeginOnRef ) )
            return rQuery.fromTo( uiBeginOnQuery, uiEndOnQuery );
        else
            return rQuery.fromToComplement( uiBeginOnQuery, uiEndOnQuery );
    } // method

    /**
     * @brief appends multiple matchTypes to the alignment
     * @details
     * This is used for appending seeds,
     * where simply size of the seed matches need to be appended.
     */
    void DLL_PORT( MA ) append( MatchType type, nucSeqIndex size );

    /**
     * @brief appends a matchType to the alignment
     */
    void append( MatchType type )
    {
        append( type, 1 );
    } // function

    /**
     * @brief returns the overlap of the alignments on the query
     * @details
     * does consider deletions correctly. For example:
     * alignment 1: ###--#---####
     * alignment 2: ---###-------
     * 1 and 2 would have overlap 1/13
     */
    inline double overlap( const Alignment& rOther ) const
    {
        // get the total area where overlaps are possible
        nucSeqIndex uiS = std::max( uiBeginOnQuery, rOther.uiBeginOnQuery );
        nucSeqIndex uiE = std::min( uiEndOnQuery, rOther.uiEndOnQuery );
        // if the total area is zero we can return 0% immediately
        if( uiS >= uiE )
            return 0;
        // indices for walking through the alignments
        nucSeqIndex uiOverlap = 0;
        size_t uiI = 0;
        nucSeqIndex uiQpos = uiBeginOnQuery;
        size_t uiIOther = 0;
        nucSeqIndex uiQposOther = rOther.uiBeginOnQuery;
        // move the index to the start positions
        while( uiQpos + data[ uiI ].second < uiS )
        {
            // all match types move forward on the query apart from a deletion
            if( data[ uiI ].first != MatchType::deletion )
                uiQpos += data[ uiI ].second;
            uiI++;
        } // while
        // move the index (of the other alignment) to the start positions
        while( uiQposOther + rOther.data[ uiIOther ].second < uiS )
        {
            // all match types move forward on the query apart from a deletion
            if( rOther.data[ uiIOther ].first != MatchType::deletion )
                uiQposOther += rOther.data[ uiIOther ].second;
            uiIOther++;
        } // while

        // compute the overlap
        while(
            // we haven't reached the end of the possible overlap area
            uiQpos < uiE && uiQposOther < uiE &&
            // we haven't reached the end of the alignment
            uiI < data.size( ) && uiIOther < rOther.data.size( ) )
        {
            // get the lengths on the query
            nucSeqIndex uiQLen = 0;
            // all match types have a query length apart from a deletion
            if( data[ uiI ].first != MatchType::deletion )
                uiQLen = data[ uiI ].second;
            nucSeqIndex uiQLenOther = 0;
            // all match types have a query length apart from a deletion
            if( rOther.data[ uiIOther ].first != MatchType::deletion )
                uiQLenOther = rOther.data[ uiIOther ].second;

            // compute the overlap
            nucSeqIndex uiS_inner = std::max( std::max( uiQpos, uiQposOther ), uiS );
            nucSeqIndex uiE_inner = std::min( std::min( uiQpos + uiQLen, uiQposOther + uiQLenOther ), uiE );
            nucSeqIndex uiCurrOverlap = 0;
            if( uiS_inner < uiE_inner )
                uiCurrOverlap = uiE_inner - uiS_inner;

            // if the type of the match is not an insertion, add to the amount of overlap
            if( data[ uiI ].first != MatchType::insertion && rOther.data[ uiIOther ].first != MatchType::insertion )
                uiOverlap += uiCurrOverlap;

            // move the correct index forward
            if( uiQpos + uiQLen < uiQposOther + uiQLenOther )
            {
                uiQpos += uiQLen;
                uiI++;
            } // if
            else
            {
                uiQposOther += uiQLenOther;
                uiIOther++;
            } // else
        } // while

        // divide by the size of the smaller alignment so that the returned overlap is in
        // percent
        nucSeqIndex uiSize = std::min( uiEndOnQuery - uiBeginOnQuery, rOther.uiEndOnQuery - rOther.uiBeginOnQuery );
        return uiOverlap / static_cast<double>( uiSize );
    } // mehtod

    /**
     * @brief appends another alignment
     */
    void append( const Alignment& rOther )
    {
        for( auto xTuple : rOther.data )
            append( xTuple.first, xTuple.second );
    } // function

    ///@brief wrapper for boost-python
    void append_boost1( MatchType type, nucSeqIndex size )
    {
        append( type, size );
    } // function

    ///@brief wrapper for boost-python
    void append_boost2( MatchType type )
    {
        append( type, 1 );
    } // function

    ///@returns the length of the alignment
    ///@brief Length of the alignment
    nucSeqIndex length( ) const
    {
        DEBUG( nucSeqIndex uiCheck = 0; for( auto xTup
                                             : data ) uiCheck += xTup.second;
               if( uiCheck != uiLength ) {
                   std::cout << "Alignment length check failed: " << uiCheck << " != " << uiLength << std::endl;
                   assert( false );
               } // if
               ) // DEBUG
        return uiLength;
    } // function

    ///@brief Start of the alignment.
    ///@returns the start of the alignment on the reference sequence.
    nucSeqIndex beginOnRef( )
    {
        return uiBeginOnRef;
    } // function

    ///@brief End of the alignment.
    ///@returns the end of the alignment on the reference sequence.
    nucSeqIndex endOnRef( )
    {
        return uiEndOnRef;
    } // function

    /**
     * @brief the NMW score for this alignment
     */
    int64_t score( ) const
    {
        // the data.size() == 0 is to allow SW to set the score directly
        assert( data.size( ) == 0 || reCalcScore( ) == iScore );
        return iScore;
    } // method

    /**
     * @brief the NMW score for this alignment
     */
    unsigned int DLL_PORT( MA ) localscore( ) const;

    /**
     * @brief returns how many nucleotides within this alignment are determined by seeds
     * as a percentage
     */
    float seedCoverage( ) const
    {
        unsigned int iCount = 0;
        for( unsigned int i = 0; i < length( ); i++ )
            if( at( i ) == MatchType::seed )
                iCount++;
        return ( (float)iCount ) / (float)length( );
    } // method

    /**
     * @brief for sorting alignment by their score
     * @details
     * When multiple alignments are created we use this function to sort them.
     */
    bool larger( const std::shared_ptr<libMS::Container> pOther ) const
    {
        const std::shared_ptr<Alignment> pAlign = std::dynamic_pointer_cast<Alignment>( pOther );
        if( pAlign == nullptr )
            return false;
        size_t uiA = 0;
        size_t uiB = 0;
        if( pAlign->bSecondary )
            uiB = 2;
        if( pAlign->bSupplementary )
            uiB = 1;
        if( bSecondary )
            uiA = 2;
        if( bSupplementary )
            uiA = 1;
        if( uiA != uiB )
            return uiA < uiB;
        auto uiS1 = score( );
        auto uiS2 = pAlign->score( );
        if( uiS1 == uiS2 )
            // if both alignments have the same score output the one with
            // the higher SoC score first (this is determined by the lower SoC index)
            return xStats.index_of_strip < pAlign->xStats.index_of_strip;
        return uiS1 > uiS2;
    } // function

    /**
     * @brief transform any alignment into a local one
     * @details
     * When an alignment is computed on the foundation of seeds it might not be local.
     * This function has a linear complexity with regard to the compressed alignment length.
     */
    void DLL_PORT( MA ) makeLocal( );

    /**
     * @brief removes dangeling Deletions
     * @details
     * When the alignment is created there might be some dangeling deletions at
     * the beginning or end. This function removes them
     */
    void DLL_PORT( MA ) removeDangeling( );

    void operator=( const std::shared_ptr<Alignment> pOther )
    {
        data = pOther->data;
        uiLength = pOther->uiLength;
        uiBeginOnRef = pOther->uiBeginOnRef;
        uiEndOnRef = pOther->uiEndOnRef;
        uiBeginOnQuery = pOther->uiBeginOnQuery;
        iScore = pOther->iScore;
        fMappingQuality = pOther->fMappingQuality;
        bSecondary = pOther->bSecondary;
        bSupplementary = pOther->bSupplementary;
        xStats = pOther->xStats;
    } // function

    inline void shiftOnRef( nucSeqIndex uiBy )
    {
        uiBeginOnRef += uiBy;
        uiEndOnRef += uiBy;
    } // method

    inline void shiftOnQuery( nucSeqIndex uiBy )
    {
        uiBeginOnQuery += uiBy;
        uiEndOnQuery += uiBy;
    } // method

    /**
     * @brief debug function
     * @details
     * Checks wether the cigar and uiBegin and EndOnQuery are consistent.
     * Does nothing in release mode.
     */
    void checkLengthOnQuery( )
    {
        assert( this->uiBeginOnQuery <= this->uiEndOnQuery );
#if DEBUG_LEVEL > 0
        int64_t uiTlen = static_cast<int64_t>( this->uiEndOnQuery ) - static_cast<int64_t>( this->uiBeginOnQuery );
        for( std::pair<MatchType, nucSeqIndex> section : this->data )
        {
            switch( section.first )
            {
                case MatchType::seed:
                case MatchType::match:
                case MatchType::insertion:
                case MatchType::missmatch:
                    uiTlen -= static_cast<int64_t>( section.second );
                    break;
                case MatchType::deletion:
                default:
                    break;
            } // switch
        } // for
        if( uiTlen != 0 )
        {
            std::cerr << "Error: Alignment length on query does not match cigar. Length is off by " << uiTlen
                      << std::endl;
            std::cerr << "Query name: " << xStats.sName << std::endl;
            std::cerr << this->toString( ) << std::endl;
            std::cerr << "Query Coords: " << this->uiBeginOnQuery << " to " << this->uiEndOnQuery << std::endl;
            assert( false );
        } // if
#endif
    } // method
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportAlignment( );
#else
void exportAlignment( libMS::SubmoduleOrganizer& xOrganizer );
#endif
#endif

#endif