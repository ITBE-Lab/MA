#include "ma/container/alignment.h"
#include "ma/container/pack.h"
#include "ma/module/fileReader.h"
#include <map>

#pragma once

namespace libMA
{

class ReadByName : public libMS::Container
{
    std::map<std::string, std::shared_ptr<NucSeq>> xQueries;

  public:
    ReadByName( std::vector<std::shared_ptr<NucSeq>> vQueries )
    {
        for( auto pQuery : vQueries )
            append( pQuery );
    } // constructor

    ReadByName( )
    {}

    void append( std::shared_ptr<NucSeq> pQuery )
    {
        if( pQuery->sName.empty( ) )
            throw std::runtime_error( "unnamed query" );
        if( xQueries.count( pQuery->sName ) > 0 )
        {
            if( !xQueries[ pQuery->sName ]->isEqual( *pQuery ) )
                throw std::runtime_error( "query with that name but different sequence already exists" );
        } // if
        else
            xQueries[ pQuery->sName ] = pQuery;
    } // method

    std::shared_ptr<NucSeq> operator[]( std::string sName )
    {
        if( xQueries.count( sName ) == 0 )
            throw std::runtime_error( "unknown query name" );
        return xQueries[ sName ];
    } // method

#ifdef WITH_PYTHON
    py::iterator iter( )
    {
        return py::make_iterator( xQueries.begin( ), xQueries.end( ) );
    }
#endif

}; // class ReadByName

class SamFileReader : public libMS::Module<Alignment, true, FileStream, Pack, ReadByName>
{
  public:
    /**
     * @brief creates a new SamFileReader.
     */
    SamFileReader( const ParameterSetManager& rParameters )
    {} // constructor

    /// http://www.martinbroadhurst.com/how-to-split-a-string-in-c.html
    template <class Container_t> Container_t splitString( const std::string& sString, char sDelim = ' ' )
    {
        std::stringstream xStream( sString );
        std::string sCurrToken;
        Container_t xCont;
        while( std::getline( xStream, sCurrToken, sDelim ) )
            xCont.push_back( sCurrToken );
        return xCont;
    } // method

    std::shared_ptr<Alignment> DLL_PORT( MA )
        execute( std::shared_ptr<FileStream> pStream, std::shared_ptr<Pack> pRef, std::shared_ptr<ReadByName> pReads )
    {
        std::string sLine = "";
        {
            std::lock_guard<std::mutex> xLock( pStream->xMutex );
            pStream->peek( ); // potentionally trigger eof
            if( pStream->eof( ) ) // eof case
                return nullptr;
            // ignore empty lines and comment/header lines (starting with '@')
            while( sLine.empty( ) || sLine[ 0 ] == '@' )
                pStream->safeGetLine( sLine );
            pStream->peek( ); // potentionally trigger eof
        } // scope for xLock
        auto vColumns = splitString<std::vector<std::string>>( sLine, '\t' );
        if( vColumns.size( ) <= 5 )
            throw std::runtime_error( "too little tab seperated columns for a SAM file!" );
        std::shared_ptr<NucSeq> pQuery;
        // only use pReads is sequence is not given in sam file
        if( vColumns[ 9 ] != "*" )
            pQuery = ( *pReads )[ vColumns[ 0 ] ];
        else
            pQuery = std::make_shared<NucSeq>( vColumns[ 9 ] );
        // vColumns[ 3 ] == POS in sam is 1 based
        nucSeqIndex uiRefStart = atoll( vColumns[ 3 ].c_str( ) ) + pRef->startOfSequenceWithName( vColumns[ 2 ] ) - 1;
        auto pRet = std::make_shared<Alignment>( uiRefStart );
        pRet->xStats.sName = vColumns[ 0 ];
        pRet->appendCigarString( vColumns[ 5 ], pQuery, *pRef );
        return pRet;
    } // method

}; // class SamFileReader

class SeedsByName : public libMS::Container
{
    std::map<std::string, std::shared_ptr<Seeds>> xSeeds;

  public:
    SeedsByName( )
    {} // constructor

    void append( std::shared_ptr<Seeds> pSeeds, std::string sName )
    {
        xSeeds[ sName ] = pSeeds;
    } // method

    std::shared_ptr<Seeds> operator[]( std::string sName )
    {
        if( xSeeds.count( sName ) == 0 )
            throw std::runtime_error( "unknown seed set name" );
        return xSeeds[ sName ];
    } // method
}; // class SeedsByName


class GetSeedsByName : public libMS::Module<Seeds, false, Alignment, SeedsByName>
{
  public:
    /**
     * @brief creates a new GetSeedsByName module.
     */
    GetSeedsByName( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<Alignment> pAlignment, std::shared_ptr<SeedsByName> pSeedsByName )
    {
        return ( *pSeedsByName )[ pAlignment->xStats.sName ];
    } // method

}; // class GetSeedsByName

class GetSeedsByReadName : public libMS::Module<Seeds, false, NucSeq, SeedsByName>
{
  public:
    /**
     * @brief creates a new GetSeedsByName module.
     */
    GetSeedsByReadName( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<Seeds> DLL_PORT( MA )
        execute( std::shared_ptr<NucSeq> pRead, std::shared_ptr<SeedsByName> pSeedsByName )
    {
        return ( *pSeedsByName )[ pRead->sName ];
    } // method

}; // class GetSeedsByName

class GetReadByName : public libMS::Module<NucSeq, false, Alignment, ReadByName>
{
  public:
    /**
     * @brief creates a new GetReadByName module.
     */
    GetReadByName( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<NucSeq> DLL_PORT( MA )
        execute( std::shared_ptr<Alignment> pAlignment, std::shared_ptr<ReadByName> pReadByName )
    {
        return ( *pReadByName )[ pAlignment->xStats.sName ];
    } // method

}; // class GetReadByName

} // namespace libMA


#ifdef WITH_PYTHON
void exportSamFileReader( libMS::SubmoduleOrganizer& xOrganizer );
#endif