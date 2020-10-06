/**
 * @file read.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "ma/container/alignment.h"
#include "msv/container/svJump.h"
#include "sql_api.h" // NEW DATABASE INTERFACE

namespace libMSV
{


template <typename DBCon>
using ReadTableType = SQLTableWithLibIncrPriKey<DBCon,
                                                PriKeyDefaultType, // sequencer id (foreign key)
                                                std::string, // read name
                                                std::shared_ptr<CompressedNucSeq> // read sequence
                                                >;
const json jReadTableDef = {
    { TABLE_NAME, "read_table" },
    { TABLE_COLUMNS,
      { { { COLUMN_NAME, "sequencer_id" } }, { { COLUMN_NAME, "_name_" } }, { { COLUMN_NAME, "sequence" } } } },
    { FOREIGN_KEY, { { COLUMN_NAME, "sequencer_id" }, { REFERENCES, "sequencer_table(id)" } } } };
/**
 * @brief this table saves reads
 */
template <typename DBCon> class ReadTable : public ReadTableType<DBCon>
{

  public:
    // std::shared_ptr<DBCon> pDB;
    SQLQuery<DBCon, PriKeyDefaultType> xGetReadId;
    SQLQuery<DBCon, uint64_t> xGetNumMatchingReads;
    SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, std::string> xGetRead;
    SQLQuery<DBCon, PriKeyDefaultType> xGetSeqId;

    ReadTable( std::shared_ptr<DBCon> pDB )
        : ReadTableType<DBCon>( pDB, jReadTableDef ),

          // pDB( pDB ),
          xGetReadId( pDB, "SELECT id FROM read_table WHERE sequencer_id = ? AND _name_ = ? " ),
          xGetNumMatchingReads( pDB, "SELECT COUNT(*) FROM read_table WHERE sequencer_id = ? AND _name_ = ? " ),
          xGetRead( pDB, "SELECT sequence, _name_ FROM read_table WHERE id = ? " ),
          xGetSeqId( pDB, "SELECT sequencer_id FROM read_table WHERE id = ? " )
    {} // default constructor


    inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
    {
        return ReadTableType<DBCon>::insert( uiSequencerId, pRead->sName, makeSharedCompNucSeq( *pRead ) );
    } // method

    inline std::shared_ptr<NucSeq> getRead( int64_t iId )
    {
        if( !xGetRead.execAndFetch( (PriKeyDefaultType)iId ) )
            throw std::runtime_error( "Read with id " + std::to_string( iId ) +
                                      " could not be found in the database." );
        auto xTuple = xGetRead.get( );
        std::get<0>( xTuple )->pUncomNucSeq->iId = iId;
        std::get<0>( xTuple )->pUncomNucSeq->sName = std::get<1>( xTuple );
        if( xGetRead.next( ) )
            assert( false ); // can never find two reads with same ID
        return std::get<0>( xTuple )->pUncomNucSeq;
    } // method
    inline int64_t getSeqId( int64_t iReadId )
    {
        return xGetSeqId.scalar( iReadId );
    } // method

    inline int64_t getReadId( int64_t iSeqId, std::string sName )
    {
        // char* pcEscapedName = PQescapeLiteral( this->pDB->pPGConn, sName.c_str( ), sName.size( ) );
        // std::string sEscapedName( pcEscapedName );
        if( xGetNumMatchingReads.scalar( (PriKeyDefaultType)iSeqId, sName ) == 0 )
            return -1;
        return xGetReadId.scalar( (PriKeyDefaultType)iSeqId, sName );
    } // method

    inline std::vector<std::shared_ptr<NucSeq>> getUsedReads( std::shared_ptr<DBCon> pDB )
    {
        SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, std::string> xGetAllUsedReads(
            pDB, "SELECT sequence, _name_ FROM read_table WHERE id IN (SELECT DISTINCT read_id FROM sv_jump_table)" );
        std::vector<std::shared_ptr<NucSeq>> vRet;
        xGetAllUsedReads.execAndForAll( [ & ]( std::shared_ptr<CompressedNucSeq> pComp, std::string sName ) {
            vRet.push_back( pComp->pUncomNucSeq );
            vRet.back( )->sName = sName;
        } );
        return vRet;
    } // method
}; // class

} // namespace libMSV

struct RangeStartInt64_t
{
    int64_t iP;
};
struct RangeEndInt64_t
{
    int64_t iP;
};

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<RangeStartInt64_t>( )
{
    return "int8range";
} // specialized method
// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<RangeEndInt64_t>( )
{
    // needs to be this type eventhough it is not used in create
    // we use it to check/set the oid of outputs & inputs...
    return "int8range";
} // specialized method

// Part1 : Specify that the end type is not used for create statements
template <> inline bool PostgreSQLDBCon::TypeTranslator::useInCreate<RangeEndInt64_t>( )
{
    return false;
} // specialized method

// Part1b : Spatial types require an indication that the argument passed at a placeholder's
//          position has the format 'WKB'.
template <>
inline std::string
PostgreSQLDBCon::TypeTranslator::getPlaceholderForType<RangeStartInt64_t>( const std::string& rsInsertedText )
{
    // closing bracket is im RangeEndInt64_t
    return "int8range(" + rsInsertedText;
} // specialized method

// Part1b : Spatial types require an indication that the argument passed at a placeholder's
//          position has the format 'WKB'.
template <>
inline std::string
PostgreSQLDBCon::TypeTranslator::getPlaceholderForType<RangeEndInt64_t>( const std::string& rsInsertedText )
{
    // opening bracket is im RangeStartInt64_t
    // comma is already provided
    return rsInsertedText + ")";
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void PostgreSQLDBCon::StmtArg::set( const RangeStartInt64_t& rRange )
{
    // forwards to the overload for int64_t
    PostgreSQLDBCon::StmtArg::set( rRange.iP );
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void PostgreSQLDBCon::StmtArg::set( const RangeEndInt64_t& rRange )
{
    // forwards to the overload for int64_t
    PostgreSQLDBCon::StmtArg::set( rRange.iP );
} // specialized method

// Part 3: Code for supporting query output: Here we can just forward to PGRowCell<int64_t>
template <>
struct /* PostgreSQLDBCon:: */ PGRowCell<RangeStartInt64_t> : public /* PostgreSQLDBCon::*/ PGRowCell<int64_t>
{}; // specialized class

// Part 3: Code for supporting query output: Here we can just forward to PGRowCell<int64_t>
template <> struct /* PostgreSQLDBCon:: */ PGRowCell<RangeEndInt64_t> : public /* PostgreSQLDBCon::*/ PGRowCell<int64_t>
{}; // specialized class

namespace libMSV
{

template <typename DBCon>
using ReadRangeTableType = SQLTableWithLibIncrPriKey<DBCon,
                                                     PriKeyDefaultType, // read id (foreign key)
                                                     RangeStartInt64_t, // read_range
                                                     RangeEndInt64_t, // read_range (dummy)
                                                     bool, // primary alignment
                                                     double // mapping quality
                                                     >;
const json jReadExtensionTableDef = {
    { TABLE_NAME, "read_range_table" },
    { TABLE_COLUMNS,
      { { { COLUMN_NAME, "read_id" } },
        { { COLUMN_NAME, "read_range" } },
        { { COLUMN_NAME, "read_range_DUMMY" } },
        { { COLUMN_NAME, "primary_alignment" } },
        { { COLUMN_NAME, "mapping_quality" } } } },
    { FOREIGN_KEY, { { COLUMN_NAME, "read_id" }, { REFERENCES, "read_table(id)" } } } };
/**
 * @brief this table saves reads
 */
template <typename DBCon> class ReadRangeTable : public ReadRangeTableType<DBCon>
{
    SQLQuery<DBCon, uint32_t> xGetCoverage;
    SQLQuery<DBCon, uint32_t> xGetPrimCoverage;

  public:
    ReadRangeTable( std::shared_ptr<DBCon> pDB )
        : ReadRangeTableType<DBCon>( pDB, jReadExtensionTableDef ),
          xGetCoverage( pDB, "SELECT COUNT(*) "
                             "FROM read_range_table "
                             "JOIN read_table ON read_table.id = read_range_table.read_id "
                             "WHERE int8range(?,?) <@ read_range " // checks if alignment encloses given range
                             "AND sequencer_id = ? "
                             "AND mapping_quality >= ? " ),
          xGetPrimCoverage( pDB, "SELECT COUNT(*) "
                                 "FROM read_range_table "
                                 "JOIN read_table ON read_table.id = read_range_table.read_id "
                                 "WHERE int8range(?,?) <@ read_range " // checks if alignment encloses given range
                                 "AND sequencer_id = ? "
                                 "AND mapping_quality >= ? "
                                 "AND primary_alignment " )
    {} // default constructor

    inline void insertAlignmentId( int64_t iReadId, std::shared_ptr<Alignment> pAlignment )
    {
        ReadRangeTable<DBCon>::insert( iReadId, RangeStartInt64_t{ (int64_t)pAlignment->beginOnRef( ) },
                                       RangeEndInt64_t{ (int64_t)pAlignment->endOnRef( ) },
                                       !pAlignment->bSecondary && !pAlignment->bSupplementary,
                                       (double)pAlignment->fMappingQuality );
    } // method

    inline void insertAlignment( std::shared_ptr<NucSeq> pRead, std::shared_ptr<Alignment> pAlignment )
    {
        insertAlignmentId( pRead->iId, pAlignment );
    } // method

    inline void genIndices( )
    {
        this->addIndex(
            json{ { INDEX_NAME, "range_index" }, { INDEX_COLUMNS, "read_range" }, { INDEX_METHOD, "GIST" } } );
    } // method

    inline void dropIndices( )
    {
        this->dropIndex( json{ { INDEX_NAME, "read_range" } } );
    } // method

    inline uint32_t coverage( int64_t from, int64_t to, int64_t iSeqId, bool bOnlyPrimary, double fMinMapQ )
    {
        return bOnlyPrimary ? xGetPrimCoverage.scalar( from, to, iSeqId, fMinMapQ )
                            : xGetCoverage.scalar( from, to, iSeqId, fMinMapQ );
    }
}; // class

} // namespace libMSV
