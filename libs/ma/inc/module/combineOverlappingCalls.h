#pragma once

#include "container/svDb.h"
namespace libMA
{
void combineOverlappingCalls( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId )
{
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, double> xQuery(
        *pDb->pDatabase,
        "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, score "
        "FROM sv_call_table "
        "WHERE sv_caller_run_id == ? " );


    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, double> xQuery2(
        *pDb->pDatabase,
        "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, score "
        "FROM sv_call_table "
        "WHERE sv_caller_run_id == ? "
        "AND from_pos + from_size >= ? "
        "AND to_pos + to_size >= ? "
        "AND from_pos <= ? "
        "AND to_pos <= ? "
        "AND id != ? "
        "AND switch_strand != ? " );


    CppSQLiteExtQueryStatement<int64_t> xDelete( *pDb->pDatabase,
                                                 "DELETE "
                                                 "FROM sv_call_table "
                                                 "WHERE id == ? " );

    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xQuerySupport(
        *pDb->pDatabase,
        "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
        "sv_jump_table.id "
        "FROM sv_call_support_table "
        "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
        "WHERE sv_call_support_table.call_id == ? " );

    auto xIt = xQuery.vExecuteAndReturnIterator( iSvCallerId );
    while( !xIt.eof( ) )
    {
        auto xTup = xIt.get( );

        auto xIt2 =
            xQuery2.vExecuteAndReturnIterator( iSvCallerId,
                                               std::get<1>( xTup ), // from_pos
                                               std::get<2>( xTup ), // to_pos
                                               std::get<1>( xTup ) + std::get<3>( xTup ), // from_pos + from_size
                                               std::get<2>( xTup ) + std::get<4>( xTup ), // to_pos + to_size
                                               std::get<0>( xTup ), // id
                                               std::get<5>( xTup ) // bSwitchStrand
            );

        if( !xIt2.eof( ) )
        {

            SvCall xPrim( std::get<1>( xTup ), // uiFromStart
                          std::get<2>( xTup ), // uiToStart
                          std::get<3>( xTup ), // uiFromSize
                          std::get<4>( xTup ), // uiToSize
                          std::get<5>( xTup ), // bSwitchStrand
                          std::get<7>( xTup ) // dScore
            );
            xPrim.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
            xPrim.iId = std::get<0>( xTup );
            auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup ) ) );
            while( !xSupportIterator.eof( ) )
            {
                auto xTup = xSupportIterator.get( );
                xPrim.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
                xPrim.vSupportingJumps.emplace_back( rParameters.getSelected( ), std::get<0>( xTup ),
                                                     std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                                                     std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ),
                                                     std::get<7>( xTup ) );
                xSupportIterator.next( );
            } // while
            while( !xIt2.eof( ) )
            {
                auto xTup2 = xIt2.get( );
                SvCall xSec( std::get<1>( xTup2 ), // uiFromStart
                             std::get<2>( xTup2 ), // uiToStart
                             std::get<3>( xTup2 ), // uiFromSize
                             std::get<4>( xTup2 ), // uiToSize
                             std::get<5>( xTup2 ), // bSwitchStrand
                             std::get<7>( xTup2 ) // dScore
                );
                xSec.pInsertedSequence = std::get<6>( xTup2 ).pNucSeq;
                xSec.iId = std::get<0>( xTup2 );
                auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup2 ) ) );
                while( !xSupportIterator.eof( ) )
                {
                    auto xTup = xSupportIterator.get( );
                    xSec.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
                    xSec.vSupportingJumps.emplace_back( rParameters.getSelected( ), std::get<0>( xTup ),
                                                        std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                                                        std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ),
                                                        std::get<7>( xTup ) );
                    xSupportIterator.next( );
                } // while

                xPrim.join( xSec );


                xDelete.bindAndExecQuery<>( xSec.iId );

                xIt2.next( );
            } // while
            xDelete.bindAndExecQuery<>( xPrim.iId );
            // @todo recompute smaller bounds
            pDb->pSvCallTable->insertCall( iSvCallerId, xPrim );
        } // if
        xIt.next( );

    } // while

} // function
}; // namespace libMA