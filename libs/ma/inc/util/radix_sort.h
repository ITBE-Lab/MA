/**
 * @file radix_sort.h
 * @brief Implements a radix sorting.
 * @author Arne Kutzner
 * @details
 * This code is adapted from KRADIX_SORT from ksw
 */
#ifndef RADIX_SORT_H
#define RADIX_SORT_H


/* Classic Insertion-Sort */
template <typename TP_ELE, typename FUN> void rs_insertsort( TP_ELE* beg, TP_ELE* end, FUN&& rskey )
{
    TP_ELE* i;
    for( i = beg + 1; i < end; ++i )
        if( rskey( *i ) < rskey( *( i - 1 ) ) )
        {
            TP_ELE *j, tmp = *i;
            for( j = i; j > beg && rskey( tmp ) < rskey( *( j - 1 ) ); --j )
                *j = *( j - 1 );
            *j = tmp;
        }
} // function

/* Radixsort
 * Algorithmic approach:
 * 1. Sort via bucketing according to the most significant column (column size in bits: BITS_PER_COL).
 * 2. Recursively sort all buckets using the next lower column or call insertion-sort in the case of few elements.
 * Algorithmically, this resembles a bucket driven Quicksort.
 * Due to the involved swapping (reordering of input vector), this sorting is not stable.
 * Further, this sorting can inly be used for unsigned integers
 */
#define BUCKETS_ON_STACK 1 // buckets on heap helps saving stackspace but is slightly slower
template <unsigned int BITS_PER_COL, // number of bits used per column in radix-sort (see Cormen etc.)
          unsigned int RADIX_SORT_MIN_SIZE, // minimum bucket size before using insertion-sort
          typename TP_ELE, // data-type of sorted elements
          typename FUN // function that delivers key-part of TP_ELE instances
          >
void rs_sort( TP_ELE* pVecBeg, TP_ELE* pVecEnd, unsigned int uiStartBit, FUN&& rskey )
{
    struct Bucket
    {
        TP_ELE *b, *e; // begin and pEnd of bucket
    }; // struct

    const unsigned int size = 1 << BITS_PER_COL; // bucket size
    const unsigned int m = size - 1; // bit mask

#if( BUCKETS_ON_STACK == 1 )
    static_assert( BITS_PER_COL <= 8 ); // if the array becomes too large, we risk a stack overflow
    Bucket pBuckVecBeg[ size ];
#else
    std::vector<Bucket> vBuckets( size );
    Bucket* const pBuckVecBeg = &vBuckets[ 0 ];
#endif
    const Bucket* const pBuckVecEnd = pBuckVecBeg + size;

    // iterate over all buckets for initializing
    for( auto k = pBuckVecBeg; k != pBuckVecEnd; ++k )
        k->b = k->e = pVecBeg;

    // iterate over the input vector and do the counting
    for( TP_ELE* i = pVecBeg; i != pVecEnd; ++i )
        ++pBuckVecBeg[ rskey( *i ) >> uiStartBit & m ].e;

    // pBuckVecBeg[n].e - pBuckVecBeg[n].pBuckVecBeg delivers the size of bucket n now
    // iterate over all buckets and sum up
    for( auto k = pBuckVecBeg + 1; k != pBuckVecEnd; ++k )
    {
        k->e += ( k - 1 )->e - pVecBeg;
        k->b = ( k - 1 )->e;
    } // for

    // reorder input-vector in-place according to the bucket sizes
    for( auto k = pBuckVecBeg; k != pBuckVecEnd; )
    {
        if( k->b != k->e )
        {
            Bucket* l;
            if( ( l = pBuckVecBeg + ( rskey( *k->b ) >> uiStartBit & m ) ) != k )
            {
                TP_ELE tmp = *k->b, swap;
                do
                {
                    swap = tmp;
                    tmp = *l->b;
                    *l->b++ = swap;
                    l = pBuckVecBeg + ( rskey( tmp ) >> uiStartBit & m );
                } while( l != k );
                *k->b++ = tmp;
            } // if
            else
                ++k->b;
        } // if
        else
            ++k;
    } // for
    {
        auto k = pBuckVecBeg + 1;
        for( pBuckVecBeg->b = pVecBeg; k != pBuckVecEnd; ++k )
            k->b = ( k - 1 )->e;
    } // var scope

    if( uiStartBit ) // still some bits left?
    {
        uiStartBit = uiStartBit > BITS_PER_COL ? uiStartBit - BITS_PER_COL : 0;
        // iterate over all buckets and sort them
        for( auto k = pBuckVecBeg; k != pBuckVecEnd; ++k )
            if( k->e - k->b > RADIX_SORT_MIN_SIZE )
                // the bucket comprises more than RS_MIN_SIZE many elements.
                // Sort them stable via a recursive call, with uiStartBit as position of the leading bit.
                rs_sort<BITS_PER_COL, RADIX_SORT_MIN_SIZE>( k->b, k->e, uiStartBit, rskey );

            else if( k->e - k->b > 1 )
                // the bucket comprises less than RS_MIN_SIZE many elements.
                // Sort them stable via insertion sort.
                rs_insertsort( k->b, k->e, rskey );
    } // if
} // function


// Example TP_ELE : uint64_t
//    fExtractKey : [](uint64_t v){return v;} // delivers the key for an instance
template <unsigned int BITS_PER_COL = 8, unsigned int RADIX_SORT_MIN_SIZE = 64, typename TP_ELE, typename FUNCTOR>
void radix_sort( TP_ELE* pBegin, TP_ELE* pEnd, FUNCTOR&& fExtractKey )
{
    if( pEnd - pBegin <= RADIX_SORT_MIN_SIZE )
        rs_insertsort( pBegin, pEnd, fExtractKey );
    else
    {
        typedef decltype( fExtractKey( *pBegin ) ) TypeOfKey;
        static_assert( !std::is_signed<TypeOfKey>::value ); // this radix sort works for unsigned keys only

        const unsigned int uiKeySizeInBits = sizeof( TypeOfKey ) * 8;
        static_assert( uiKeySizeInBits % BITS_PER_COL == 0 ); // size of key must pBuckVecEnd multiple of column size

        rs_sort<BITS_PER_COL, RADIX_SORT_MIN_SIZE>( pBegin, pEnd, uiKeySizeInBits - BITS_PER_COL, fExtractKey );
    } // else
} // function

#endif