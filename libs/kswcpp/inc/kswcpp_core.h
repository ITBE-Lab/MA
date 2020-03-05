/* Core of the ksw module.
 * The core does not include the SIMD classes. This has to be done by the core's user.
 */
#pragma once

#include <assert.h>
#include <limits>

/* TO DO: remove */
#define kmalloc( km, size ) malloc( ( size ) )
#define kcalloc( km, count, size ) calloc( ( count ), ( size ) )
#define krealloc( km, ptr, size ) realloc( ( ptr ), ( size ) )
#define kfree( km, ptr ) free( ( ptr ) )

static inline void ksw_reset_extz( kswcpp_extz_t* ez )
{
    ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
    ez->max = 0, ez->score = ez->mqe = ez->mte = std::numeric_limits<int>::min( ); // KSW_NEG_INF;
    ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}

static inline int ksw_apply_zdrop( kswcpp_extz_t* ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e )
{
    int r, t;
    if( is_rot )
        r = a, t = b;
    else
        r = a + b, t = a;
    if( H > (int32_t)ez->max )
    {
        ez->max = H, ez->max_t = t, ez->max_q = r - t;
    }
    else if( t >= ez->max_t && r - t >= ez->max_q )
    {
        int tl = t - ez->max_t, ql = ( r - t ) - ez->max_q, l;
        l = tl > ql ? tl - ql : ql - tl;
        if( zdrop >= 0 && ( int32_t )( ez->max - H ) > zdrop + l * e )
        {
            ez->zdropped = 1;
            return 1;
        }
    }
    return 0;
}

static inline uint32_t* ksw_push_cigar( void* km, int* n_cigar, int* m_cigar, uint32_t* cigar, uint32_t op, int len )
{
    //// printf("ksw_push_cigar in\n");
    if( *n_cigar == 0 || op != ( cigar[ ( *n_cigar ) - 1 ] & 0xf ) )
    {
        //// printf("ksw_push_cigar realloc start, %d, %d\n", *n_cigar, *m_cigar);
        if( *n_cigar == *m_cigar )
        {
            *m_cigar = *m_cigar ? ( *m_cigar ) << 1 : 4;
            cigar = (uint32_t*)krealloc( km, cigar, ( *m_cigar ) << 2 );
        }
        //// printf("ksw_push_cigar realloc start, %d\n", *n_cigar );
        cigar[ ( *n_cigar )++ ] = len << 4 | op;
        //// printf("ksw_push_cigar realloc end, %d\n", cigar == NULL);
    }
    else
        cigar[ ( *n_cigar ) - 1 ] += len << 4;
    //// printf("ksw_push_cigar out\n");
    return cigar;
}

// In the backtrack matrix, value p[] has the following structure:
//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
template <typename TP_DIFF_SCORE_UNSIGNED>
static inline void ksw_backtrack__( void* km, int is_rot, int is_rev, int min_intron_len,
                                    const TP_DIFF_SCORE_UNSIGNED* p, const int* off, const int* off_end, int n_col,
                                    int i0, int j0, int* m_cigar_, int* n_cigar_, uint32_t** cigar_ )
{ // p[] - lower 3 bits: which type gets the max; bit
    int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
    uint32_t *cigar = *cigar_, tmp;
    while( i >= 0 && j >= 0 )
    { // at the beginning of the loop, _state_ tells us which state to check
        int force_state = -1;
        if( is_rot )
        {
            r = i + j;
            if( i < off[ r ] )
                force_state = 2;
            if( off_end && i > off_end[ r ] )
                force_state = 1;
            tmp = force_state < 0 ? p[ r * n_col + i - off[ r ] ] : 0;
        }
        else
        {
            if( j < off[ i ] )
                force_state = 2;
            if( off_end && j > off_end[ i ] )
                force_state = 1;
            tmp = force_state < 0 ? p[ i * n_col + j - off[ i ] ] : 0;
        }

        if( state == 0 )
            state = tmp & 7; // if requesting the H state, find state one maximizes it.
        else if( !( tmp >> ( state + 2 ) & 1 ) )
            state =
                0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H

        if( state == 0 )
            state =
                tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
        if( force_state >= 0 )
            state = force_state;
        if( state == 0 )
            cigar = ksw_push_cigar( km, &n_cigar, &m_cigar, cigar, 0, 1 ), --i, --j; // match
        else if( state == 1 || ( state == 3 && min_intron_len <= 0 ) )
            cigar = ksw_push_cigar( km, &n_cigar, &m_cigar, cigar, 2, 1 ), --i; // deletion
        else if( state == 3 && min_intron_len > 0 )
            cigar = ksw_push_cigar( km, &n_cigar, &m_cigar, cigar, 3, 1 ), --i; // intron
        else
            cigar = ksw_push_cigar( km, &n_cigar, &m_cigar, cigar, 1, 1 ), --j; // insertion
    } // while
    if( i >= 0 )
        cigar = ksw_push_cigar( km, &n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len ? 3 : 2,
                                i + 1 ); // first deletion
    if( j >= 0 )
        cigar = ksw_push_cigar( km, &n_cigar, &m_cigar, cigar, 1, j + 1 ); // first insertion
    if( !is_rev )
        for( i = 0; i<n_cigar>> 1; ++i ) // reverse CIGAR
            tmp = cigar[ i ], cigar[ i ] = cigar[ n_cigar - 1 - i ], cigar[ n_cigar - 1 - i ] = tmp;
    *m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
} // function

#define USE_HMAX // compute the maximum on SIMD-vector level

/* Codepart responsible for processing the difference values in the context of scoring
 */
template <typename T_SIMD_VEC, typename T_ScoreSmall, typename T_SCORES>
static inline bool calcMaxScore( const int& approx_max, // boolean switch
                                 const int& r, // row
                                 const int& st0, // aligned start of SIMD-Vec
                                 const int& en0, // aligned end of SIMD-Vec
                                 const T_ScoreSmall* const& u8, // difference value
                                 const T_ScoreSmall* const& v8, // difference value
                                 const int& en, // reference end (index), >= 0
                                 const int& qe, // query end (index)
                                 const int& qlen, // query len
                                 const int& tlen, // (target) reference len
                                 const int& zdrop, // zdrop parameter
                                 const int8_t& e2, // gap-extend 2 parameter
                                 const int flag, // control flags ( only used in the case of approximated scores)
                                 T_SCORES& H0, // only used in the case of approximated scores
                                 T_SCORES& last_H0_t, // only used in the case of approximated scores
                                 T_SCORES* const& H, // H-vector with scores (32 bit etc., size equal to reference size)
                                 kswcpp_extz_t* const& ez )
{
    static_assert( ( sizeof( T_SIMD_VEC ) / sizeof( T_SCORES ) == T_SIMD_VEC::SIZE ) &&
                   ( sizeof( T_SIMD_VEC ) % sizeof( T_SCORES ) == 0 ) );

    if( !approx_max )
    { // find the exact max with a 32-bit score array
      // This is very expensive and consumes a lot of runtime
        T_SCORES max_H, max_t; // type is the signed scoring type (standard is 32 bit)
                               // compute H[], max_H and max_t
        if( r > 0 )
        {
#ifdef USE_HMAX
            T_SIMD_VEC v_HH, v_tt;
#else
            int32_t HH[ 4 ], tt[ 4 ];
#endif // USE_HMAX
            size_t en1 = st0 + ( ( en0 - st0 ) / T_SIMD_VEC::SIZE ) * T_SIMD_VEC::SIZE; // compute aligned end point

            T_SIMD_VEC v_max_H_, v_max_t_;

            max_H = H[ en0 ] =
                en0 > 0 ? H[ en0 - 1 ] + u8[ en0 ] : H[ en0 ] + v8[ en0 ]; // special casing the last element
            max_t = en0;

            v_max_H_.set1( max_H );
            v_max_t_.set1( max_t );

            T_SCORES t; // column index; used continuously after for loop; However, this index must fit into T_SCORES!
            for( t = static_cast<T_SCORES>( st0 ); t < static_cast<T_SCORES>( en1 ); t += T_SIMD_VEC::SIZE )
            { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;

                T_SIMD_VEC v_H1, v_tmp; //  , v_t_;

                v_H1.load( &H[ t ] );
                T_SIMD_VEC v_t_ =
                    T_SIMD_VEC::expand_int8_t( &v8[ t ] ); // v_t_.set4( v8[t], v8[t+1], v8[t+2], v8[t+3] );
                v_H1 = v_H1 + v_t_;
                v_H1.store( &H[ t ] );
                v_t_.set1( t );
                v_tmp = v_H1 > v_max_H_;
                v_max_H_ = v_max_H_.blendv( v_H1, v_tmp );
                v_max_t_ = v_max_t_.blendv( v_t_, v_tmp );

            } // for t

#ifdef USE_HMAX
            v_HH = v_max_H_;
            v_tt = v_max_t_;
#else
            v_max_H_.store( HH );
            v_max_t_.store( tt );
#endif // USE_HMAX

            // iterate over elements of SIMD vector (horizontal maximum)
#ifdef USE_HMAX
            max_H = v_HH.hmax( );
            max_t = v_tt.hmax( );
#else
            for( int i = 0; i < 4; ++i )
                if( max_H < HH[ i ] )
                    max_H = HH[ i ], max_t = tt[ i ] + i;
#endif // USE_HMAX

            // for the rest of values that haven't been computed with SSE
            for( ; t < static_cast<T_SCORES>( en0 ); ++t )
            {
                H[ t ] += (T_SCORES)v8[ t ];
                if( H[ t ] > max_H )
                    max_H = H[ t ], max_t = t;
            } // for
        } // if
        else
        { // special case r==0 (first row); update ez
            H[ 0 ] = v8[ 0 ] - qe;
            max_H = H[ 0 ];
            max_t = 0;
        } // else

        // adapt mte (max score when reaching the end of target == reference) and
        // adapt mqe (max score when reaching the end of query)
        if( en0 == tlen - 1 && H[ en0 ] > ez->mte )
            ez->mte = H[ en0 ], ez->mte_q = r - en;
        if( r - st0 == qlen - 1 && H[ st0 ] > ez->mqe )
            ez->mqe = H[ st0 ], ez->mqe_t = st0;

        // check for break due to zdrop
        if( ksw_apply_zdrop( ez, 1, max_H, r, max_t, zdrop, e2 ) )
            return true;
        if( r == qlen + tlen - 2 && en0 == tlen - 1 )
            ez->score = H[ tlen - 1 ];
    } // if (approximate max)
    else
    { /* find approximate max; Z-drop might be inaccurate, too. */
        if( r > 0 )
        {
            if( last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0 )
            {
                T_SCORES d0 = v8[ last_H0_t ]; // orig. int32_t
                T_SCORES d1 = u8[ last_H0_t + 1 ]; // orig. int32_t
                if( d0 > d1 )
                    H0 += d0;
                else
                    H0 += d1, ++last_H0_t;
                // DEBUG std::cout << r << ' ' << H0 << ' ' << d0 << ' ' << d1 << std::endl;
            }
            else if( last_H0_t >= st0 && last_H0_t <= en0 )
            {
                H0 += v8[ last_H0_t ];
            }
            else
            {
                ++last_H0_t, H0 += u8[ last_H0_t ];
            }
        } // if
        else
            H0 = v8[ 0 ] - qe, last_H0_t = 0;

        /* check for break due to zdrop */
        if( ( flag & KSW_EZ_APPROX_DROP ) && ksw_apply_zdrop( ez, 1, H0, r, last_H0_t, zdrop, e2 ) )
            return true;
        if( r == qlen + tlen - 2 && en0 == tlen - 1 )
            ez->score = H0;
    } // else

    return false; // do not break
} // function

template <typename TP_SCORE, // type for scoring
          typename TP_SCORE_VEC, // SIMD vector-type for scoring
          typename TP_DIFF_SCORE, // signed type for differences
          typename TP_DIFF_SCORE_UNSIGNED, // unsigned type for differences
          typename TP_DIFF_VEC, // SIMD-vector type for differences
          bool LEFT_GAP_ALIGN = true,
          bool WITH_CIGAR = true>
static inline void kswcpp_inner_core( const int qlen, // query length
                                      const uint8_t* query, // query
                                      const int tlen, // reference length
                                      const uint8_t* target, // reference
#ifdef USE_CPP_PARAM
                                      const KswCppParam<5>& xParam,
#else
                                      int8_t m, const int8_t* mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
#endif
                                      int w, // bandwidth (if -1, we do an auto setting of the bandwidth)
                                      const int zdrop, //
                                      const int flag,
                                      kswcpp_extz_t* ez,
                                      AlignedMemoryManager& rxMemManager )
{
#ifdef USE_CPP_PARAM
    const int8_t& m = xParam.m;
    const int8_t* mat = &xParam.mat[ 0 ];
    int8_t q = xParam.q;
    int8_t e = xParam.e;
    int8_t q2 = xParam.q2;
    int8_t e2 = xParam.e2;
    // int w = xParam.w; // not part of params anymore
    // const int zdrop = xParam.zdrop;
    // const int flag = xParam.flag;
#endif

    /* Code start. */
    constexpr const TP_SCORE iKSW_NEG_INF = std::numeric_limits<TP_SCORE>::min( );

    int t, qe = q + e, n_col_, *off = 0, *off_end = 0, qlen_, /* last_st, */ last_en, wl, wr, max_sc, min_sc,
           long_thres, long_diff;

    // int with_cigar = !(flag & KSW_EZ_SCORE_ONLY);
    int approx_max = !!( flag & KSW_EZ_APPROX_MAX );

    TP_SCORE H0 = 0, last_H0_t = 0; // used in the context of approximative score computation (orig. int32_t)

    TP_DIFF_SCORE_UNSIGNED* qr; // reverse of the query
    TP_DIFF_SCORE_UNSIGNED* sf; // copy of reference
#ifdef OLD_MEMORY_MANAGMENT
    T_VecEleUnsigned* mem;
#endif
    TP_DIFF_SCORE_UNSIGNED* mem2 = 0;

    ksw_reset_extz( ez );
    if( m <= 1 || qlen <= 0 || tlen <= 0 )
        return;

    /* make sure q + e no larger than q2 + e2 */
    if( q2 + e2 < q + e )
    { /* swap q and q2 as well as q2 and e2 */
        t = q;
        q = q2;
        q2 = t;
        t = e;
        e = e2;
        e2 = t;
    } // if

    TP_DIFF_VEC q_struct;
    TP_DIFF_VEC q2_struct;
    TP_DIFF_VEC qe_struct;
    TP_DIFF_VEC qe2_struct;
    TP_DIFF_VEC zero_struct;
    TP_DIFF_VEC sc_mch_struct;
    TP_DIFF_VEC sc_mis_struct;

    zero_struct.setzero( );
    q_struct.set1( q ); // FIX ME: could be part of the parameter
    q2_struct.set1( q2 );
    qe_struct.set1( q + e );
    qe2_struct.set1( q2 + e2 );
    sc_mch_struct.set1( mat[ 0 ] ); // vector keeping match score
    sc_mis_struct.set1( mat[ 1 ] ); // vector keeping mismatch score

    // Auto-bandwidth setting
    if( w < 0 )
        w = tlen > qlen ? tlen : qlen;

    wl = wr = w;

    /* tlen_ = size of reference expressed in SIMD-arrays */
    const int tlen_ = ( tlen + ( TP_DIFF_VEC::SIZE - 1 ) ) / TP_DIFF_VEC::SIZE;
    n_col_ = qlen < tlen ? qlen : tlen;
    n_col_ = ( ( n_col_ < w + 1 ? n_col_ : w + 1 ) + ( TP_DIFF_VEC::SIZE - 1 ) ) / TP_DIFF_VEC::SIZE + 1;
    qlen_ = ( qlen + ( TP_DIFF_VEC::SIZE - 1 ) ) / TP_DIFF_VEC::SIZE;
    // DEBUG std::cout << "n_col_:" << n_col_ << std::endl;

    for( t = 1, max_sc = mat[ 0 ], min_sc = mat[ 1 ]; t < m * m; ++t )
    {
        max_sc = max_sc > mat[ t ] ? max_sc : mat[ t ];
        min_sc = min_sc < mat[ t ] ? min_sc : mat[ t ];
    }
    if( -min_sc > 2 * ( q + e ) )
        return; // otherwise, we won't see any mismatches

    long_thres = e != e2 ? ( q2 - q ) / ( e - e2 ) - 1 : 0;
    if( q2 + e2 + long_thres * e2 > q + e + long_thres * e )
        ++long_thres;
    long_diff = long_thres * ( e - e2 ) - ( q2 - q ) - e2;

#ifdef OLD_MEMORY_MANAGMENT
    mem = (uint8_t*)calloc( tlen_ * 8 + qlen_ + 1, 16 );
    u_pstruct = (struct__mXXXi*)( ( (size_t)mem + 15 ) >> 4 << 4 ); // 16-byte aligned
#else
    TP_DIFF_VEC* u_pstruct;
    TP_DIFF_VEC* v_pstruct;
    TP_DIFF_VEC* x_pstruct;
    TP_DIFF_VEC* y_pstruct;
    TP_DIFF_VEC* x2_pstruct;
    TP_DIFF_VEC* y2_pstruct;
    TP_DIFF_VEC* s_pstruct;
    TP_DIFF_VEC* p_pstruct = NULL;

    u_pstruct = rxMemManager.reserveMemMatrix<TP_DIFF_VEC>( ( tlen_ * 8 + qlen_ ) );
    v_pstruct = u_pstruct + tlen_; // aligned to SIMD-type

#endif // OLD_MEMORY_MANAGMENT
    x_pstruct = v_pstruct + tlen_; // aligned to SIMD-type
    y_pstruct = x_pstruct + tlen_; // aligned to SIMD-type
    x2_pstruct = y_pstruct + tlen_; // aligned to SIMD-type
    y2_pstruct = x2_pstruct + tlen_; // aligned to SIMD-type
    s_pstruct = y2_pstruct + tlen_; // aligned to SIMD-type
    sf = (TP_DIFF_SCORE_UNSIGNED*)( s_pstruct + tlen_ ); // sf is a reference copy in T_VecEleUnsigned
    qr = sf + tlen_ * TP_DIFF_VEC::SIZE; // qr is the query reverse in T_VecEleUnsigned

#if 0
										 /* Initialize gap lengths */
	struct__mXXXi vGap; vGap.set1( -q - e );
	struct__mXXXi vGap2; vGap2.set1( -q2 - e2 );

	for(size_t uiItr = 0; uiItr < (size_t)tlen_; uiItr++)
	{
		u_pstruct[uiItr] = vGap;
		v_pstruct[uiItr] = vGap;
		x_pstruct[uiItr] = vGap;
		y_pstruct[uiItr] = vGap;
		x2_pstruct[uiItr] = vGap2;
		y2_pstruct[uiItr] = vGap2;
	} // for
#else
    /* In order to use memset, difference type must be 8-bit.
     * memset is a bit faster than above for loop.
     */
    static_assert( sizeof( TP_DIFF_VEC ) == TP_DIFF_VEC::SIZE );

    memset( &u_pstruct[ 0 ], -q - e, tlen_ * TP_DIFF_VEC::SIZE );
    memset( &v_pstruct[ 0 ], -q - e, tlen_ * TP_DIFF_VEC::SIZE );
    memset( &x_pstruct[ 0 ], -q - e, tlen_ * TP_DIFF_VEC::SIZE );
    memset( &y_pstruct[ 0 ], -q - e, tlen_ * TP_DIFF_VEC::SIZE );
    memset( &x2_pstruct[ 0 ], -q2 - e2, tlen_ * TP_DIFF_VEC::SIZE );
    memset( &y2_pstruct[ 0 ], -q2 - e2, tlen_ * TP_DIFF_VEC::SIZE );
#endif

    /* Initialize H-Vector, if not approximate max.
     */
    TP_SCORE* H = 0;
    if( !approx_max )
    {
#ifdef OLD_MEMORY_MANAGMENT
        H = (T_ScrLrg*)malloc( tlen_ * struct__mXXXi::SIZE * sizeof( T_ScrLrg ) ); // 4 = sizeof(int32_t)
#else
        H = rxMemManager.reserveMemH_Vector<TP_SCORE>( tlen_ * TP_DIFF_VEC::SIZE );
#endif
        /* Could be done via SIMD as well... */
        for( t = 0; t < tlen_ * (int)TP_DIFF_VEC::SIZE; ++t )
            H[ t ] = iKSW_NEG_INF; // Set all H values to min()
    } // if

    if( WITH_CIGAR )
    { /* Allocate aligned memory for p */
        mem2 = (TP_DIFF_SCORE_UNSIGNED*)malloc( ( ( qlen + tlen - 1 ) * n_col_ + 1 ) * sizeof( TP_DIFF_VEC ) );
        p_pstruct = (TP_DIFF_VEC*)( ( ( (size_t)mem2 + ( sizeof( TP_DIFF_VEC ) - 1 ) ) / sizeof( TP_DIFF_VEC ) ) *
                                    sizeof( TP_DIFF_VEC ) );

        /* Allocate memory for the cigar itself */
        off = (int*)malloc( ( qlen + tlen - 1 ) * sizeof( int ) * 2 );
        off_end = off + qlen + tlen - 1;
    } // if

    /* Create the reverse of the query.
     * Performs appropriate type conversion to T_VecEleUnsigned.
     */
    for( t = 0; t < qlen; ++t )
        qr[ t ] = query[ qlen - 1 - t ];

    /* Copy the reference to sf.
     * WARNING: If the difference type is not 8 bit we have to do a type conversion here.
     */
    memcpy( sf, target, tlen ); // copy the reference to target

    /* Central iteration */
    for( int r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r )
    {
        int st = 0; // start reference
        int en = tlen - 1; // end
        int st0; // SIMD-Vec aligned start
        int en0; // SIMD-Vec aligned end
        int st_, en_;

        TP_DIFF_SCORE x1; // originally int8_t
        TP_DIFF_SCORE x21; // originally int8_t
        TP_DIFF_SCORE v1; // originally int8_t

        TP_DIFF_SCORE_UNSIGNED* qrr = qr + ( qlen - 1 - r ); // originally uint8_t

        TP_DIFF_SCORE* u8 = (TP_DIFF_SCORE*)u_pstruct; // originally int8_t
        TP_DIFF_SCORE* v8 = (TP_DIFF_SCORE*)v_pstruct; // originally int8_t
        TP_DIFF_SCORE* x8 = (TP_DIFF_SCORE*)x_pstruct; // originally int8_t
        TP_DIFF_SCORE* x28 = (TP_DIFF_SCORE*)x2_pstruct; // originally int8_t

        TP_DIFF_VEC x1_struct;
        TP_DIFF_VEC x21_struct;
        TP_DIFF_VEC v1_struct;

        /* find the boundaries */
        if( st < r - qlen + 1 )
            st = r - qlen + 1;
        if( en > r )
            en = r;
        if( st<( r - wr + 1 )>> 1 )
            st = ( r - wr + 1 ) >> 1; // take the ceil
        if( en > ( r + wl ) >> 1 )
            en = ( r + wl ) >> 1; // take the floor
        if( st > en )
        { /* Out of band drop */
            ez->zdropped = 1;
            break;
        } // if

        st0 = st; // st is start reference
        en0 = en;

        st = ( st / TP_DIFF_VEC::SIZE ) * TP_DIFF_VEC::SIZE; // align st with respect to the SIMD-vector
        en = ( en + TP_DIFF_VEC::SIZE ) / TP_DIFF_VEC::SIZE * TP_DIFF_VEC::SIZE - 1;

        /* set boundary conditions */
        if( st > 0 )
        {
            if( st - 1 >= last_st && st - 1 <= last_en )
            {
                x1 = x8[ st - 1 ], x21 = x28[ st - 1 ], v1 = v8[ st - 1 ]; // (r-1,s-1) calculated in the last round
            } // if
            else
            {
                x1 = -q - e, x21 = -q2 - e2;
                v1 = -q - e;
            } // else
        } // if
        else
        {
            x1 = -q - e;
            x21 = -q2 - e2;
            v1 = r == 0 ? -q - e : r < long_thres ? -e : r == long_thres ? long_diff : -e2; // otherwise
        } // else
        if( en >= r )
        {
            ( (TP_DIFF_SCORE*)y_pstruct )[ r ] = -q - e;
            ( (TP_DIFF_SCORE*)y2_pstruct )[ r ] = -q2 - e2; // DEBUG: y2 is the troublemaker
            u8[ r ] = r == 0 ? -q - e : r < long_thres ? -e : r == long_thres ? long_diff : -e2;
        } // if

        /* loop fission: set scores first
         * WARNING: Here is a hidden bug in the case of no generic scoring (if branch)
         *          The final execution of the for-loop can write beyond the end s_pstruct
         */
        if( !( flag & KSW_EZ_GENERIC_SC ) )
        { /* Scoring only by match and mismatch */
            TP_DIFF_VEC m1_struct;
            TP_DIFF_VEC sc_N_struct;
            m1_struct.set1( m - 1 ); // m - 1 is the code of the N wildcard
            sc_N_struct.set1( -e2 ); // score for N wildcard, gap extension

            for( t = st0; t <= en0; t += TP_DIFF_VEC::SIZE )
            {
                TP_DIFF_VEC sq_struct;
                TP_DIFF_VEC st_struct;
                TP_DIFF_VEC tmp_struct;
                TP_DIFF_VEC mask_struct;

                sq_struct.load( &sf[ t ] ); // unaligned read (load reference at position t)
                st_struct.load( &qrr[ t ] ); // unaligned read
                mask_struct = ( sq_struct == m1_struct ) | ( st_struct == m1_struct );
                tmp_struct = ( sq_struct == st_struct );

                /* Blend in match and mismatch score */
                tmp_struct = sc_mis_struct.blendv( sc_mch_struct, tmp_struct );
                tmp_struct = tmp_struct.blendv( sc_N_struct, mask_struct );

                tmp_struct.store( (TP_DIFF_SCORE*)s_pstruct + t );
            } // for
        } // if
        else
        { /* Use the matrix values for the scoring profile */
            for( t = st0; t <= en0; ++t )
                ( (TP_DIFF_SCORE_UNSIGNED*)s_pstruct )[ t ] =
                    mat[ sf[ t ] * m + qrr[ t ] ]; // strange; here is cast to unsigned element type (uint8_t)
        } // else

        /* loop core */
        x1_struct.moveAndZero32( (uint8_t)x1 ); // uint8_t is needed for plain bit copy (no evaluation of sign)
        x21_struct.moveAndZero32( (uint8_t)x21 );
        v1_struct.moveAndZero32( (uint8_t)v1 );

        st_ = st / TP_DIFF_VEC::SIZE;
        en_ = en / TP_DIFF_VEC::SIZE;

        assert( en_ - st_ + 1 <= n_col_ );

        TP_DIFF_VEC* pr_pstruct = nullptr; // used in the context of cigar construction merely
        if( WITH_CIGAR )
        { /* Only requied for cigar */
            pr_pstruct = (TP_DIFF_VEC*)( p_pstruct + ( ( r * n_col_ ) - st_ ) );
            off[ r ] = st;
            off_end[ r ] = en;
        } // if

        TP_DIFF_VEC d_struct; // used in the context of cigar construction merely
        for( t = st_; t <= en_; ++t )
        { /* DP Block 1 */
            TP_DIFF_VEC z_struct = s_pstruct[ t ];
            TP_DIFF_VEC xt1_struct = x_pstruct[ t ].template shiftl<1>( ) | x1_struct;
            x1_struct = x_pstruct[ t ].template shiftr<TP_DIFF_VEC::SIZE - 1>( );
            TP_DIFF_VEC vt1_struct = v_pstruct[ t ].template shiftl<1>( ) | v1_struct;
            v1_struct = v_pstruct[ t ].template shiftr<TP_DIFF_VEC::SIZE - 1>( );
            TP_DIFF_VEC a_struct = xt1_struct + vt1_struct;
            TP_DIFF_VEC ut_struct = u_pstruct[ t ];
            TP_DIFF_VEC b_struct = y_pstruct[ t ] + u_pstruct[ t ];
            TP_DIFF_VEC x2t1_struct = x21_struct | x2_pstruct[ t ].template shiftl<1>( );
            x21_struct = x2_pstruct[ t ].template shiftr<TP_DIFF_VEC::SIZE - 1>( );
            TP_DIFF_VEC a2_struct = x2t1_struct + vt1_struct;
            TP_DIFF_VEC b2_struct = y2_pstruct[ t ] + u_pstruct[ t ];

            if( WITH_CIGAR )
            {
                if( LEFT_GAP_ALIGN )
                    d_struct = TP_DIFF_VEC( 1 ) & ( a_struct > z_struct );
                else
                    d_struct = ( z_struct > a_struct ).andnot( TP_DIFF_VEC( 1 ) );
            } // if

            z_struct = z_struct.max( a_struct );
            if( WITH_CIGAR )
            {
                if( LEFT_GAP_ALIGN )
                    d_struct = d_struct.blendv( TP_DIFF_VEC( 2 ), b_struct > z_struct );
                else
                    d_struct = TP_DIFF_VEC( 2 ).blendv( d_struct, z_struct > b_struct );
            } // if
            z_struct = z_struct.max( b_struct );
            if( WITH_CIGAR )
            {
                if( LEFT_GAP_ALIGN )
                    d_struct = d_struct.blendv( TP_DIFF_VEC( 3 ), a2_struct > z_struct );
                else
                    d_struct = TP_DIFF_VEC( 3 ).blendv( d_struct, z_struct > a2_struct );
            } // if
            z_struct = z_struct.max( a2_struct );
            if( WITH_CIGAR )
            {
                if( LEFT_GAP_ALIGN )
                    d_struct = d_struct.blendv( TP_DIFF_VEC( 4 ), b2_struct > z_struct );
            } // if
            else
                d_struct = TP_DIFF_VEC( 4 ).blendv( d_struct, z_struct > b2_struct );
            z_struct = z_struct.max( b2_struct );

            z_struct = z_struct.min( sc_mch_struct ); //  match score

            /* DP Block 2*/
            u_pstruct[ t ] = z_struct - vt1_struct;
            v_pstruct[ t ] = z_struct - ut_struct;
            TP_DIFF_VEC tmp_struct = z_struct - q_struct; // Debug: q OK
            a_struct = a_struct - tmp_struct;

            b_struct = b_struct - tmp_struct;
            tmp_struct = z_struct - q2_struct;
            a2_struct = a2_struct - tmp_struct;
            b2_struct = b2_struct - tmp_struct;

            if( WITH_CIGAR )
            { /* Prepare cigar */
                if( LEFT_GAP_ALIGN )
                { /* Left aligned gaps */
                    tmp_struct = a_struct > zero_struct;
                    x_pstruct[ t ] = ( tmp_struct & a_struct ) - qe_struct;

                    d_struct = d_struct | ( tmp_struct & TP_DIFF_VEC( 0x08 ) );
                    tmp_struct = LEFT_GAP_ALIGN ? b_struct > zero_struct : zero_struct > b_struct;
                    y_pstruct[ t ] = ( tmp_struct & b_struct ) - qe_struct;
                    d_struct = d_struct | ( tmp_struct & TP_DIFF_VEC( 0x10 ) );
                    tmp_struct = LEFT_GAP_ALIGN ? a2_struct > zero_struct : zero_struct > a2_struct;
                    x2_pstruct[ t ] = ( tmp_struct & a2_struct ) - qe2_struct;
                    d_struct = d_struct | ( tmp_struct & TP_DIFF_VEC( 0x20 ) );
                    tmp_struct = LEFT_GAP_ALIGN ? b2_struct > zero_struct : zero_struct > b2_struct;
                    y2_pstruct[ t ] = ( tmp_struct & b2_struct ) - qe2_struct;
                    d_struct = d_struct | ( tmp_struct & TP_DIFF_VEC( 0x40 ) );
                    pr_pstruct[ t ] = d_struct;
                } // if
                else
                { /* Right aligned gaps */
                    tmp_struct = zero_struct > a_struct;
                    x_pstruct[ t ] = ( tmp_struct.andnot( a_struct ) ) - qe_struct;
                    d_struct = d_struct | ( tmp_struct.andnot( TP_DIFF_VEC( 0x08 ) ) );
                    tmp_struct = zero_struct > b_struct;
                    y_pstruct[ t ] = ( tmp_struct.andnot( b_struct ) ) - qe_struct;
                    d_struct = d_struct | ( tmp_struct.andnot( TP_DIFF_VEC( 0x10 ) ) );
                    tmp_struct = zero_struct > a2_struct;
                    x2_pstruct[ t ] = ( tmp_struct.andnot( a2_struct ) ) - qe2_struct;
                    d_struct = d_struct | ( tmp_struct.andnot( TP_DIFF_VEC( 0x20 ) ) );
                    tmp_struct = zero_struct > b2_struct;
                    y2_pstruct[ t ] = ( tmp_struct.andnot( b2_struct ) ) - qe2_struct;
                    d_struct = d_struct | ( tmp_struct.andnot( TP_DIFF_VEC( 0x40 ) ) );
                    pr_pstruct[ t ] = d_struct;
                } // else
            } // if
            else
            { /* Compute the scores merely */
                x_pstruct[ t ] = a_struct.max( zero_struct ) - qe_struct;
                y_pstruct[ t ] = b_struct.max( zero_struct ) - qe_struct;
                x2_pstruct[ t ] = a2_struct.max( zero_struct ) - qe2_struct;
                y2_pstruct[ t ] = b2_struct.max( zero_struct ) - qe2_struct;
            } // else
        } // for

        bool bStop = calcMaxScore<TP_SCORE_VEC>( approx_max, // approx_max, // boolean switch
                                                 r, // row
                                                 st0, // SIMD-Vec aligned start
                                                 en0, // SIMD-Vec aligned end
                                                 u8, // difference values
                                                 v8, // difference values
                                                 en, // reference end (index)
                                                 qe, // q + e; gap_open + gep_extend
                                                 qlen, // query len
                                                 tlen, // reference len
                                                 zdrop, // zdrop parameter
                                                 e2, // gap-extend 2 param
                                                 flag, // control flags
                                                 H0, // only used in the case of approximated scores
                                                 last_H0_t, // only used in the case of approximated scores
                                                 H, // H-vector with 32 bit scores (modifiable)
                                                 ez ); // receives alignment outcome

        if( bStop )
            break;

        last_st = st, last_en = en;
    } // for r
#ifdef OLD_MEMORY_MANAGMENT
    kfree( km, mem );
    if( !approx_max )
        kfree( km, H );
#endif
    if( WITH_CIGAR )
    { /* Cigar construction by backtracking */
        int rev_cigar = !!( flag & KSW_EZ_REV_CIGAR );
        if( !ez->zdropped && !( flag & KSW_EZ_EXTZ_ONLY ) )
        { /* Backtracking up to tlen - 1, qlen - 1 */
            ksw_backtrack__( nullptr, 1, rev_cigar, 0, (TP_DIFF_SCORE_UNSIGNED*)p_pstruct, off, off_end,
                             n_col_ * TP_DIFF_VEC::SIZE, tlen - 1, qlen - 1, &ez->m_cigar, &ez->n_cigar, &ez->cigar );
        } // if
        else if( !ez->zdropped && ( flag & KSW_EZ_EXTZ_ONLY ) && ez->mqe > (int)( ez->max ) ) // AK: int cast added
        { /* Backtracking up to mqe_t, qlen - 1 */
            ez->reach_end = 1;
            ksw_backtrack__( nullptr, 1, rev_cigar, 0, (TP_DIFF_SCORE_UNSIGNED*)p_pstruct, off, off_end,
                             n_col_ * TP_DIFF_VEC::SIZE, ez->mqe_t, qlen - 1, &ez->m_cigar, &ez->n_cigar, &ez->cigar );
        } // else if
        else if( ez->max_t >= 0 && ez->max_q >= 0 )
        { /* Backtracking up to max_t, max_q */
            ksw_backtrack__( nullptr, 1, rev_cigar, 0, (TP_DIFF_SCORE_UNSIGNED*)p_pstruct, off, off_end,
                             n_col_ * TP_DIFF_VEC::SIZE, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar );
        } // else if
        kfree( nullptr, mem2 );
        kfree( nullptr, off );
    } // if
} // function


template <typename TP_SCORE, // type for scoring
          typename TP_SCORE_VEC, // SIMD vector-type for scoring
          typename TP_DIFF_SCORE, // signed type for differences
          typename TP_DIFF_SCORE_UNSIGNED, // unsigned type for differences
          typename TP_DIFF_VEC // SIMD-vector type for differences
          >
void kswcpp_core( const int qlen, // query length
                  const uint8_t* query, // query
                  const int tlen, // reference length
                  const uint8_t* target, // reference
#ifdef USE_CPP_PARAM
                  const KswCppParam<5>& xParam,
#else
                  int8_t m, const int8_t* mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
#endif
                  int w, // bandwidth (if -1, we do an auto setting of the bandwidth)
                  const int zdrop, //
                  const int flag,
                  kswcpp_extz_t* ez,
                  AlignedMemoryManager& rxMemManager )
{
    bool bWithCigar = !( flag & KSW_EZ_SCORE_ONLY );
    bool bLeftGapAlign = !( flag & KSW_EZ_RIGHT );

    if( bWithCigar )
        if( bLeftGapAlign )
            kswcpp_inner_core<TP_SCORE, TP_SCORE_VEC, TP_DIFF_SCORE, TP_DIFF_SCORE_UNSIGNED, TP_DIFF_VEC, true, true>(
                qlen, query, tlen, target, xParam, w, zdrop, flag, ez, rxMemManager );
        else
            kswcpp_inner_core<TP_SCORE, TP_SCORE_VEC, TP_DIFF_SCORE, TP_DIFF_SCORE_UNSIGNED, TP_DIFF_VEC, false, true>(
                qlen, query, tlen, target, xParam, w, zdrop, flag, ez, rxMemManager );
    else
        /* Without cigar implies that the gap alignments is not significant */
        kswcpp_inner_core<TP_SCORE, TP_SCORE_VEC, TP_DIFF_SCORE, TP_DIFF_SCORE_UNSIGNED, TP_DIFF_VEC, true, false>(
            qlen, query, tlen, target, xParam, w, zdrop, flag, ez, rxMemManager );
} // function
