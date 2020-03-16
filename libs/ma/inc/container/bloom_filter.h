#if 0
/**
 * @file bloom_filter.h
 * @brief implements a k-mer bloom filter as is required for Yuansheng Liu et al. 2019
 * @details
 * Paper title:
 * Fast detection of maximal exact matches via fixed sampling of query K-mers and Bloom filtering of index K-mers
 */
#pragma once

#include "container/container.h"
#include <bitset>

namespace libMA
{

/**
 * @brief a bloom filter for k-mers
 * @details
 * @param uiH number of hash functions
 * @param uiS size of the bit array
 * @param uiK k-mer size
 */
template <size_t uiH, size_t uiS, size_t uiK> class BloomFilter : public Container
{
  private:
    std::bitset<uiS> vContent;
    std::array<int64_t, 5> xSeedTable;

  public:
    BloomFilter( )
    {} // default constructor


    int64_t rollingHash( int64_t iLastHash, uint8_t uiNew, uint8_t uiOld )
    {
        return ( iLastHash << 1 ) ^ ( xSeedTable[ uiOld ] << uiK ) ^ xSeedTable[ uiNew ];
    } // method

    int64_t hash( uint8_t* puiNucs )
    {
        int64_t iRet = 0;
        for( size_t uiI = 0; uiI < uiK; uiI++ )
            iRet ^= xSeedTable[ puiNucs[ uiI ] ] << ( uiK - 1 - uiI );
        return iRet;
    }



}; // class

} // namespace libMA
#endif