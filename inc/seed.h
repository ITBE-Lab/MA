/** 
 * @file seed.h
 * @brief Implements Seeds.
 * @author Markus Schmidt
 */
#ifndef SEED_H
#define SEED_H

#include "container.h"
#include "interval.h"

///@brief any index on the query or reference nucleotide sequence is given in this datatype
typedef uint64_t nucSeqIndex;

/**
 * @brief A seed.
 * @details
 * A extracted seed, that comprises two intervals, one on the query one on the reference.
 * Both intervals are equal in size.
 * @note the overloaded functions of Interval refer to the Interval on the query.
 */
class Seed: public Container, public Interval<nucSeqIndex>
{
private:
    ///@brief the beginning of the match on the reference
    nucSeqIndex uiPosOnReference;

public:

    /**
     * @brief Creates a new Seed.
     */
    Seed(
            const nucSeqIndex uiPosOnQuery, 
            const nucSeqIndex uiLenght, 
            const nucSeqIndex uiPosOnReference
        )
            :
        Interval(uiPosOnQuery, uiLenght),
        uiPosOnReference(uiPosOnReference)
    {}//constructor
    
    /**
     * @brief Returns the beginning of the seed on the reference.
     */
    nucSeqIndex start_ref() const
    {
        return uiPosOnReference;
    }//function
    
    /**
     * @brief Returns the end of the seed on the reference.
     */
    nucSeqIndex end_ref() const
    {
        return uiPosOnReference + size();
    }//function
    
    /**
     * @brief Returns the value of the seed.
     * @details
     * A seeds value corresponds to its size.
     */
    nucSeqIndex getValue() const
    {
        return size();
    }//function

    /**
     * @brief Copys from another Seed.
     */
    inline Seed& operator=(const Seed& rxOther)
    {
        Interval::operator=(rxOther);
        uiPosOnReference = rxOther.uiPosOnReference;
        return *this;
    }// operator
}; //class

#endif