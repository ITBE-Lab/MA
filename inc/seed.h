#ifndef SEED_H
#define SEED_H

#include "container.h"
#include "interval.h"

//any index on the query or reference nucleotide sequence is given in this datatype
typedef uint64_t nucSeqIndex;


class Seed: public Container, public Interval<nucSeqIndex>
{
private:
    ////the beginning of the match on the reference
    nucSeqIndex uiPosOnReference;

public:
    
    Seed(
            const nucSeqIndex uiPosOnQuery, 
            const nucSeqIndex uiLenght, 
            const nucSeqIndex uiPosOnReference
        )
            :
        Interval(uiPosOnQuery, uiLenght),
        uiPosOnReference(uiPosOnReference)
    {}//constructor
    
    nucSeqIndex start_ref() const
    {
        return uiPosOnReference;
    }//function
    
    nucSeqIndex end_ref() const
    {
        return uiPosOnReference + size();
    }//function
    
    nucSeqIndex getValue() const
    {
        return size();
    }//function

    inline Seed& operator=(const Seed& rxOther)
    {
        Interval::operator=(rxOther);
        uiPosOnReference = rxOther.uiPosOnReference;
        return *this;
    }// operator
}; //class

#endif