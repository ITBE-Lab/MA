#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <boost/python.hpp>
#include "container.h"
#include "intervalTree.h"

class Alignment : public Container
{
public:
    enum MatchType{
        match,
        missmatch,
        insertion,
        deletion
    };//enum
private:
    std::vector<std::tuple<MatchType, nucSeqIndex>> data;
    nucSeqIndex uiLength;
    nucSeqIndex uiBeginOnRef;
    nucSeqIndex uiEndOnRef;

public:
    Alignment()
            :
        data(),
        uiLength(0),
        uiBeginOnRef(0)
    {}//constructor
    Alignment(nucSeqIndex uiBeginOnRef, nucSeqIndex uiEndOnRef)
            :
        data(),
        uiLength(0),
        uiBeginOnRef(uiBeginOnRef),
        uiEndOnRef(uiEndOnRef)
    {}//constructor


    ContainerType getType(){return ContainerType::alignment;}

    MatchType at(nucSeqIndex i)
    {
        //everything after the query is a deletion
        if(i >= uiLength)
            return MatchType::deletion;

        //the MatchType match type is stored in a compressed format -> extract it
        nucSeqIndex j = 0;
        unsigned int k = 0;
        while(k < data.size() && (j += std::get<1>(data[k++])) <= i);

        return std::get<0>(data[k-1]);
    }//function

    MatchType operator[](nucSeqIndex i)
    {
        return at(i);
    }//operator

    void append(MatchType type, nucSeqIndex size)
    {
        if(data.size() != 0 && std::get<0>(data.back()) == type)
            std::get<1>(data.back()) += size;
        else
            data.push_back(std::make_tuple(type, size));
        uiLength += size;
    }//function

    void append(MatchType type)
    {
        append(type, 1);
    }//function

    void append_boost1(MatchType type, nucSeqIndex size)
    {
        append(type, size);
    }//function

    void append_boost2(MatchType type)
    {
        append(type, 1);
    }//function

    nucSeqIndex lenght()
    {
        return uiLength;
    }//function

    nucSeqIndex beginOnRef()
    {
        return uiBeginOnRef;
    }//function

    nucSeqIndex endOnRef()
    {
        return uiEndOnRef;
    }//function
};//class

/* function to export this module to boost python */
void exportAlignment();

#endif