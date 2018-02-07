#include "container/alignment.h"
using namespace libMA;

//The scores from needlemanWunsch.cpp
extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;

void Alignment::append(MatchType type, nucSeqIndex size)
{
    //adjust the score of the alignment
    if(type == MatchType::seed || type == MatchType::match)
        iScore += iMatch * size;
    else if(type == MatchType::missmatch)
        //iMissMatch is a penalty not a score
        iScore -= iMissMatch * size;
    else if(type == MatchType::insertion || type == MatchType::deletion)
    {
        //iGap & iExtend is a penalty not a score
        nucSeqIndex s = size;
        if(length() == 0 || (at(length()-1) != type) )
        {
            iScore -= iGap;
            //s--;
        }//if
        iScore -= iExtend * s;
    }//if
    /*
    * we are storing in a compressed format 
    * since it actually makes quite a lot of things easier
    * same thing here: just check weather the last symbol is the same as the inserted one
    *      if so just add the amount of new symbols
    *      else make a new entry with the correct amount of symbols
    */
    if(data.size() != 0 && std::get<0>(data.back()) == type)
        std::get<1>(data.back()) += size;
    else
        data.push_back(std::make_tuple(type, size));
    uiLength += size;
}//function

void EXPORTED Alignment::makeLocal()
{
    if(uiLength == 0)
        return;
    std::vector<int> vScores;
    int iMaxScore = 0;
    unsigned int iMaxStart = 0;
    unsigned int iMaxEnd = data.size();
    unsigned int iLastStart = 0;
    int iScoreCurr = 0;

    for(unsigned int index = 0; index < data.size(); index++)
    {
        switch (std::get<0>(data[index]))
        {
            case MatchType::deletion :
            case MatchType::insertion :
                iScoreCurr -= iGap;
                iScoreCurr -= iExtend * std::get<1>(data[index]);
                break;
            case MatchType::missmatch :
                iScoreCurr -= iMissMatch * std::get<1>(data[index]);
                break;
            case MatchType::match :
            case MatchType::seed :
                iScoreCurr += iMatch * std::get<1>(data[index]);
                break;
        }//switch
        if(iScoreCurr <= 0)
        {
            iScoreCurr = 0;
            iLastStart = index+1;
        }//if
        DEBUG_2(
            std::cout << std::get<0>(data[index]) << "," << std::get<1>(data[index]) <<
            " (" << iScoreCurr << ") | ";
        )
        if(iScoreCurr >= iMaxScore)
        {
            iMaxScore = iScoreCurr;
            iMaxStart = iLastStart;
            iMaxEnd = index+1;
        }//if
    }//for
    DEBUG_2(
        std::cout << std::endl;
        std::cout << iMaxStart << " " << iMaxEnd << std::endl;
    )
    //adjust the begin/end on ref/query
    if(iMaxStart <= iMaxEnd)
    {
        for(unsigned int index = 0; index < iMaxStart; index++)
            switch (std::get<0>(data[index]))
            {
                case MatchType::deletion :
                    uiBeginOnRef += std::get<1>(data[index]);
                    break;
                case MatchType::insertion :
                    uiBeginOnQuery += std::get<1>(data[index]);
                    break;
                default :
                    uiBeginOnRef += std::get<1>(data[index]);
                    uiBeginOnQuery += std::get<1>(data[index]);
                    break;
            }//switch
        for(unsigned int index = iMaxEnd; index < data.size(); index++)
            switch (std::get<0>(data[index]))
            {
                case MatchType::deletion :
                    uiEndOnRef -= std::get<1>(data[index]);
                    break;
                case MatchType::insertion :
                    uiEndOnQuery -= std::get<1>(data[index]);
                    break;
                default :
                    uiEndOnRef -= std::get<1>(data[index]);
                    uiEndOnQuery -= std::get<1>(data[index]);
                    break;
            }//switch
    }//if
    //erase everything before and after
    if(iMaxEnd < data.size())
        data.erase(data.begin()+iMaxEnd, data.end());
    if(iMaxStart > 0)
        data.erase(data.begin(), data.begin()+iMaxStart);
    //adjust score accordingly
    iScore = iMaxScore;
    //adjust length variable accordingly
    uiLength = 0;
    for(unsigned int index = 0; index < data.size(); index++)
        uiLength += std::get<1>(data[index]);
    DEBUG(
        if(uiEndOnRef < uiBeginOnRef)
        {
            std::cout << "---" << std::endl;
            for(unsigned int index = 0; index < data.size(); index++)
                std::cout << std::get<0>(data[index]) << "," << std::get<1>(data[index]) << std::endl;
                exit(0);
        }
        if(reCalcScore() != iScore)
            std::cerr << "WARNING set wrong score or removed wrong elements in makeLocal" 
                << std::endl;
    )
}//function

int Alignment::reCalcScore() const
{
    int iScore = 0;
    for(unsigned int index = 0; index < data.size(); index++)
        switch (std::get<0>(data[index]))
        {
            case MatchType::deletion :
            case MatchType::insertion :
                iScore -= iGap;
                iScore -= iExtend * std::get<1>(data[index]);
                break;
            case MatchType::missmatch :
                iScore -= iMissMatch * std::get<1>(data[index]);
                break;
            case MatchType::match :
            case MatchType::seed :
                iScore += iMatch * std::get<1>(data[index]);
                break;
        }//switch
    return iScore;
}//function


void exportAlignment()
{
    boost::python::class_<
            Alignment, 
            boost::python::bases<Container>, 
            std::shared_ptr<Alignment>
        >(
                "Alignment",
                "contains the final output of the aligner\n"
            )
        .def(
                "at", 
                &Alignment::at,
                "arg1: self\n"
                "arg2: index at which to look\n"
                "returns: match type at the given position\n"
            )
        .def(boost::python::init<nucSeqIndex>())
        .def(
                "__getitem__", 
                &Alignment::at,
                "arg1: self\n"
                "arg2: index at which to look\n"
                "returns: match type at the given position\n"
            )
        .def(
                "append", 
                &Alignment::append_boost1,
                "arg1: self\n"
                "arg2: the matchtype to append\n"
                "arg3: how many times shall the matchtype be appended\n"
                "returns: nil\n"
            )
        .def(
                "append", 
                &Alignment::append_boost2,
                "arg1: self\n"
                "arg2: the matchtype to append\n"
                "returns: nil\n"
            )
        .def(
                "begin_on_ref", 
                &Alignment::beginOnRef,
                "arg1: self\n"
                "returns: starting position of the alignment on the reference\n"
            )
        .def(
                "end_on_ref", 
                &Alignment::endOnRef,
                "arg1: self\n"
                "returns: ending position of the alignment on the reference\n"
            )
        .def(
                "__len__", 
                &Alignment::length,
                "arg1: self\n"
                "returns: length of the alignmen\n"
            )
        .def(
                "length", 
                &Alignment::length,
                "arg1: self\n"
                "returns: length of the alignmen\n"
            )
        .def(
                "seed_coverage", 
                &Alignment::seedCoverage,
                "arg1: self\n"
            )
        .def(
                "get_score", 
                &Alignment::score,
                "arg1: self\n"
            )
        .def(
                "num_by_seeds", 
                &Alignment::numBySeeds,
                "arg1: self\n"
            )
        .def(
                "make_local", 
                &Alignment::makeLocal,
                "arg1: self\n"
            )
        .def_readonly("stats", &Alignment::xStats)
        .def_readonly("begin_on_query", &Alignment::uiBeginOnQuery)
        .def_readonly("end_on_query", &Alignment::uiEndOnQuery)
        .def_readwrite("mapping_quality", &Alignment::fMappingQuality)
    ;

    
    //export the matchType enum
    boost::python::enum_<MatchType>("MatchType")
        .value("match", MatchType::match)
        .value("seed", MatchType::seed)
        .value("missmatch", MatchType::missmatch)
        .value("insertion", MatchType::insertion)
        .value("deletion", MatchType::deletion);

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
        std::shared_ptr<Alignment>,
        std::shared_ptr<Container> 
    >();

}