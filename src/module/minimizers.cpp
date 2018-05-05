/** 
 * @file minimizers.cpp
 * @author Markus Schmidt
 */
#include "module/minimizers.h"

using namespace libMA;


//#include "BWTmem.h"
#include <vector>
//#include "analysis.h"
#include <memory>
#include <atomic>
#include <chrono>
//#include "assembly.h"


ContainerVector Minimizers::getInputType() const
{
    return ContainerVector{
            //the query sequence
            std::shared_ptr<Container>(new NucSeq()),
        };
}
std::shared_ptr<Container> Minimizers::getOutputType() const
{
    return std::shared_ptr<Container>(new MinimizersVector<w,k>());
}


std::shared_ptr<Container> Minimizers::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<NucSeq> pQuerySeq = std::dynamic_pointer_cast<NucSeq>((*vpInput)[0]);
    MinimizersVector<w,k> vCur;
    MinimizersVector<w,k> vLast;
    auto pRet = std::make_shared<MinimizersVector<w,k>>();
    /*
     * we collect all (u,k)-start-minimizers for 1 <= u < w
     */
    for(nucSeqIndex u = 1; u < w; u++)
    {
        vCur.clear();
        Minimizer<k> xFirst(pQuerySeq, 0);
        //save the minimum minimizers in vCur (initially the first element is the minimum)
        vCur.push_back(std::make_pair(xFirst, 0));
        for(nucSeqIndex j = 1; j < u; j++)
        {
            Minimizer<k> xAlt(pQuerySeq, j);
            assert(!vCur.empty());

            //if we found a new minimum clear vCur and save the new minimum
            if(xAlt <= vCur[0].first)
            {
                if(xAlt < vCur[0].first)
                    vCur.clear();
                vCur.push_back(std::make_pair(xAlt, j));
            }//if
        }//for
        //save the minimizer(s)
        for(auto& xMini : vCur)
        {
            assert(xMini.second + k < pQuerySeq->length());
            //check if the minimizer was in vLast if so ignore it...
            bool bNew = true;
            for(auto& xMini2 : vLast)
            {
                if(xMini2.second == xMini.second)
                    bNew = false;
                if(xMini2.second >= xMini.second)
                    break;
            }//for
            //save the new minimizer(s)
            if(bNew)
                pRet->push_back(xMini);
        }//for
        //save the current minimizers so that we can check for duplicates in the nxt iteration
        //(this required constant time...)
        vLast.swap(vCur);
    }//for
    /*
     * we collect all (w,k)-minimizers
     */
    for(nucSeqIndex i = 0; i < pQuerySeq->length() - w - k; i++)
    {
        if(bPrint && i % 1000000 == 0)
            std::cout << i/1000000 << "/" << ( pQuerySeq->length() - w - k)/1000000 << std::endl;

        vCur.clear();
        Minimizer<k> xFirst(pQuerySeq, i);
        //save the minimum minimizers in vCur (initially the first element is the minimum)
        vCur.push_back(std::make_pair(xFirst, i));
        for(nucSeqIndex j = i+1; j < i+w; j++)
        {
            Minimizer<k> xAlt(pQuerySeq, j);
            assert(!vCur.empty());

            //if we found a new minimum clear vCur and save the new minimum
            if(xAlt <= vCur[0].first)
            {
                if(xAlt < vCur[0].first)
                    vCur.clear();
                vCur.push_back(std::make_pair(xAlt, j));
            }//if
        }//for
        //save the minimizer(s)
        for(auto& xMini : vCur)
        {
            assert(xMini.second + k < pQuerySeq->length());
            //check if the minimizer was in vLast if so ignore it...
            bool bNew = true;
            for(auto& xMini2 : vLast)
            {
                if(xMini2.second == xMini.second)
                    bNew = false;
                if(xMini2.second >= xMini.second)
                    break;
            }//for
            //save the new minimizer(s)
            if(bNew)
                pRet->push_back(xMini);
        }//for
        //save the current minimizers so that we can check for duplicates in the nxt iteration
        //(this required constant time...)
        vLast.swap(vCur);
    }//for
    /*
     * we collect all (u,k)-end-minimizers for 1 <= u < w
     */
    for(nucSeqIndex u = 1; u < w; u++)
    {
        vCur.clear();
        auto uiOffset = pQuerySeq->length() - w - k;
        Minimizer<k> xFirst(pQuerySeq, uiOffset);
        //save the minimum minimizers in vCur (initially the first element is the minimum)
        vCur.push_back(std::make_pair(xFirst, uiOffset));
        for(nucSeqIndex j = 1; j < u; j++)
        {
            Minimizer<k> xAlt(pQuerySeq, j + uiOffset);
            assert(!vCur.empty());

            //if we found a new minimum clear vCur and save the new minimum
            if(xAlt <= vCur[0].first)
            {
                if(xAlt < vCur[0].first)
                    vCur.clear();
                vCur.push_back(std::make_pair(xAlt, j + uiOffset));
            }//if
        }//for
        //save the minimizer(s)
        for(auto& xMini : vCur)
        {
            assert(xMini.second + k < pQuerySeq->length());
            //check if the minimizer was in vLast if so ignore it...
            bool bNew = true;
            for(auto& xMini2 : vLast)
            {
                if(xMini2.second == xMini.second)
                    bNew = false;
                if(xMini2.second >= xMini.second)
                    break;
            }//for
            //save the new minimizer(s)
            if(bNew)
                pRet->push_back(xMini);
        }//for
        //save the current minimizers so that we can check for duplicates in the nxt iteration
        //(this required constant time...)
        vLast.swap(vCur);
    }//for

     DEBUG(
        for(auto xPair : *pRet)
        {
            auto xCheck = Minimizer<k>(pQuerySeq, xPair.second);
            if(xPair.first != xCheck )
            {
                std::cout << k << "-Minimizer wrong; is complement:" 
                    << (xPair.first.bRevComp ? "true" : "false") << std::endl;
                exit(0);
            }//if
        }//for
    )//DEBUG


    return pRet;
}//function


ContainerVector MinimizersToSeeds::getInputType() const
{
    return ContainerVector
    {
        //the hashindex
        std::shared_ptr<Container>(new MinimizersHash<Minimizers::w,Minimizers::k>()),
        //the query sequence
        std::shared_ptr<Container>(new MinimizersVector<Minimizers::w,Minimizers::k>())
        //the reference
        DEBUG_PARAM(std::shared_ptr<Container>(new Pack()))
    };
}
std::shared_ptr<Container> MinimizersToSeeds::getOutputType() const
{
    return std::shared_ptr<Container>(new Seeds());
}

std::shared_ptr<Container> MinimizersToSeeds::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    auto pIndex = std::dynamic_pointer_cast<
        MinimizersHash<Minimizers::w,Minimizers::k>>((*vpInput)[0]);
    assert(!pIndex->vKeys.empty());
    auto pMinimizers = std::dynamic_pointer_cast<
        MinimizersVector<Minimizers::w,Minimizers::k>>((*vpInput)[1]);
    DEBUG(
        std::shared_ptr<Pack> pRefSeq = 
            std::static_pointer_cast<Pack>((*vpInput)[2]);
    )
    return pIndex->toSeeds(
            pMinimizers
            DEBUG_PARAM(pRefSeq)
        );
}//function

#ifdef WITH_PYTHON
void exportMinimizers()
{
    //export the BinarySeeding class
    boost::python::class_<
            Minimizers, 
            boost::python::bases<Module>,
            std::shared_ptr<Minimizers>
        >(
            "Minimizers"
        )
        .def_readwrite("print", &Minimizers::bPrint)
        ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<Minimizers>,
        std::shared_ptr<Module> 
    >();
    //export the BinarySeeding class
    boost::python::class_<
            MinimizersToSeeds, 
            boost::python::bases<Module>,
            std::shared_ptr<MinimizersToSeeds>
        >(
            "MinimizersToSeeds"
        )
        ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<MinimizersToSeeds>,
        std::shared_ptr<Module> 
    >();

}//function
#endif
