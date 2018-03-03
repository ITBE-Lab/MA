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
    auto pRet = std::make_shared<MinimizersVector<w,k>>();
    MinimizersVector<w,k> vCur;
    MinimizersVector<w,k> vLast;
    /*
     * we collect all minimizers
     */
    for(nucSeqIndex i = 0; i < pQuerySeq->length() - w - k; i++)
    {
        vCur.clear();
        //save the minimum minimizers in vCur (initially the first element is the minimum)
        vCur.push_back(std::make_pair(Minimizer<k>(*pQuerySeq, i), i));
        for(unsigned int j = 1; j < w; j++)
        {
            Minimizer<k> xAlt(*pQuerySeq, i+j);
            //if we found a new minimum clear vCur and save the new minimum
            if(xAlt <= vCur[0].first)
            {
                if(xAlt < vCur[0].first)
                    vCur.clear();
                vCur.push_back(std::make_pair(xAlt, i+j));
            }//if
        }//for
        //save the minimizer(s)
        for(auto& xMini : vCur)
        {
            //check if the minimizer was in vLast if so ignore it...
            bool bNew = false;
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
    };
}
std::shared_ptr<Container> MinimizersToSeeds::getOutputType() const
{
	return std::shared_ptr<Container>(new MinimizersVector<w,k>());
}


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
		;
	boost::python::implicitly_convertible< 
		std::shared_ptr<Minimizers>,
		std::shared_ptr<Module> 
	>();

}//function