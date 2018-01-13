#include "module/pipe.h"

using namespace libMABS;

std::mutex xLock;

ContainerVector Pipe::getInputType() const
{
    ContainerVector xRet;
    for(unsigned int i=0; i< uiN; i++)
        xRet.push_back(std::shared_ptr<Container>(new ContainerVector(aIn[i]->getType())));
    return xRet;
}//function

std::shared_ptr<Container> Pipe::getOutputType() const
{
    return std::shared_ptr<Container>(new ContainerVector(pOut->getType()));
}//function

std::shared_ptr<Container> Pipe::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<ContainerVector> pFirst = std::static_pointer_cast<ContainerVector>((*vpInput)[0]);
    std::shared_ptr<ContainerVector> pRet(new ContainerVector(pOut->getType()));
    while(true)
    {
        {
            std::lock_guard<std::mutex> xGuard(xLock);
            if(pFirst->empty())
                break;
            for(unsigned int i=0; i< uiN; i++)
            {
                aIn[i]->set( std::static_pointer_cast<ContainerVector>((*vpInput)[i])->back() );
                std::static_pointer_cast<ContainerVector>((*vpInput)[i])->pop_back();
            }//for
        }//end of scope xGuard lock
        pRet->push_back(pOut->get());
    }//while
    return pRet;
}//function

