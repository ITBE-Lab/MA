#include "module/execOnVector.h"
using namespace libMA;

ContainerVector ExecOnVec::getInputType() const
{
    std::shared_ptr<ContainerVector> pVec(new ContainerVector(pModule->getInputType()[0]));

    ContainerVector vRet = pModule->getInputType();
    vRet[0] = pVec;

    return vRet;
}//function

std::shared_ptr<Container> ExecOnVec::getOutputType() const
{
    return std::shared_ptr<Container>(new ContainerVector(pModule->getOutputType()));
}//function

std::shared_ptr<Container> ExecOnVec::execute(std::shared_ptr<ContainerVector> vpInput)
{
    // vp input is organized to following way: first a vector then single elements.
    // each element of the first vector shall be executed with all following elements
    std::shared_ptr<ContainerVector> pInVec = std::shared_ptr<ContainerVector>(
        std::static_pointer_cast<ContainerVector>((*vpInput)[0]));

    DEBUG_2(
        std::cout << "executing on: " << pInVec->size() << std::endl;
    )

    std::shared_ptr<ContainerVector> pResults = std::shared_ptr<ContainerVector>(
            new ContainerVector(pModule->getOutputType(), pInVec->size())
        );

    if(pInVec->empty())
        return pResults;
    {
        ThreadPool xPool( NUM_THREADS_ALIGNER );
        for(unsigned int i = 0; i < pResults->size(); i++)
        {
            (*pResults)[i] = std::shared_ptr<Container>();
            xPool.enqueue(
                [&pResults, &vpInput]
                (
                    size_t, 
                    std::shared_ptr<ContainerVector> pInVec,
                    std::shared_ptr<Module> pModule,
                    unsigned int i
                )
                {
                    // create a input vector and add the first element to it
                    // the appropriate element of the vector that is given as first input
                    std::shared_ptr<ContainerVector> vInput(new ContainerVector { (*pInVec)[i] });
                    // add all following elements
                    for(unsigned int j = 1; j < vpInput->size(); j++)
                        vInput->push_back((*vpInput)[j]);
                    try
                    {
                        (*pResults)[i] = pModule->execute(vInput);
                    }
                    catch(NullPointerException e) 
                    {
                        std::cerr << e.what() << std::endl;
                    }
                    catch(std::exception e) 
                    {
                        std::cerr << e.what() << std::endl;
                    }
                    catch(...) 
                    {
                        std::cerr << "unknown exception when executing" << std::endl;
                    }
                    if((*pResults)[i] == nullptr)
                        std::cerr << "linesweep deleviered nullpointer as result" << std::endl;
                    if((*pResults)[i] == nullptr)
                        throw NullPointerException("linesweep deleviered nullpointer as result");
                },//lambda
                pInVec, pModule, i
            );
        }//for
    }//scope xPool

    if(sort)
    {
        std::sort(
            pResults->begin(), pResults->end(),
            []
            (std::shared_ptr<Container> a, std::shared_ptr<Container> b)
            {
                return a->smaller(b);
            }//lambda
        );//sort function call
        assert(pResults->size() <= 1 || !pResults->back()->smaller(pResults->front()));
    }//if

    if(nMany != 0 && pResults->size() > nMany)
    {
        auto end = pResults->end();
        for(unsigned int i = 0; i < nMany; i++)
            end--;
        pResults->erase(pResults->begin(), end);
        assert(pResults->size() == nMany);
    }//if

    return pResults;
}//function

ContainerVector Tail::getInputType() const
{
    return ContainerVector{
                std::shared_ptr<ContainerVector>(new ContainerVector(type)),
            };
}//function

std::shared_ptr<Container> Tail::getOutputType() const
{
    return type;
}//function


std::shared_ptr<Container> Tail::execute(std::shared_ptr<ContainerVector> vpInput)
{
    std::shared_ptr<ContainerVector> pInVec = std::shared_ptr<ContainerVector>(
        std::static_pointer_cast<ContainerVector>((*vpInput)[0]));

    std::shared_ptr<Container> pRet(nullptr);
    if(pInVec->empty())
        return type;

    /*
     * we need to transfer the ownership of the returned element to a new pointer,
     * since python might delete the input after the function call returned.
     * Since we do not want to copy anything we replace the element.
     *
     * This should be stated in the description of Tail... 
     * sth like: invalidates the input vector
     */
    pRet.swap(pInVec->back());
    pInVec->pop_back();

    return pRet;
}//function

void exportExecOnVector()
{
    //export the ExecOnVec class
    boost::python::class_<
            ExecOnVec, 
            boost::python::bases<Module>,
            std::shared_ptr<ExecOnVec>
        >(
        "ExecOnVec",
        boost::python::init<std::shared_ptr<Module>, bool, unsigned int>()
            [boost::python::with_custodian_and_ward_postcall<0,1>()]
    )
        .def(boost::python::init<
            std::shared_ptr<Module>, bool>()
            [boost::python::with_custodian_and_ward_postcall<0,1>()])
        .def(boost::python::init<std::shared_ptr<Module>>()
            [boost::python::with_custodian_and_ward_postcall<0,1>()])
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<ExecOnVec>,
        std::shared_ptr<Module> 
    >();

    //export the Tail class
    boost::python::class_<
            Tail, 
            boost::python::bases<Module>,
            std::shared_ptr<Tail>
        >(
        "Tail",
        boost::python::init<std::shared_ptr<Container>>()[
            boost::python::with_custodian_and_ward_postcall<0,1>()
        ]
    );
    boost::python::implicitly_convertible< 
        std::shared_ptr<Tail>,
        std::shared_ptr<Module> 
    >();
}//function