/** 
 * @file execOnVector.cpp
 * @author Markus Schmidt
 */
#include "module/sw_gpu.h"
using namespace libMA;

// @fixme this is not supposed to be here...
class GPUReturn
{
public:
    int iMaxScore;
    std::vector<size_t> vMaxPos;
    GPUReturn(int iMaxScore, std::vector<size_t> vMaxPos)
            :
        iMaxScore(iMaxScore),
        vMaxPos(vMaxPos)
    {}// default constructor
    GPUReturn(){}

    bool operator==(const GPUReturn& rOther)
    {
        return iMaxScore == rOther.iMaxScore;
    }// operator
};// class

//from the separately compiled file...
extern std::vector<GPUReturn> cudaAlign
(   std::vector<char> &rvRefSeq, // reference sequence
	std::vector<std::vector<char>> &rvQuerySeqs // vector of query sequences
);

ContainerVector SW_GPU::getInputType() const
{
    return ContainerVector{
                std::shared_ptr<ContainerVector>(new ContainerVector(std::make_shared<NucSeq>())),
                std::make_shared<NucSeq>(),
            };
}//function

std::shared_ptr<Container> SW_GPU::getOutputType() const
{
    return nullptr;
}//function


std::shared_ptr<Container> SW_GPU::execute(std::shared_ptr<ContainerVector> vpInput)
{
    //const std::shared_ptr<ContainerVector>& pQueries = std::static_pointer_cast<ContainerVector>((*vpInput)[0]);
    //const std::shared_ptr<NucSeq>& pRef = std::static_pointer_cast<NucSeq>((*vpInput)[1]);

    //AUFRUF
    return nullptr;
}//function


std::vector<GPUReturn> testGPUSW(std::shared_ptr<ContainerVector> pQueries, std::shared_ptr<NucSeq> b)
{

    std::vector<std::vector<char>> vQueries;
    for( const auto& pQuery : *pQueries )
    {
        const auto& a = std::static_pointer_cast<NucSeq>(pQuery);
        std::vector<char> vQuery( a->pGetSequenceRef(), a->pGetSequenceRef() + a->length() );
        vQueries.push_back(vQuery);
    }// for
    std::vector<char> vRef( b->pGetSequenceRef(), b->pGetSequenceRef() + b->length() );

    auto vResults = cudaAlign
    (
        vRef,
        vQueries
    );

    return vResults;
}

void exportSW_GPU()
{
    boost::python::def("testGPUSW", &testGPUSW);

    //export the Tail class
    boost::python::class_<GPUReturn>("GPUReturn")
        .def_readwrite("vMaxPos", &GPUReturn::vMaxPos)
        .def_readwrite("iMaxScore", &GPUReturn::iMaxScore);

        
    boost::python::class_<std::vector<GPUReturn>>("GPUReturnList")
        .def(boost::python::vector_indexing_suite<
            std::vector<GPUReturn>,
            /*
            *    true = noproxy this means that the content of the vector is already exposed by
            *    boost python. 
            *    if this is kept as false, Container would be exposed a second time.
            *    the two Containers would be different and not inter castable.
            */
            true
        >());
        
    boost::python::class_<std::vector<size_t>>("GPUReturnMaxPosList")
        .def(boost::python::vector_indexing_suite<
            std::vector<size_t>,
            /*
            *    true = noproxy this means that the content of the vector is already exposed by
            *    boost python. 
            *    if this is kept as false, Container would be exposed a second time.
            *    the two Containers would be different and not inter castable.
            */
            true
        >());


    //export the Tail class
    //boost::python::class_<
    //        SW_GPU, 
    //        boost::python::bases<Module>,
    //        std::shared_ptr<SW_GPU>
    //    >(
    //    "SW_GPU"
    //);
    //boost::python::implicitly_convertible< 
    //    std::shared_ptr<SW_GPU>,
    //    std::shared_ptr<Module> 
    //>();
}//function