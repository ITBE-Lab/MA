#include "execOnVector.h"



std::vector<std::shared_ptr<Container>> ExecOnVec::getInputType()
{
	return std::vector<std::shared_ptr<Container>>{
			std::shared_ptr<Container>(new ContainerVector(pModule->getInputType())),
		};
}
std::shared_ptr<Container> ExecOnVec::getOutputType()
{
	return std::shared_ptr<Container>(new ContainerVector(pModule->getInputType()));
}

void exportExecOnVector()
{
    //export the LineSweepContainer class
	boost::python::class_<
            ExecOnVec, 
            boost::python::bases<CppModule>,
            std::shared_ptr<ExecOnVec>
        >(
        "ExecOnVec",
        "Uses linesweeping to remove contradicting "
        "matches within several strips of consideration.\n",
        boost::python::init<std::shared_ptr<CppModule>>()
    );
	boost::python::implicitly_convertible< 
		std::shared_ptr<ExecOnVec>,
		std::shared_ptr<CppModule> 
	>();
}//function