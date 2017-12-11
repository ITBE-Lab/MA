#include "module/module.h"
using namespace libLAuS;
std::mutex xPython;

bool libLAuS::typeCheck(
        std::shared_ptr<Container> pData, 
        std::shared_ptr<Container> pExpected
    )
{
    if(pExpected->getType()->canCast(pData->getType()))
        return true;
    std::cerr << "Types did not match. Got: " << pData->getType()->getTypeName();
    std::cerr << " but expected: " << pExpected->getType()->getTypeName() << std::endl;
    return false;
}//function

bool libLAuS::typeCheck(
        ContainerVector vData, 
        ContainerVector vExpected
    )
{
    if(vData.size() != vExpected.size())
    {
        std::cerr << "Types of vectors did not match. " << std::endl;
        std::cerr << "Expected " << vExpected.size() << " elements but got " << vData.size();
        std::cerr << std::endl;
        std::cerr << "Expected: [";
        for(unsigned int i = 0; i < vExpected.size(); i++)
            std::cerr << vExpected[i]->getType()->getTypeName() << ", ";
        std::cerr << "] but got: [";
        for(unsigned int i = 0; i < vData.size(); i++)
            std::cerr << vData[i]->getType()->getTypeName() << ", ";
        std::cerr << "]" << std::endl;
        return false;
    }//for
    for(unsigned int i = 0; i < vData.size(); i++)
        if(!typeCheck(vData[i], vExpected[i]))
        {
            std::cerr << "Types of vectors did not match.";
            std::cerr << std::endl;
            std::cerr << "Element " << i << "s type is incorrect." << std::endl;
            std::cerr << "Expected: [";
            for(unsigned int i = 0; i < vExpected.size(); i++)
                std::cerr << vExpected[i]->getType()->getTypeName() << ", ";
            std::cerr << "] but got: [";
            for(unsigned int i = 0; i < vData.size(); i++)
                std::cerr << vData[i]->getType()->getTypeName() << ", ";
            std::cerr << "]" << std::endl;
            return false;
        }//if
    return true;
}//function

std::shared_ptr<Pledge> Module::promiseMe(
        std::shared_ptr<Module> pThis, 
        std::vector<std::shared_ptr<Pledge>> vInput
    )
{
    ContainerVector vCastInput(vInput.begin(), vInput.end());
    if(!typeCheck(vCastInput, pThis->getInputType()))
    {
        std::cerr << "promise of module " << pThis->getName() << " had the wrong type" << std::endl;
        throw new ModuleIO_Exception("Input type and expected input type did not match.");
    }//if
    return Pledge::makePledge(pThis, vInput);
}//function

void exportModule()
{
    //module is an abstract class and should never be initialized
    boost::python::class_<Module>(
            "Module",
            boost::python::no_init
        )
        .def(
                "execute",
                &Module::pyExecute
            )
        .def(
                "get_input_type",
                &Module::getInputType
            )
        .def(
                "get_output_type",
                &Module::getOutputType
            )
        .def(
                "get_name",
                &Module::getName
            )
        .def(
                "promise_me",
                &Module::promiseMe,
                //boost::python::with_custodian_and_ward_postcall<0,1,
                boost::python::with_custodian_and_ward_postcall<0,2>()
                //>()
            )
    ;
    boost::python::class_<
            Pledge, 
            boost::noncopyable,
            boost::python::bases<Container>, 
            std::shared_ptr<Pledge>
        >(
            "Pledge",
            boost::python::init<std::shared_ptr<Container>>()
        )
            .def(
                    "make_pledge",
                    &Pledge::makePyPledge,
                    boost::python::with_custodian_and_ward_postcall<0,1,
                    boost::python::with_custodian_and_ward_postcall<0,3>
                    >()
                )
            .staticmethod("make_pledge")
            .def(
                    "set",
                    &Pledge::set
                )
            .def(
                    "get",
                    &Pledge::get
                    //THIS WARD CAUSES THE MEMORY LEAKS
                    //,boost::python::with_custodian_and_ward_postcall<1,0>()
                )
            .def(
                    "simultaneous_get",
                    &Pledge::simultaneousGet
                )
            .staticmethod("simultaneous_get")
            .def(
                    "get_pledger",
                    &Pledge::getPledger
                )
            .def_readwrite("exec_time", &Pledge::execTime)
        ;

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
                                            std::shared_ptr<Pledge>,
                                            std::shared_ptr<Container> 
                                        >();

    iterable_converter()
        .from_python<std::vector<std::shared_ptr<Pledge>>>();
}