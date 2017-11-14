#include "cppModule.h"

std::mutex xPython;

bool typeCheck(
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

bool typeCheck(
        std::vector<std::shared_ptr<Container>> vData, 
        std::vector<std::shared_ptr<Container>> vExpected
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
        std::cerr << std::endl;
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

std::shared_ptr<Pledge> CppModule::promiseMe(
        std::shared_ptr<CppModule> pThis, 
        std::vector<std::shared_ptr<Pledge>> vInput
    )
{
    std::vector<std::shared_ptr<Container>> vCastInput(vInput.begin(), vInput.end());
    if(!typeCheck(vCastInput, pThis->getInputType()))
        throw new ModuleIO_Exception("Input type and expected input type did not match.");
    return Pledge::makePledge(pThis, vInput);
}//function

void exportModule()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<CppModule>(
            "CppModule",
            "Any class implementing a algorithm should inherit Module\n",
            boost::python::no_init
        )
        .def(
                "execute",
                &CppModule::pyExecute,
                "arg1: self\n"
                "arg2: tuple/vector of containers. "
                "Their types must match the input type of this module\n"
                "returns: a container holding the results of the "
                "computations done by this module. Must match the output type of the module.\n"
                "\n"
                "perform the task implemented by the respective instance.\n"
            )
        .def(
                "get_input_type",
                &CppModule::getInputType,
                "arg1: self\n"
                "returns: the type of containers that is expected as input by this module\n"
            )
        .def(
                "get_output_type",
                &CppModule::getOutputType,
                "arg1: self\n"
                "returns: the type of containers that is expected as output from this module\n"
            )
        .def(
                "promise_me",
                &CppModule::promiseMe,
                boost::python::with_custodian_and_ward_postcall<0,1,
                boost::python::with_custodian_and_ward_postcall<0,2>
                >(),
                "arg1: self\n"
                "arg3: a list of pledges that are required as inputs\n"
                "returns: a promise to create a container of the type getOutputType\n"
                "\n"
                "The module will only hold the pledge if all inputs are delivered."
            )
    ;
    boost::python::class_<
            Pledge, 
            boost::noncopyable,
            boost::python::bases<Container>, 
            std::shared_ptr<Pledge>
        >(
            "Pledge",
            "Represents the pledge to deliver some container.\n"
            "Content may be provided by a Module (use promiseMe function) or by calling set.\n",
            boost::python::init<std::shared_ptr<Container>>(
                "arg1: self\n"
                "arg2: Type of the promised container\n"
            )
        )
            .def(
                    "make_pledge",
                    &Pledge::makePyPledge,
                    boost::python::with_custodian_and_ward_postcall<0,1,
                    boost::python::with_custodian_and_ward_postcall<0,3>
                    >(),
                    "arg1: self\n"
                    "arg3: a list of pledges that are required as inputs\n"
                    "returns: a promise to create a container of the type getOutputType\n"
                    "\n"
                    "The module will only hold the pledge if all inputs are delivered."
                )
            .staticmethod("make_pledge")
            .def(
                    "set",
                    &Pledge::set,
                    "arg1: self\n"
                    "arg2: the container that shall be stored\n"
                    "returns: nil\n"
                )
            .def(
                    "get",
                    &Pledge::get,
                    boost::python::with_custodian_and_ward_postcall<1,0>(),
                    "arg1: self\n"
                    "returns: the pledged container\n"
                    "/n"
                    "Every part of the computational graph that is needed to fullfill " 
                    "this pledge will be computed.\n"
                )
            .def(
                    "simultaneous_get",
                    &Pledge::simultaneousGet,
                    "arg1: class\n"
                    "arg2: pledges\n"
                    "arg1: num_threads\n"
                    "/n"
                )
            .staticmethod("simultaneous_get")
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