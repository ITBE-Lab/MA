#include "module.h"


bool typeCheck(
        std::shared_ptr<Container> pData, 
        ContainerType xExpected
    )
{
    if(xExpected == ContainerType::any)
        return true;
    if(xExpected == ContainerType::unknown)
        return false;
    if(xExpected == ContainerType::nothing)
        return pData == nullptr;
    return pData->getType() == xExpected;
}//function

bool typeCheck(
        std::vector<std::shared_ptr<Container>> vData, 
        std::vector<ContainerType> vExpected
    )
{
    if(vData.size() != vExpected.size())
        return false;
    for(unsigned int i = 0; i < vData.size(); i++)
        if(!typeCheck(vData[i], vExpected[i]))
            return false;
    return true;
}//function

Pledge Module::promiseMe(std::vector<std::shared_ptr<Pledge>> vInput)
{
    std::vector<std::shared_ptr<Container>> vCastInput(vInput.begin(), vInput.end());
    if(!typeCheck(vCastInput, getInputType()))
        throw new ModuleIO_Exception("Input type and expected input type did not match.");
    return Pledge(this, getOutputType(), vInput);
}//function

void exportModule()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<Module>(
            "Module",
            "Any class implementing a algorithm should inherit Module\n",
            boost::python::no_init
        )
        .def(
                "execute",
                &Module::saveExecute,
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
                &Module::getInputType,
                "arg1: self\n"
                "returns: the type of containers that is expected as input by this module\n"
            )
        .def(
                "get_output_type",
                &Module::getOutputType,
                "arg1: self\n"
                "returns: the type of containers that is expected as output from this module\n"
            )
    ;
    /*boost::python::class_<
            Pledge, 
            boost::python::bases<Container>, 
            std::shared_ptr<Pledge>
        >(
            "Pledge",
            "Represents the pledge to deliver some container.\n"
            "Content may be provided by a Module (use promiseMe function) or by calling set.\n",
            boost::python::init<ContainerType>(
                "arg1: self\n"
                "arg2: Type of the promised container\n"
            )
        )
            .def(
                    "set",
                    &Pledge::set,
                    "arg1: self\n"
                    "arg2: the container that shall be stored\n"
                    "arg3: the number of the call the container shall be stored for\n"
                )
        ;
                */
}