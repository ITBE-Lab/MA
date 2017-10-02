#include "module.h"

unsigned int Pledge::iCurrentCallNum = 0;


bool typeCheck(
        ContainerType xData, 
        ContainerType xExpected
    )
{
    if(xExpected == ContainerType::any)
        return true;
    if(xExpected == ContainerType::unknown)
        return false;
    return xData == xExpected;
}//function

bool typeCheck(
        std::shared_ptr<Container> pData, 
        ContainerType xExpected
    )
{
    if(xExpected == ContainerType::nothing && pData == nullptr)
        return true;
    return typeCheck(pData->getType(), xExpected);
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

std::shared_ptr<Pledge> Module::promiseMe(std::vector<std::shared_ptr<Pledge>> vInput)
{
    std::vector<std::shared_ptr<Container>> vCastInput(vInput.begin(), vInput.end());
    if(!typeCheck(vCastInput, getInputType()))
        throw new ModuleIO_Exception("Input type and expected input type did not match.");
    return std::shared_ptr<Pledge>(new Pledge(this, getOutputType(), vInput));
}//function

void exportModule()
{
    //module is an abstract class and should never be initialized
	boost::python::class_<Module>(
            "Module__",
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
        .def(
                "promise_me",
                &Module::promiseMe,
                boost::python::with_custodian_and_ward_postcall<0,1>(),
                "arg1: self\n"
                "arg2: a list of pledges that are required as inputs\n"
                "returns: a promise to create a container of the type getOutputType\n"
                "\n"
                "The module will only hold the pledge if all inputs are delivered."
            )
    ;
    boost::python::class_<
            Pledge, 
            boost::python::bases<Container>, 
            std::shared_ptr<Pledge>
        >(
            "Pledge",
            "Represents the pledge to deliver some container.\n"
            "Content may be provided by a Module (use promiseMe function) or by calling set.\n",
            boost::python::init<
                    boost::python::object, 
                    ContainerType, 
                    std::vector<std::shared_ptr<Pledge>>
                >
            (
                "arg1: self\n"
                "arg2: py_pledger\n"
                "arg3: type\n"
                "arg4: vPredecessors\n"
            )
            [boost::python::with_custodian_and_ward_postcall<0,1>()]
        )
            .def(boost::python::init<ContainerType>(
                "arg1: self\n"
                "arg2: Type of the promised container\n"
            ))
            .def(
                    "set",
                    &Pledge::set,
                    "arg1: self\n"
                    "arg2: the container that shall be stored\n"
                    "returns: nil\n"
                )
            .def(
                    "next",
                    &Pledge::next,
                    "arg1: self\n"
                    "returns: the pledged container\n"
                    "/n"
                    "This call will trigger all containers to invalidate their content.\n"
                    "Then every part of the computational graph that is needed to fullfill " 
                    "this pledge will be computed.\n"
                )
            .def(
                    "get",
                    &Pledge::get,
                    "arg1: self\n"
                    "returns: the pledged container\n"
                    "/n"
                    "Every part of the computational graph that is needed to fullfill " 
                    "this pledge will be computed.\n"
                )
        ;

	//tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
                                            std::shared_ptr<Pledge>,
                                            std::shared_ptr<Container> 
                                        >();

    iterable_converter()
        .from_python<std::vector<std::shared_ptr<Pledge>>>();
}