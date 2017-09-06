#include "container.h"

std::vector<std::string> vContainerTypeNames =
{
    "fM_index",
    "nucSeq",
    "packedNucSeq",
    "segmentList",
    "segment",
    "stripOfConsideration",
    "stripOfConsiderationList",
    "vector",
    "unknown",
    "nothing",
    "any",
};//array

std::string Container::getTypeInfo()
{
    if(vContainerTypeNames.size() > getType())
        return std::string(vContainerTypeNames[getType()]);
    return std::string("unknown");
}//function

//Taken from : https://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python
/// @brief Type that allows for registration of conversions from
///        python iterable types.
struct iterable_converter
{
    /// @note Registers converter from a python interable type to the
    ///       provided type.
    template <typename Container>
    iterable_converter&
    from_python()
    {
        boost::python::converter::registry::push_back(
        &iterable_converter::convertible,
        &iterable_converter::construct<Container>,
        boost::python::type_id<Container>());

        // Support chaining.
        return *this;
    }

    /// @brief Check if PyObject is iterable.
    static void* convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    /// @brief Convert iterable PyObject to C++ container type.
    ///
    /// Container Concept requirements:
    ///
    ///   * Container::value_type is CopyConstructable.
    ///   * Container can be constructed and populated with two iterators.
    ///     I.e. Container(begin, end)
    template <typename Container>
    static void construct(
            PyObject* object,
            boost::python::converter::rvalue_from_python_stage1_data* data
        )
    {
        namespace python = boost::python;
        // Object is a borrowed reference, so create a handle indicting it is
        // borrowed for proper reference counting.
        python::handle<> handle(python::borrowed(object));

        // Obtain a handle to the memory block that the converter has allocated
        // for the C++ type.
        typedef python::converter::rvalue_from_python_storage<Container>
                                                                    storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        typedef python::stl_input_iterator<typename Container::value_type>
                                                                        iterator;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        new (storage) Container(
        iterator(python::object(handle)), // begin
        iterator());                      // end
        data->convertible = storage;
    }
};

void exportContainer()
{
    //contianer is an abstract class and should never be initialized
	boost::python::class_<Container, std::shared_ptr<Container>>(
            "Container", 
            "class: Container abstract\n"
            "   Any class holding data should inherit Container.\n",
            boost::python::no_init
        )
        .def(
                "get_type_info", 
                &Container::getTypeInfo,
                "method: get_type_info()\n"
                "   returns: an enum describing the type of data stored.\n"
                "\n"
                "   Used to check weather a module can work with the given containers.\n"
            );

    boost::python::register_ptr_to_python< std::shared_ptr<Container> >();

    //make vectors of container-pointers a thing
    iterable_converter()
        .from_python<std::vector<std::shared_ptr<Container>>>()
        .from_python<std::vector<ContainerType>>();

    //export the containertype enum
    boost::python::enum_<ContainerType>("ContainerType")
        .value("unknown", ContainerType::unknown)
        .value("vector", ContainerType::vector)
        .value("FM_index", ContainerType::fM_index)
        .value("nucSeq", ContainerType::nucSeq)
        .value("segmentList", ContainerType::segmentList)
        .value("segment", ContainerType::segment)
        .value("packedNucSeq", ContainerType::packedNucSeq)
        .value("nothing", ContainerType::nothing);
}//function