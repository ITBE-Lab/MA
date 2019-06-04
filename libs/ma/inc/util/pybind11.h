#ifdef WITH_PYTHON

#ifndef PYBIND11_EXT_H
#define PYBIND11_EXT_H

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
//#include <pybind11/stl.h> <- @note NEVER INCLUDE THIS HERE IT LEADS TO SEGFAULTS WITHIN PYBIND

NAMESPACE_BEGIN(PYBIND11_NAMESPACE)


// @todo request this feature for pybind lib
// std::vector
//
template <typename Vector, typename super_type, typename holder_type = std::unique_ptr<Vector>, typename... Args>
class_<Vector, holder_type> bind_vector_ext(handle scope, std::string const &name, Args&&... args) {
    using Class_ = class_<Vector, super_type, holder_type>;

    // If the value_type is unregistered (e.g. a converting type) or is itself registered
    // module-local then make the vector binding module-local as well:
    using vtype = typename Vector::value_type;
    auto vtype_info = detail::get_type_info(typeid(vtype));
    bool local = !vtype_info || vtype_info->module_local;

    Class_ cl(scope, name.c_str(), pybind11::module_local(local), std::forward<Args>(args)...);

    // Declare the buffer interface if a buffer_protocol() is passed in
    detail::vector_buffer<Vector, Class_, Args...>(cl);

    cl.def(init<>());

    // Register copy constructor (if possible)
    detail::vector_if_copy_constructible<Vector, Class_>(cl);

    // Register comparison-related operators and functions (if possible)
    detail::vector_if_equal_operator<Vector, Class_>(cl);

    // Register stream insertion operator (if possible)
    detail::vector_if_insertion_operator<Vector, Class_>(cl, name);

    // Modifiers require copyable vector value type
    detail::vector_modifiers<Vector, Class_>(cl);

    // Accessor and iterator; return by value if copyable, otherwise we return by ref + keep-alive
    detail::vector_accessor<Vector, Class_>(cl);

    cl.def("__bool__",
        [](const Vector &v) -> bool {
            return !v.empty();
        },
        "Check whether the list is nonempty"
    );

    cl.def("__len__", &Vector::size);

    return cl;
}

NAMESPACE_END(PYBIND11_NAMESPACE)


#endif
#endif
