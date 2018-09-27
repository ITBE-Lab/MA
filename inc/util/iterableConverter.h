/**
 * @file iterableConverter.h
 * @brief implements IterableConverter
 * @note: Taken from :
 * https://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python
 */

#ifdef WITH_PYTHON

#ifndef ITERABLE_CONVERTER_H
#define ITERABLE_CONVERTER_H

/**
 * @brief Type that allows for registration of conversions from python iterable types.
 */
class IterableConverter
{
  private:
    /// @brief Check if PyObject is iterable.
    static void *convertible( PyObject *object )
    {
        return PyObject_GetIter( object ) ? object : NULL;
    } // function

    /// @brief Convert iterable PyObject to C++ container type.
    ///
    /// C Concept requirements:
    ///
    ///   * C::value_type is CopyConstructable.
    ///   * C can be constructed and populated with two iterators.
    ///     I.e. C(begin, end)
    template <typename C>
    static void construct( PyObject *object,
                           boost::python::converter::rvalue_from_python_stage1_data *data )
    {
        // Object is a borrowed reference, so create a handle indicting it is
        // borrowed for proper reference counting.
        boost::python::handle<> handle( boost::python::borrowed( object ) );

        // Obtain a handle to the memory block that the converter has allocated
        // for the C++ type.
        typedef boost::python::converter::rvalue_from_python_storage<C> storage_type;
        void *storage = reinterpret_cast<storage_type *>( data )->storage.bytes;

        typedef boost::python::stl_input_iterator<typename C::value_type> iterator;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.  The C++
        // c is populated by passing the begin and end iterators of
        // the python object to the c's constructor.
        new( storage ) C( iterator( boost::python::object( handle ) ), // begin
                          iterator( ) ); // end
        data->convertible = storage;
    } // function
  public:
    /// @brief Registers converter from a python interable type to the
    ///       provided type.
    template <typename C> IterableConverter &from_python( )
    {
        boost::python::converter::registry::push_back( &IterableConverter::convertible,
                                                       &IterableConverter::construct<C>,
                                                       boost::python::type_id<C>( ) );

        // Support chaining.
        return *this;
    } // function
}; // class

#endif
#endif
