/**
 * @file exception.h
 * @brief Defines various exception classes that are used throughout the aligner.
 * @author Arne Kutzner
 * @author Markus Schmidt
 * @details
 * The exceptions are exposed to python.
 */

/**
 * @defgroup exception
 * @brief several classes that deal with exceptions within the cpp code.
 * @details
 * Annotated_exception provides a way for printing the specific error.
 * @{
 */

#pragma once

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <cmath>
#include <string>

#ifdef WITH_PYTHON
// Bug in Python 3.7 and 3.8 that breaks nlohmann::json
// See:(https://bugs.python.org/issue36020)
#if defined( WIN32 ) && !defined( HAVE_SNPRINTF ) && defined( _MSC_VER ) && _MSC_VER >= 1900
#define HAVE_SNPRINTF
#endif

#include <pybind11/pybind11.h> // Question AK: Isn't it better to include our pybind11.h here.
namespace py = pybind11;
#endif
/// @endcond

/**
 * @brief An annotated exception class on the foundation of std::exception.
 * @details
 * can be overloaded to provide different tpyes of exeptions,
 * where each exception object can have it's unique description.
 */
class AnnotatedException : public std::exception
{
  private:
    /* automated memory deallocation !
     */
    std::string sText;

  public:
    /**
     * @brief takes the string that shall be printed in case the exception is thrown.
     * @details
     * prepends information about the exception type to the string.
     */
    AnnotatedException( std::string sText ) : sText( sText )
    {} // constructor

    ~AnnotatedException( ) throw( )
    {} // destructor

    /**
     * @brief Information about the exception.
     * @returns instance specific information about the exception.
     */
    virtual const char* what( ) const throw( )
    {
        return sText.c_str( );
    } // method
};