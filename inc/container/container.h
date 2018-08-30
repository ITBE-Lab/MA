/**
 * @file container.h
 * @brief Implements Container.
 * @author Markus Schmidt
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include "util/support.h"
#include "util/exception.h"
#include "util/iterableConverter.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <memory>
#ifdef WITH_PYTHON
    #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif
/// @endcond

namespace libMA
{
    class Pledge;

    /**
     * @defgroup container
     * @brief All classes containing data elements.
     * @details
     * All classes that store data should inherit from Container.
     * These classes are then added to this Group.
     */

    /**
     * @brief Abstract class intended to hold data objects used by @ref Module "modules".
     * @details
     * All classes containing data should inherit from this class.
     * @ingroup container
     */
    class Container
    {
    public:
        /** 
        * @returns the type of the container as int.
        * @brief Used by @ref Module "module" for type checking its Inputs.
        */
        volatile bool canCast(std::shared_ptr<Container> c) const
        {
            return true;
        }//function

        virtual std::string getTypeName() const
        {
            return "Container";
        }//function

        virtual std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new Container());
        }//function

        /**
         * @details
         * used by ExecOnVec to sort containers in order
         */
        virtual bool larger(const std::shared_ptr<Container> other) const
        {
            return false;
        }// operator
    };//class


    /**
     * @brief Container that shall represent Null
     * @details
     * Should be used if a module does not require input or produces no ouput.
     * @ingroup container
     */
    class Nil : public Container
    {
    public:
        static const std::shared_ptr<Nil> pEoFContainer;

        Nil()
        {}//default constructor

        bool canCast(std::shared_ptr<Container> c) const
        {
            return false;
        }//function

        std::string getTypeName() const
        {
            return "Nil";
        }//function

        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new Nil());
        }//function
    };//class

    /**
     * @brief A vector of containers, that itself is a Container
     * @details
     * This can be used to define a vector as input and output of Modules.
     * @note:
     * @todo this requires type checking! While the class holds a dummy container as type 
     * it does not yet guarantee that all added containers are of the defined type.
     * Could be solved using a template but then how should python deal with this class
     * @ingroup container
     */
    class ContainerVector: public Container
    {
    public:
        typedef std::shared_ptr<Container> TP_PTR_CONT;
        typedef std::vector<TP_PTR_CONT> TP_VEC;
        TP_VEC vContent;

        typedef typename TP_VEC::value_type value_type;
        typedef typename TP_VEC::size_type size_type;
        typedef typename TP_VEC::difference_type difference_type;
        typedef typename TP_VEC::iterator iterator;


        std::shared_ptr<Container> contentType;

        ContainerVector(std::initializer_list<std::shared_ptr<Container>> init) :
            vContent(init)
        {
            contentType = front()->getType();
        }//initializer list constructor

        template< class InputIt >
        ContainerVector(InputIt xBegin, InputIt xEnd) :
            vContent(xBegin, xEnd)
        {}//iterator constructor

        ContainerVector() :
            contentType(new Container())
        {}//container vector

        ContainerVector(std::shared_ptr<Container> contentType) :
            vContent(),
            contentType(contentType)
        {}//container vector

        ContainerVector(const std::shared_ptr<ContainerVector> pOther) :
            vContent(pOther->vContent),
            contentType(pOther->contentType)
        {}//container vector

        ContainerVector(std::shared_ptr<std::vector<std::shared_ptr<Container>>> pContent) :
            vContent(*pContent),
            contentType(pContent->front()->getType())
        {}//container vector

        ContainerVector(std::shared_ptr<Container> contentType, size_t numElements) :
            vContent(numElements),
            contentType(contentType)
        {}//container vector

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            std::shared_ptr<ContainerVector> casted = std::dynamic_pointer_cast<ContainerVector>(c);
            if(casted == nullptr)
                return false;
            return casted->contentType->getType()->canCast(contentType->getType());
        }//function

        //overload
        std::string getTypeName() const
        {
            return "ContainerVector(" + contentType->getTypeName() + ")";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            assert(contentType != nullptr);
            return std::shared_ptr<Container>(new ContainerVector(contentType));
        }//function

        std::vector<std::shared_ptr<Container>> get()
        {
            return vContent;
        }//function

        // setter
        inline TP_PTR_CONT & operator[] (size_t uiI)
        {
            return vContent[uiI];
        }// operator

        // getter
        inline const TP_PTR_CONT & operator[] (size_t uiI) const
        {
            return vContent[uiI];
        }// operator

        inline void push_back( const TP_PTR_CONT& value )
        {
            vContent.push_back( value );
        }// method

        inline void pop_back( void )
        {
            vContent.pop_back( );
        }// method
        
        template< class... Args >
        inline void emplace_back( Args&&... args )
        {
            vContent.emplace_back( args... );
        }// method
        
        inline size_t size( void ) const
        {
            return vContent.size();
        }// method
        
        inline bool empty( void ) const
        {
            return vContent.empty();
        }// method
        
        inline TP_PTR_CONT & front( void )
        {
            return vContent.front();
        }// method
        
        inline TP_PTR_CONT & back( void )
        {
            return vContent.back();
        }// method
        
        inline const TP_PTR_CONT & front( void ) const
        {
            return vContent.front();
        }// method
        
        inline const TP_PTR_CONT & back( void ) const
        {
            return vContent.back();
        }// method
        
        inline TP_VEC::iterator begin( void ) noexcept
        {
            return vContent.begin();
        }// method
        
        inline TP_VEC::iterator end( void ) noexcept
        {
            return vContent.end();
        }// method
        
        inline TP_VEC::const_iterator begin( void ) const noexcept
        {
            return vContent.begin();
        }// method
        
        inline TP_VEC::const_iterator end( void ) const noexcept
        {
            return vContent.end();
        }// method

        inline void erase( TP_VEC::iterator pos )
        {
            vContent.erase( pos );
        }// method

        inline void erase( 
                TP_VEC::iterator first,
                TP_VEC::iterator last 
            )
        {
            vContent.erase( first, last );
        }// method

        inline TP_VEC::iterator insert( TP_VEC::const_iterator pos, const TP_PTR_CONT& value )
        {
            return vContent.insert( pos, value );
        }// method
    
        template< class InputIt >
        inline TP_VEC::iterator insert( TP_VEC::const_iterator pos, InputIt first, InputIt last )
        {
            return vContent.insert( pos, first, last );
        }// method

    };//class
}//namespace libMA

#ifdef WITH_PYTHON
/** 
 * @brief Function to export Container to boost python.
 * @ingroup export
 */
void exportContainer();
#endif

#endif