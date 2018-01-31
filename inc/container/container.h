/**
 * @file container.h
 * @brief Implements Container.
 * @author Markus Schmidt
 */

#ifndef CONTAINER_H
#define CONTAINER_H

#include "util/support.h"
#include <vector>
#include "util/exception.h"
#include "util/iterableConverter.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

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
        //@todo document me
        ///to remember if one of the inputs is dry
        bool bDry = false;
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
        virtual bool smaller(const std::shared_ptr<Container> other) const
        {
            return true;
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
        Nil()
        {}//default constrctor
        /** 
        * @returns the type of the container as int.
        * @brief Used by @ref Module "module" for type checking its Inputs.
        */
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
     * @ingroup container
     */
    class ContainerVector:
        public Container,
        public std::vector<std::shared_ptr<Container>>
    {
	public:
        std::shared_ptr<Container> contentType;

		ContainerVector(std::initializer_list<std::shared_ptr<Container>> init)
			:
			vector(init)
		{}//initializer list constructor

		template< class InputIt >
		ContainerVector(InputIt xBegin, InputIt xEnd
			)
			:
			vector(xBegin, xEnd)
		{}//iterator constructor

        ContainerVector()
            :
            contentType(new Container())
        {}//container vector

        ContainerVector(std::shared_ptr<Container> contentType)
                :
            vector(),
            contentType(contentType)
        {}//container vector

        ContainerVector(const std::shared_ptr<ContainerVector> pOther)
                :
            vector(*pOther),
            contentType(pOther->contentType)
        {}//container vector

        ContainerVector(std::shared_ptr<std::vector<std::shared_ptr<Container>>> pContent)
                :
            vector(*pContent),
            contentType(pContent->front()->getType())
        {}//container vector

        ContainerVector(std::shared_ptr<Container> contentType, size_t numElements)
            :
            vector(numElements),
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
            return "ContainerVector";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            assert(contentType != nullptr);
            return std::shared_ptr<Container>(new ContainerVector(contentType));
        }//function

        std::vector<std::shared_ptr<Container>> get()
        {
            return *this;
        }//function



    };//class
}//namespace libMA

/** 
 * @brief Function to export Container to boost python.
 * @ingroup export
 */
void exportContainer();
#endif