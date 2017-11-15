/**
 * TODO:
 * @file execOnVector.h
 */

#ifndef EXEC_ON_VECTOR_H
#define EXEC_ON_VECTOR_H

#include "cppModule.h"
#include "threadPool.h"

/**
 * @brief TODO:
 */
class ExecOnVec: public CppModule
{
private:
    std::shared_ptr<CppModule> pModule;
    bool sort;
    unsigned int nMany;//0 means return all

public:
    ExecOnVec(std::shared_ptr<CppModule> pModule, bool sort, unsigned int nMany)
            :
        CppModule(),
        pModule(pModule),
        sort(sort),
        nMany(nMany)
    {}//constructor

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);
    //overload
	std::vector<std::shared_ptr<Container>> getInputType() const;
    //overload
    std::shared_ptr<Container> getOutputType() const;
    
    std::string getName() const
    {
        return "ExecOnVec( " + pModule->getName() + ")";
    }
};//class

/**
 * @brief TODO:
 */
class Tail: public CppModule
{
private:
    std::shared_ptr<Container> type;
public:
    Tail(std::shared_ptr<Container> type)
            :
        CppModule(),
        type(type)
    {}//constructor

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);
    //overload
	std::vector<std::shared_ptr<Container>> getInputType() const;
    //overload
    std::shared_ptr<Container> getOutputType() const;
    
    std::string getName() const
    {
        return "Tail";
    }
};//class

/**
 * @brief Exposes the SweepAllReturnBest @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportExecOnVector();

#endif
