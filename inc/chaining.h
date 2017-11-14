/**
 * @file chaining.h
 * @author Markus Schmidt
 * @details
 * check out the respective paper.
 */


#ifndef CHAINING_H
#define CHAINING_H

#include "rmq.h"
#include "sequence.h"
#include "pack.h"
#include "cppModule.h"



class Chain
{
public:
    Seed s;
    int score;
	std::shared_ptr<Chain> pred;
	RMQ<int64_t>::RMQData *t1, *t2;
	

    Chain(Seed s)
            :
        s(s),
        score(s.size() * SCORE_MATCH),
        pred()
    {}//constructor

    const bool operator< (const Chain& other) const
    {
        return score < other.score;
    }//function
};//class

class Chaining : public CppModule
{
private:

    inline int64_t gc1_end(const Seed s) const
    {
        int64_t x = (int64_t)s.end_ref()-(int64_t)s.end();
        int64_t y = s.end();
        return COST_INS_DEL * x + COST_POSS_MATCH * y;
    }//function

    inline int64_t gc2_end(const Seed s) const
    {
        int64_t x = s.end_ref();
        int64_t y = (int64_t)s.end()-(int64_t)s.end_ref();
        return COST_POSS_MATCH * x + COST_INS_DEL * y;
    }//function

    inline int64_t gc1_start(const Seed s) const
    {
        int64_t x = (int64_t)s.start_ref()-(int64_t)s.start();
        int64_t y = s.start();
        return COST_INS_DEL * x + COST_POSS_MATCH * y;
    }//function

    inline int64_t gc2_start(const Seed s) const
    {
        int64_t x = s.start_ref();
        int64_t y = (int64_t)s.start()-(int64_t)s.start_ref();
        return COST_POSS_MATCH * x + COST_INS_DEL * y;
    }//function

    inline RMQ<int64_t>::RMQData t1(Seed seed, std::shared_ptr<Chain> chain) const
    {
        return RMQ<int64_t>::RMQData(
                (int64_t)seed.end_ref()-(int64_t)seed.end(),
                seed.end(),
                chain,
                chain->score + gc1_end(seed)
            );
    }//function

    inline RMQ<int64_t>::RMQData t2(Seed seed, std::shared_ptr<Chain> chain) const
    {
        return RMQ<int64_t>::RMQData(
                seed.end_ref(),
                (int64_t)seed.end()-(int64_t)seed.end_ref(),
                chain,
                chain->score + gc2_end(seed)
            );
    }//function

public:
    //overload
	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput);

    //overload
    std::vector<std::shared_ptr<Container>> getInputType();

    //overload
    std::shared_ptr<Container> getOutputType();



};//class
/**
 * @brief Exposes the Chaining @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportChaining();

#endif