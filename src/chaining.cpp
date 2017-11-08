#include "chaining.h"

std::vector<ContainerType> Chaining::getInputType()
{
	return std::vector<ContainerType>{
			//the strip of consideration
			ContainerType::seeds,
		};
}//function

ContainerType Chaining::getOutputType()
{
	return ContainerType::seeds;
}//function


std::shared_ptr<Container> Chaining::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{
    std::shared_ptr<Seeds> pSeeds = std::static_pointer_cast<Seeds>(vpInput[0]);

    std::vector<RMQ<int64_t>::RMQData> data1;
    std::vector<RMQ<int64_t>::RMQData> data2;
    data1.push_back(RMQ<int64_t>::RMQData(-100000,0,nullptr, -1000000));
    data2.push_back(RMQ<int64_t>::RMQData(0,-100000,nullptr, -1000000));

    std::vector<std::shared_ptr<Chain>> chains;
    for(Seed seed : *pSeeds)
    {
        std::shared_ptr<Chain> chain = std::shared_ptr<Chain>(new Chain(seed));
        chains.push_back(chain);
        data1.push_back(t1(seed, chain));
        data2.push_back(t2(seed, chain));
    }//for
    for(RMQ<int64_t>::RMQData& d : data1)
        if(d.pChain != nullptr)
            d.pChain->t1 = &d;
    for(RMQ<int64_t>::RMQData& d : data2)
        if(d.pChain != nullptr)
            d.pChain->t2 = &d;
    std::shared_ptr<Chain> bestChain = chains[0];

    //first octant
    RMQ<int64_t> d1 = RMQ<int64_t>(data1);
    //2nd octant
    RMQ<int64_t> d2 = RMQ<int64_t>(data2);


    std::sort(
        chains.begin(), chains.end(),
        []
        (const std::shared_ptr<Chain> a, const std::shared_ptr<Chain> b)
        {
            if(a->s.end_ref() == b->s.end_ref())
                return a->s.end() < b->s.end();
            return a->s.end_ref() < b->s.end_ref();
        }//lambda
    );//function call

    //the actual chaining
    for(std::shared_ptr<Chain> chain : chains)
    {
        DEBUG_2(
            std::cout << "computing chain for (ref/query/size): " << chain->s.end_ref() << " "
                    << chain->s.end() << " "
                    << chain->s.size() << std::endl;
        )
        RMQ<int64_t>::RMQData& a = d1.rmq(
                -1000000,0,
                (int64_t)chain->s.end_ref()-(int64_t)chain->s.end(), chain->s.end()//TODO: replace with function
            );
        RMQ<int64_t>::RMQData& b = d2.rmq(
                0,(int64_t)chain->s.end()-(int64_t)chain->s.end_ref(),
                chain->s.end_ref(),-1000000 //TODO: replace with function
            );
        //using the RMQdata's to check for smaller insted of the chains
        //since we do not have to check for nullptrs this way
        if(a < b)
            chain->pred = b.pChain;
        else
            chain->pred = a.pChain;

        if(chain->pred == nullptr)
            continue;
        assert(chain->pred != chain);
        DEBUG_2(
            std::cout << "candidate: " << chain->pred->s.end_ref() << " "
                    << chain->pred->s.end() << " "
                    << chain->pred->s.size() << std::endl;
        )
        assert(chain->pred->s.end_ref() < chain->s.end_ref());
        assert(chain->pred->s.end() < chain->s.end());
        //update score
        int64_t x = chain->s.end_ref() - chain->pred->s.end_ref();
        int64_t y = chain->s.end() - chain->pred->s.end();

        int64_t possibleMatches = std::min(x,y);
        int64_t certainMatches = std::min(possibleMatches, (int64_t)chain->s.size());
        possibleMatches -= certainMatches;
        int64_t insOrDls = std::max(x,y) - std::min(x,y);

        //set the score for the new chain
        int64_t addScore = certainMatches * SCORE_MATCH
                - possibleMatches * COST_POSS_MATCH
                - insOrDls * COST_INS_DEL;

        DEBUG_2(
            std::cout << "best candidate would require " << insOrDls << " insertions or deletions "
            << possibleMatches << " possible matches and " << certainMatches << " certain matches for a score of " << chain->pred->score << "; add score is: " << addScore << std::endl;
        )

        //if chaining this seed results in a worse chain then don't...
        //this makes our global chaining into a local chaining
        if(chain->score > addScore + chain->pred->score)
        {
            DEBUG(
                std::cout << "best chain not worth"<< std::endl;
            )
            chain->pred = nullptr;
        }//if
        else{
            chain->score = addScore + chain->pred->score;
            chain->t1->score = chain->score - gc1(chain->s);
            chain->t2->score = chain->score - gc2(chain->s);
        }//else


        DEBUG(
            std::cout << "current score: "<< chain->score << std::endl;
        )


        //remember best chain
        if(*bestChain < *chain)
            bestChain = chain;
    }//for

    std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
    for(std::shared_ptr<Chain> chain : chains)
    {
        if(chain->pred != nullptr)
        {
            pRet->push_front(chain->s);
            pRet->push_front(chain->pred->s);
            pRet->push_front(Seed(0,0,0));
        }
    }//for
    while(bestChain != nullptr)
    {
        DEBUG_2(
            std::cout << "solution:" << std::endl;
            std::cout << bestChain->s.start() << " "
                    << bestChain->s.start_ref() << " "
                    << bestChain->s.size() << std::endl;
        )
        pRet->push_front(bestChain->s);
        bestChain = bestChain->pred;
    }//while
    DEBUG_2(
        std::cout << "done" << std::endl;
    )
    return pRet;
}//function

void exportChaining()
{
    //export the LineSweepContainer class
	boost::python::class_<Chaining, boost::python::bases<CppModule>>(
        "Chaining",
        "Uses chaining to remove contradicting "
        "matches within one strip of consideration.\n"
        "\n"
        "Execution:\n"
        "	Expects query, ref, strip_vec as input.\n"
        "		query: the query as NucleotideSequence\n"
        "		ref: the reference sequence as Pack\n"
        "		strip_vec: the areas that shall be evaluated as StripOfConsiderationVector\n"
        "	returns strip_vec.\n"
        "		strip_vec: the evaluated areas\n"
    );
}//function
