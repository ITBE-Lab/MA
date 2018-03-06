/** 
 * @file alignment.h
 * @brief Implements a Container that holds a finished alignment.
 * @author Markus Schmidt
 */

#ifndef MINIMIZER_HASH_H
#define MINIMIZER_HASH_H

#include "container/segment.h"
#include "container/seed.h"
#include "util/support.h"
#include <map>
#include <fstream>
#include <queue>
#include <vector>

#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

namespace libMA
{

    /*std::function<int(uint8_t a, uint8_t b, unsigned int iPos)>*/ 
    template<size_t k>
    class Minimizer
    {
    public:
        /**
         * @brief The minimizer sequence
         */
        uint8_t xSeq[k];
        bool bRevComp = false;

        static uint32_t maxIndex()
        {
            return std::numeric_limits<uint32_t>::max();
        }//function

        bool isAmbiguous()
        {
            for(unsigned int i = 0; i < k; i++)
                if(xSeq[i] >= 4)//found ambiguous character..
                    return true;
            return false;
        }//function

        uint32_t toIndex(const uint8_t* pSeq) const
        {
            //check for potential overflows
            assert(k*2 < 32);
            /*
            * ordering in nucSeq is 'A', 'C', 'G', 'T'
            * we want to have C, A, T, G however therefore translation table:
            */
            static const uint8_t translate[4] = {1, 0, 3, 2};

            uint32_t uiRet = 0;
            /*
            * In DNA sequences, the letters C and G often occur less frequently than
            * A and T.We assign the values 0, 1, 2, 3 to C, A, T, G, respectively,
            * for the odd numbered bases of k-mers, and reverse the
            * ordering for even numbered bases. This tends to start minimizers
            * with the valuable (in the sense of the significance
            * of a match) letters C and G, and makes the minimum k-mer
            * CGCGCG.... There are many other possibilities.
            */
            for(unsigned int i = 0; i < k; i++)
            {
                if(pSeq[i] >= 4)
                    return maxIndex();
                uint8_t uiNuc = translate[pSeq[i]];
                uiRet <<= 2;
                if(i % 2 == 0)
                    uiRet |= uiNuc;
                else
                    uiRet |= 3 - uiNuc;
            }//for
            return uiRet;
        }//function

        uint32_t toIndex() const
        {
            return toIndex(xSeq);
        }//function

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator<(const Minimizer& rOther) const
        {
            return toIndex() < rOther.toIndex();
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator<=(const Minimizer& rOther) const
        {
            return toIndex() <= rOther.toIndex();
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator>(const Minimizer& rOther) const
        {
            return toIndex() > rOther.toIndex();
        }//operator

        void operator=(const uint8_t* pOtherSeq)
        {
            for(unsigned int i = 0; i < k; i++)
                xSeq[i] = pOtherSeq[i];
        }//operator

        void operator=(const Minimizer& rOther)
        {
            operator=(rOther.xSeq);
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator==(const Minimizer& rOther) const
        {
            return toIndex() == rOther.toIndex();
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator!=(const Minimizer& rOther) const
        {
            return toIndex() != rOther.toIndex();
        }//operator

        Minimizer(std::shared_ptr<NucSeq> pInit, nucSeqIndex uiStart)
        {
            assert(pInit->length() >= k + uiStart);
            //this makes sure that a sequence and it's reverse complement have the same minimizer
            uint8_t xSeq2[k];
            for(int i = 0; i < (int)k; i++)
            {
                xSeq[i] = (*pInit)[i + uiStart];
                assert( ((int)k)-(i+1) >= 0);
                assert( ((int)k)-(i+1) < (int)k);
                xSeq2[((int)k)-(i+1)] = (uint8_t)complement((*pInit)[i + uiStart]);
            }//for

            //if the reverse complement minimizer is smaller switch to that
            if(toIndex() < toIndex(xSeq2))
            {
                *this = xSeq2;
                bRevComp = true;
            }//if
        }//constructor

        Minimizer(std::shared_ptr<NucSeq> pInit)
                :
            Minimizer(pInit, 0)
        {}//constructor

        Minimizer()
        {}//default constructor
    };//class

    template<size_t w, size_t k>
    class MinimizersHash;

    template<size_t w, size_t k>
    class MinimizersVector 
            :
        public std::vector<std::pair<Minimizer<k>, nucSeqIndex>>,
        public Container
    {
    public:

        MinimizersVector()
        {}//constructor

        /**
         * creates a hash table from the vector
         */
        std::shared_ptr<MinimizersHash<w,k>> toHash(
            nucSeqIndex uiRefSize 
            DEBUG_PARAM(std::shared_ptr<Pack> pPack)
        );

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<MinimizersVector>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "MinimizersVector";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new MinimizersVector());
        }//function
    };//class

    template<size_t w, size_t k>
    class MinimizersHash : public Container
    {
    public:
        static const size_t FILE_VERSION = 1;
        const nucSeqIndex uiRefSize;

        //it might be better to actually use a hash table here...
        std::vector<nucSeqIndex> vValues;
        //pair: key , index in value of the key
        std::vector<std::pair<uint32_t,uint32_t>> vKeys;


        MinimizersHash(nucSeqIndex uiRefSize = 0, uint32_t vValuesSize = 0, uint32_t vKeysSize = 0)
                :
            uiRefSize(uiRefSize),
            vValues(),
            vKeys()
        {
            vValues.resize(vValuesSize);
            vKeys.resize(vKeysSize);
        }//constructor

        //this should never be copied
        MinimizersHash(MinimizersHash& rOther) = delete;


        unsigned int keyLen() const
        {
            return vKeys.size();
        }//function

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<MinimizersHash>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "MinimizersHash";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new MinimizersHash());
        }//function
    public:
        static std::shared_ptr<MinimizersHash> fromFile(std::string sFileName)
        {
            // create and open an archive for input
            std::ifstream ifs(sFileName.c_str(), std::ios::binary);

            size_t check;
            ifs.read((char*)&check, sizeof(size_t));
            if(check != FILE_VERSION)
            {
                std::cerr << "cant read file 1: " << check << std::endl;
                throw "@todo proper exception";
            }
            nucSeqIndex uiRefSize;
            ifs.read((char*)&uiRefSize, sizeof(nucSeqIndex));
            ifs.read((char*)&check, sizeof(size_t));
            if(k != check)
            {
                std::cerr << "cant read file 2" << std::endl;
                throw "@todo proper exception";
            }
            ifs.read((char*)&check, sizeof(size_t));
            if(w != check)
            {
                std::cerr << "cant read file 3" << std::endl;
                throw "@todo proper exception";
            }

            //actual data
            size_t uiKeySize;
            size_t uiValueSize;
            ifs.read((char*)&uiKeySize, sizeof(size_t));
            ifs.read((char*)&uiValueSize, sizeof(size_t));
            auto pRet = std::make_shared<MinimizersHash>(uiRefSize, uiValueSize, uiKeySize);

            ifs.read((char*)&pRet->vKeys[0], pRet->vKeys.size() * sizeof(uint32_t)*2);
            ifs.read((char*)&pRet->vValues[0], pRet->vValues.size() * sizeof(nucSeqIndex));

            ifs.close();

            assert(!pRet->vKeys.empty());
            assert(!pRet->vValues.empty());
            return pRet;
        }//function

        void toFile(std::string sFileName)
        {
            std::ofstream ofs(sFileName.c_str(), std::ios::trunc | std::ios::binary);
            // write class instance to archive
            size_t out = FILE_VERSION;
            ofs.write((char*)&out, sizeof(size_t));
            ofs.write((char*)&uiRefSize, sizeof(nucSeqIndex));
            out = k;
            ofs.write((char*)&out, sizeof(size_t));
            out = w;
            ofs.write((char*)&out, sizeof(size_t));
            out = vKeys.size();
            ofs.write((char*)&out, sizeof(size_t));
            out = vValues.size();
            ofs.write((char*)&out, sizeof(size_t));
            //write the keys; they are uint32_t pairs in a vector therefore we can mass write
            ofs.write((char*)&vKeys[0], vKeys.size() * sizeof(uint32_t)*2);
            //write the values; they are nucSeqIndex in a vector so we can mass write
            ofs.write((char*)&vValues[0], vValues.size() * sizeof(nucSeqIndex));

            ofs.close();
        }//function

        //returns a pointer to the position and the size
        std::pair<nucSeqIndex*, uint32_t> operator[](Minimizer<k>& rKey)
        {
            assert(!vKeys.empty());
            auto xItLower = std::lower_bound(
                vKeys.begin(),
                vKeys.end(),
                rKey.toIndex(),
                []
                (std::pair<uint32_t,uint32_t> xKey , uint32_t xVal)
                {
                    return xKey.first < xVal;
                }//lambda
            );//lower bound function call
            if(xItLower == vKeys.end() || xItLower->first != rKey.toIndex())
                return std::make_pair(nullptr, 0);
            auto xItNext = xItLower + 1;
            assert(xItNext->first != xItLower->first);
            uint32_t uiSize = xItNext->second - xItLower->second;
            assert(uiSize > 0);
            return std::make_pair(&vValues[xItLower->second], uiSize);
        }//function

        class iterator : public Container
        {
        public:
            const nucSeqIndex uiRefSize;
            typedef std::tuple< //content
                    nucSeqIndex*, // reference position
                    nucSeqIndex, // query position
                    uint32_t, // remaining reference positions
                    bool //reverse complement ref positions
                > queueContent;
            typedef std::function<bool(queueContent,queueContent)> queueCompFunc;
            std::priority_queue<queueContent, std::vector<queueContent>, queueCompFunc> xQueue;
            iterator(
                    std::shared_ptr<MinimizersVector<w,k>> pQueryMinimizers, 
                    MinimizersHash& rHash
                )
                    :
                uiRefSize(rHash.uiRefSize),
                xQueue(
                    [&]
                    (
                        queueContent a,
                        queueContent b
                    )
                    {
                        auto aRefPos = 
                            std::get<3>(a) ? uiRefSize - *std::get<0>(a) : *std::get<0>(a);
                        auto bRefPos = 
                            std::get<3>(b) ? uiRefSize - *std::get<0>(b) : *std::get<0>(b);
                        //@todo remove me (temporary sorting for SOCs...)
                        return aRefPos + std::get<1>(b) > bRefPos + std::get<1>(a);

                        if(aRefPos == bRefPos)
                            //sort so that larger query positions come first
                            return std::get<1>(a) > std::get<1>(b);
                        //sort so that smaller reference positions come first (with priority)
                        return aRefPos < bRefPos;
                    }//lambda
                )//constructor for priority queue
            {
                auto xIter = pQueryMinimizers->begin();
                while(xIter != pQueryMinimizers->end())
                {
                    auto xHashPair = rHash[xIter->first];
                    if(xHashPair.second > 0)
                        xQueue.push(std::make_tuple(
                            xIter->first.bRevComp ? 
                                 // we want to reverse the order for complemented minimizers
                                xHashPair.first + xHashPair.second - 1
                                : xHashPair.first, //reference positions pointer
                            xIter->second, // query position
                            xHashPair.second,// amount reference positions
                            xIter->first.bRevComp// was the minimizer reverse complemented
                        ));
                    xIter++;
                }//while
            }//constructor

            iterator()
                    :
                uiRefSize(0),
                xQueue()
            {}//default constructor
            
            //overload
            bool canCast(std::shared_ptr<Container> c) const
            {
                return std::dynamic_pointer_cast<iterator>(c) != nullptr;
            }//function

            //overload
            std::string getTypeName() const
            {
                return "MinimizersHash::iterator";
            }//function

            //overload
            std::shared_ptr<Container> getType() const
            {
                return std::shared_ptr<Container>(new iterator());
            }//function

            bool empty()
            {
                return xQueue.empty();
            }//function

            Seed operator*() const
            {
                assert(!xQueue.empty());
                return Seed(
                    std::get<1>(xQueue.top()), 
                    k, 
                    //check weather we need to reverse complement the reference position
                    std::get<3>(xQueue.top()) ? 
                        uiRefSize - (*std::get<0>(xQueue.top()) + k)
                        :
                        *std::get<0>(xQueue.top())
                );
            }//operator

            iterator& operator++()
            {
                //copy the top element
                queueContent xTuple;
                xTuple = xQueue.top();
                //remove the top element from the priority queue
                xQueue.pop();
                //increment the vector iterator of the top element
                if(std::get<3>(xQueue.top()))
                    //remember: we need to reverse the order for complemented minimizers
                    std::get<0>(xTuple)--;
                else
                    std::get<0>(xTuple)++;
                //if there are still elements left in the reference vector
                std::get<2>(xTuple) = std::get<2>(xTuple)-1;
                if(std::get<2>(xTuple) > 0)
                    //readd the previous top element
                    xQueue.push(xTuple);
                return *this;
            }//operator
        };//class

        iterator begin(std::shared_ptr<MinimizersVector<w,k>> pQueryMinimizers)
        {
            return iterator(pQueryMinimizers, *this);
        }//function

        std::shared_ptr<Seeds> toSeeds(std::shared_ptr<MinimizersVector<w,k>> pQueryMinimizers)
        {
            auto it = begin(pQueryMinimizers);
            auto pRet = std::make_shared<Seeds>();
            while(!it.empty())
            {
                pRet->push_back(*it);
                ++it;
            }//while
            return pRet;
        }//function

    };//inner class

    template<size_t w, size_t k>
    std::shared_ptr<MinimizersHash<w,k>> MinimizersVector<w,k>::toHash(
        nucSeqIndex uiRefSize
        DEBUG_PARAM(std::shared_ptr<Pack> pPack)
        )
    {
        assert(!this->empty());
        DEBUG(
            for(auto xElement : *this)
                xElement.first.toIndex();
        )//DEBUG
        //sort now so that there is no need to sort for the SOC
        //sideeffect: lets us fill the hashtable much faster and easier
        std::cout << "sorting..." << std::endl;
        std::sort(
            this->begin(), 
            this->end(),
            [&]
            (std::pair<Minimizer<k>, nucSeqIndex>& a, 
                std::pair<Minimizer<k>, nucSeqIndex>& b)
            {
                if(a.first.toIndex() == b.first.toIndex())
                    //the soc order
                    return (a.first.bRevComp ? uiRefSize - (a.second + k) : a.second) < 
                        (b.first.bRevComp ? uiRefSize - (b.second + k) : b.second);
                //the fill order (this takes priority)
                return a.first.toIndex() < b.first.toIndex();
            }//lambda
        );//sort function call
        std::cout << "done sorting" << std::endl;
        auto pRet = std::make_shared<MinimizersHash<w,k>>(uiRefSize, this->size());
        //remember the vector in the hash table that we are currently filling
        auto pCurrent = this->begin();
        //this is in order to check for duplicates
        // set to max value so that the first element does not trigger the check...
        nucSeqIndex uiLast = (nucSeqIndex)-1;
        //for the progress print
        unsigned int i = 0;
        //fill in the first key
        pRet->vKeys.emplace_back(pCurrent->first.toIndex(), i);
        pRet->vValues[i] = pCurrent->first.bRevComp ? 
            uiRefSize - (pCurrent->second + k) : pCurrent->second;
        DEBUG(
            auto xCheck = Minimizer<k>(pPack->vExtract(pRet->vValues[i], pRet->vValues[i]+k), 0);
            if(pCurrent->first != xCheck )
            {
                std::cout << k << "-Minimizer complement:" << pCurrent->first.bRevComp << " wrong ";
                for(size_t f = 0; f < k; f++)
                    std::cout << (int)pCurrent->first.xSeq[f];
                std::cout << " != ";
                for(size_t f = 0; f < k; f++)
                    std::cout << (int)xCheck.xSeq[f];
                std::cout << std::endl;
                std::cout << pPack->vExtract(pRet->vValues[i]-k, pRet->vValues[i])->toString() 
                    << std::endl;
                std::cout << pPack->vExtract(pCurrent->second, pCurrent->second+k)->toString() 
                    << std::endl;
                exit(0);
            }//if
        )//DEBUG
        //fill in all other elements
        while(++pCurrent != this->end())
        {
            if(i % 1000000 == 0)
                std::cout << i/1000000 << "/" << this->size()/1000000 << std::endl;
            //check that there are no duplicates and ignore all minimizers with Ns 
            if(uiLast == pCurrent->second || pCurrent->first.isAmbiguous())
            {
                //we skipped one element so we need to remove that element from the vector
                pRet->vValues.pop_back();
                continue;
            }//if
            uiLast = pCurrent->second;
            //in this case we need to create a new vector in the hash table
            if(pCurrent->first.toIndex() != pRet->vKeys.back().first)
                pRet->vKeys.emplace_back(pCurrent->first.toIndex(), i);
            //save the current element
            pRet->vValues[++i] = pCurrent->first.bRevComp ? 
                uiRefSize - (pCurrent->second + k) : pCurrent->second;
            DEBUG(
                auto xCheck = Minimizer<k>(pPack->vExtract(pRet->vValues[i], pRet->vValues[i]+k), 0);
                if(pCurrent->first != xCheck )
                {
                    std::cout << k << "-Minimizer complement:" 
                        << pCurrent->first.bRevComp << " wrong ";
                    for(size_t f = 0; f < k; f++)
                        std::cout << (int)pCurrent->first.xSeq[f];
                    std::cout << " != ";
                    for(size_t f = 0; f < k; f++)
                        std::cout << (int)xCheck.xSeq[f];
                    std::cout << std::endl;
                    std::cout << pPack->vExtract(pRet->vValues[i]-k, pRet->vValues[i])->toString() 
                        << std::endl;
                    std::cout << pPack->vExtract(pCurrent->second, pCurrent->second+k)->toString() 
                        << std::endl;
                    exit(0);
                }//if
            )//DEBUG
        }//while
        pRet->vKeys.emplace_back(Minimizer<k>::maxIndex(), i);
        return pRet;
    }//function

}//namespace libMA

/**
 * @brief Exposes the MinimizersHash container to boost python.
 * @ingroup export
 */
void exportMinimizersHash();

#endif