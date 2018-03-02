#include<vector>
#include<cmath>
#include <boost/random.hpp>

#define MAXLOOP 1000

static boost::mt19937_64 g_gen;

#include "hashfunctions.h"
#include "tools/timer.h"
#include "tools/papi.h"

//#define DEBUG 0


inline size_t rand_int(int bound)
{
    boost::uniform_int<size_t> dis(0, bound - 1);
    return dis(g_gen);
}

namespace cuckoohashing {

    uint32_t* t1;
    uint32_t* t2;

    uint32_t m;

    HashFunction* h;
    std::vector<uint32_t> stash;

    void init(uint32_t _m, HashFunction* _h) 
    {
        h = _h;
        m = _m;

        t1 = new uint32_t[m];
        t2 = new uint32_t[m];


        for (uint32_t i = 0; i < m; i++)
        {
            t1[i] = 0;
            t2[i] = 0;
        }
    }

    bool lookup(uint32_t key)
    {
        if (t1[h->h1(key) % m] == key)
            return true;
        if (t2[h->h2(key) % m] == key)
            return true;
        for (uint32_t i = 0; i < stash.size(); i++)
	{
            if (stash[i] == key)
	    {
                return true;
	    }
        }
	return false;
    }

    void remove(uint32_t key)
    {
        if (t1[h->h1(key) % m] == key)
            t1[h->h1(key) % m] = 0;
        if (t2[h->h2(key) % m] == key)
            t2[h->h2(key) % m] = 0;
        for (uint32_t i = 0; i < stash.size(); i++)
            if (stash[i] == key)
                stash.erase(stash.begin() + i); 
    }

    void insert(uint64_t key)
    {
        uint64_t tmp = 0;
        uint64_t hash = 0;
        uint8_t i = 1;
        uint16_t c = 0;
#ifdef DEBUG            
        std::cout <<
            "Key: " << key << 
            " h1: " << h->h1(key) % m << 
            " h2: " << h->h2(key) % m <<
            std::endl;
#endif

        while (c < MAXLOOP)
        {
                
            if (i == 1)
            {
                hash = h->h1(key) % m;
                tmp = t1[hash];
                t1[hash] = key;
            }
            else
            {
                hash = h->h2(key) % m;
                tmp = t2[hash];
                t2[hash] = key;
            }
            key = tmp;
            if (key == 0)
                break;
            c++;
            i = 3 - i;
        }
        if (key != 0)
        {
            stash.push_back(key);
        }
    }
}

std::vector<uint32_t> create_hypercube(int l)
{
    std::vector<uint32_t> keys;

    for (int i1 = 0; i1 < l; i1++)
    {
        for (int i2 = 0; i2 < l; i2++)
        {
            for (int i3 = 0; i3 < l; i3++)
            {
                for (int i4 = 0; i4 < l; i4++)
                {
                    keys.push_back((i1 << 24) | (i2 << 16) | (i3 << 8) | i4);                
                }
            }
        }
    }
    return keys;
}

std::vector<uint32_t> create_keys(int n)
{
    std::vector<uint32_t> keys;
    keys.reserve(n);
    for (int i = 1; i < n + 1; i++)
    {
        keys.push_back(i);
    }
    return keys;
}

int main(int argc, char** argv)
{
    if (argc < 3 || argc > 4)
    {
        std::cout << "Usage: seed method [n]" << std::endl;
	std::cout << "If [n] is not given, input will be the hypercube [32]^4" << std::endl;
	std::cout << "Available Methods: \n" 
		  << "\t 0 - simple tabulation 8-bit char \n" 
		  << "\t 1 - simple tabulation 16-bit char \n" 
		  << "\t 2 - Murmur3\n" 
		  << "\t 3 - degree 3 polynomial \n" 
		  << "\t 4 - degree 20 polynomial \n" 
		  << "\t 5 - Z, 3 tables with sqrt(n) entries\n" 
		  << "\t 6 - Z, 4 tables with n^{1/4} entries\n" 
		  << "\t 7 - Z, 8 tables with sqrt(n) entries\n" 
		  << "\t 8 - Z, 16 tables with n^{1/4} entries\n" 
		  << "\t 9 - Z, 1 table, 6-wise independence, with sqrt(n) entries\n" 
		  << "\t 10 - Z, 1 table, 12-wise independence, with n^{1/4} entries\n" 
		  << "\t 11 - Z, 1 table, 16-wise independence, with sqrt(n) entries\n" 
		  << "\t 12 - fully random, just returns random hash values. Warning: Does not store the key-value mapping!\n" << std::endl;
        return 0;
    }

    uint32_t seed = atoi(argv[1]);
    int method = atoi(argv[2]);
    std::vector<uint32_t> keys;
    HashFunction* h;
    uint32_t n, m;
    
    g_gen.seed(seed);


    if (argc == 4)
    {
        keys = create_keys(atoi(argv[3]));
    }
    else
    {
        keys = create_hypercube(32);
    }
    
    std::random_shuffle(keys.begin(), keys.end(), rand_int);
    n = keys.size();
    m = 1.005 * n;

    int l1 = (int) ceil(log2(std::sqrt(n)));
    int l2 = (int) ceil(log2(std::pow(n, 0.25)));

    switch (method)
    {
        case 0:
            h = new SimpleTab8();
            break;
        case 1:
            h = new SimpleTab16();
            break;
        case 2:
            h = new Murmur3();
            break;
        case 3:
            h =  new PolK(3);
            break;
        case 4:
            h = new PolK(20);
            break;
        case 5:
            // fail prob. 1/n^{1/2}
            h = new ADW(3, l1);
            break;
        case 6:
            //fail prob. 1/n^{1/3}
            h = new ADW(4, l2);
            break;
        case 7:
            //fail prob. 1/n^{3}
            h = new ADW(8, l1);
            break;
        case 8:
            // fail prob 1/n^3
            h = new ADW(16, l2);
            break;
        case 9:
            // fail prob 1/n^{1/3}
            h = new ADWunfixed(6, 1, l1);
            break;
        case 10:
            // fail prob 1/n^{1/3}
            h = new ADWunfixed(12, 1, l2);
            break;
        case 11:
            // fail prob 1/n^3
            h = new ADWunfixed(16, 1, l1);
            break;
        case 12:
            h = new FullyRandom();
            break;
        default:
            std::cerr << " Method not supported " << std::endl;
            return 0;
    }
    

    cuckoohashing::init(m, h);
    
    ClockIntervalBase<CLOCK_MONOTONIC> timer;
    ClockIntervalBase<CLOCK_PROCESS_CPUTIME_ID> cpu_timer;
    PApiWrapper papi;

    papi.add_event(PAPI_TOT_INS); // Total Instructions
    papi.add_event(PAPI_TOT_CYC); // Total Cycles

    papi.add_event(PAPI_L2_TCM); // L2 Total Cache Misses
    papi.add_event(PAPI_L1_TCM); // L2 Total Cache Misses

    //papi.add_event(PAPI_BR_CN); // Conditional branch instructions executed
    //papi.add_event(PAPI_BR_TKN); // Conditional branch instructions taken
    //papi.add_event(PAPI_BR_NTK); // Conditional branch instructions not taken
    //papi.add_event(PAPI_BR_MSP); // Conditional branch instructions mispred
    //papi.add_event(PAPI_BR_PRC); // Conditional branch instructions correctly predicted
    

    


    papi.start(), cpu_timer.start(), timer.start();
    for (std::vector<uint32_t>::iterator it = keys.begin() ; it != keys.end(); it++)
    {
        cuckoohashing::insert(*it);
    }
    timer.stop(), cpu_timer.stop(), papi.stop();

    std::cout <<
                " m=" << m <<
                " n=" << n <<
                " seed=" << seed <<
                " h=" << method << 
                " name=" << h->getDescription() << 
                " time=" << timer.delta() <<
                " cpu_time=" << cpu_timer.delta() <<
                " stash_size=" << cuckoohashing::stash.size();
    
    
    for (size_t i = 0 ; i < papi.get_num_counter(); ++i)
    {
        std::cout << " " << papi.get_counter_name(i) << "=" << papi.get_counter_result(i);
    }

    std::cout << std::endl;

//    for (std::vector<uint32_t>::iterator it = keys.begin() ; it != keys.end(); it++)
//    {
//        if (!cuckoohashing::lookup(*it))
//        {
//            std::cout << "Couldn't find " << *it;
//            break;
//        }
//    }
    switch (method)
    {
        case 0:
            delete (SimpleTab8*) h;
            break;
        case 1:
            delete (SimpleTab16*) h;
            break;
        case 2:
            delete (Murmur3*)h;
            break;
        case 3:
            delete (PolK*)h;
            break;
        case 4:
            delete (PolK*) h;
            break;
        case 5:
        case 6:
        case 7:
        case 8:
            delete (ADW*) h;
            break;
        case 9:
        case 10:
        case 11:
            delete (ADWunfixed*) h;
            break;
    }
    delete[] cuckoohashing::t1;
    delete[] cuckoohashing::t2;

    return 0;
}
