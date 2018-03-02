#include <stdint.h>
#include <string>
#include <iostream>
#include <vector>

static inline uint32_t rotl32 ( uint32_t x, int8_t r )
{
      return (x << r) | (x >> (32 - r));
}

#define ROTL32(x,y) rotl32(x,y)

namespace multshift32 {
    
    inline uint32_t hash(uint32_t x, uint32_t l, uint32_t a)
    {
        return (a * x) >> (32 - l);
    }
}


namespace multshift64 {

    inline uint64_t hash(uint32_t x, uint64_t l, uint64_t a)
    {
        return (a * x) >> (64 - l);
    }
}

namespace multshift2wise {

    inline uint32_t hash(uint32_t x, uint32_t l, uint64_t a, uint64_t b)
    {
        return (a * x + b) >> (64 - l);
    }
}

class HashFunction {
    public:
        virtual uint32_t h1(uint32_t x) = 0;
        virtual uint32_t h2(uint32_t x) = 0;
        virtual std::string getDescription() = 0;
};



class Pol3: public HashFunction {
    public:
        
        Pol3()
        {
            p = (1LL<<48) - 1; 
            boost::uniform_int<uint64_t> dis(0, p);
            boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint64_t> > rand (g_gen, dis);

            a1 = rand();
            a2 = rand();
            b1 = rand();
            b2 = rand();
            c1 = rand();
            c2 = rand();
            
        }

        uint32_t h1(uint32_t x)
        {
            return ((((a1 * x) % p) * x + b1 * x + c1) % p);
        }

        uint32_t h2(uint32_t x)
        {
            return ((((a2 * x) % p) * x + b2 * x + c2) % p);
        }

        std::string getDescription()
        {
            return "Pol3";
        }

    private:
        uint64_t a1, a2, b1, b2, c1, c2, p;

};

class PolK: public HashFunction {
    public:
        PolK(uint32_t _k)
        {
            k = _k;
            p = (1LL<<61) - 1; 
            
            boost::uniform_int<uint64_t> dis(0, p);
            boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint64_t> > rand (g_gen, dis);
            
            a1 = new uint64_t[k];
            a2 = new uint64_t[k];

            for (uint32_t i = 0; i < k; i++)
            {
                a1[i] = rand();
                a2[i] = rand();
            }
        }

        virtual ~PolK()
        {
            delete[] a1;
            delete[] a2;
        }

        uint64_t cwtrick(uint32_t x, uint64_t a, uint64_t b)
        {
            uint64_t a0, a1, c0, c1;

            //multiply x with the lower 32bit of a
            a0 = (a & 0xFFFFFFFF) * x;
            // multiply x with the upper 32bit of a
            a1 = (a >> 32) * x;
            //the first 64 bit of the result
            c0 = a0 + (a1 << 32);
            //bits 33..72 of the result
            c1 = (a0 >> 32) + a1;
            // the modulo operation for this mersenne prime and addition of b
            return ((c0 & p) + (c1 >> 29) + b);
        }

        uint32_t h1(uint32_t x)
        {
            uint64_t res = a1[0];
            for (uint32_t i = 1; i < k; i++)
            {
                res = cwtrick(x, res, a1[i]);
            }
            res = (res & p) + (res >> 61);
            if (res >= p)
            {
                res -= p;
            }
            return (uint32_t) res;
        }

        uint32_t h2(uint32_t x)
        {
            uint64_t res = a2[0];
            for (uint32_t i = 1; i < k; i++)
            {
                res = cwtrick(x, res, a2[i]);
            }
            res = (res & p) + (res >> 61);
            if (res >= p)
            {
                res -= p;
            }
            return (uint32_t) res;
        }

        std::string getDescription()
        {
            return "k-ind-cw";
        }

    private:
        uint32_t k;
        uint64_t p;
        uint64_t* a1;
        uint64_t* a2;
};

class ADWunfixed: public HashFunction 
{
    public:

        ADWunfixed(unsigned short _k, unsigned short _c, uint32_t _l)
        {
            c = _c;
            l = _l;
            k = _k;
            size = 1 << l;

            z = new uint32_t[2 * c * size];
            
            boost::uniform_int<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
            boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint32_t> > rand (g_gen, dis);

            f = new PolK(k);
            g.reserve(k);

            for (uint32_t i = 0; i < c; i++)
            {
                g.push_back(new PolK(k));
            }

            for (uint32_t i = 0; i < (2 * c * size); i++)
            {
                // fill table with random values
                z[i] = rand();    
            }

        }

        virtual ~ADWunfixed()
        {
            for (uint32_t i = 0; i < c; i++)
            {
                delete g[i];
            }
            delete f;
            delete[] z;
        }

        uint32_t h1(uint32_t x) 
        {
            uint32_t res = f->h1(x); 
            for (uint32_t i = 0; i < c; i++)
            {
                res += z[i * size + (g[i]->h1(x) >> (32 - l))];
            }

            return (uint32_t) res;
        }
        
        uint32_t h2(uint32_t x) 
        {
            uint32_t res = f->h2(x); 
            for (uint32_t i = 0; i < c; i++)
            {
                res += z[(c + i) * size + (g[i]->h1(x) >> (32 - l))];
            }

            return (uint32_t) res;
        }

        
        std::string getDescription()
        {
            std::ostringstream convert; 
            convert << "ADW-unfixed-" << k << "-" << c << "-" << l;
            return convert.str();
        }

    private:
        unsigned short c,k;
        uint32_t l;
        uint32_t size;
        std::vector<PolK*> g;
        PolK* f;

        uint32_t* z;
};

boost::uniform_int<uint32_t> g_dis(0, std::numeric_limits<uint32_t>::max());

class FullyRandom: public HashFunction
{
    public:
        FullyRandom()
        {
            boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint32_t> > rand(g_gen, g_dis);
        }
	
        uint32_t h1(uint32_t __attribute__((__unused__)) x) 
        {
            return rand();
        }

        uint32_t h2(uint32_t __attribute__((__unused__)) x)
        {
            return rand();
        }
        
        std::string getDescription()
        {
            return "fully-random";
        }

    private:

};

class ADW: public HashFunction 
{

    public:

        ADW(unsigned short _c, uint32_t _l)
        {
            c = _c;
            l = _l;
            size = 1 << l;

            g = new uint32_t[c];
            z = new uint32_t[2 * c * size];
            
            boost::uniform_int<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
            boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint32_t> > rand (g_gen, dis);

            for (uint32_t i = 0; i < c; i++)
            {
                // choose odd numbers for c
                uint32_t a = 0;
                while (a % 2 == 0)
                {
                    a = rand();
                }
                g[i] = a;
            }

            for (uint32_t i = 0; i < (2 * c * size); i++)
            {
                // fill table with random values
                z[i] = rand();    
            }

            f1_a = rand();
            f1_b = rand();
            f2_a = rand();
            f2_b = rand();
        }

        virtual ~ADW()
        {
            delete[] g;
            delete[] z;
        }

        uint32_t h1(uint32_t x) 
        {
            uint32_t res = multshift2wise::hash(x, 32, f1_a, f1_b);
            for (uint32_t i = 0; i < c; i++)
            {
                res += z[i * size + multshift32::hash(x, l, g[i])];
            }

            return (uint32_t) res;
        }

        uint32_t h2(uint32_t x)
        {
            uint32_t res = multshift2wise::hash(x, 32, f2_a, f2_b);
            for (uint32_t i = 0; i < c; i++)
            {
                res += z[ (c + i) * size + multshift32::hash(x, l, g[i])];
            }

            return res;
        }
        
        std::string getDescription()
        {
            std::ostringstream convert;
            convert << "ADW-" << c << "-" << l;
            return convert.str();
        }

    private:
        unsigned short c;
        uint32_t l;
        uint32_t size;
        uint32_t* g;

        uint32_t* z;
        uint32_t f1_a;
        uint32_t f1_b;
        uint32_t f2_a;
        uint32_t f2_b;


};

class SimpleTab8: public HashFunction
{
    public:
    SimpleTab8()
    {
        z1 = new uint32_t[1<<10];
        z2 = new uint32_t[1<<10];
        
        boost::uniform_int<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
        boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint32_t> > rand (g_gen, dis);

        for (uint32_t i = 0; i < 1<<10; i++)
        {
            z1[i] = rand();
        }
        for (uint32_t i = 0; i < 1<<10; i++)
        {
            z2[i] = rand();
        }
    }

    virtual ~SimpleTab8()
    {
        delete[] z1;
        delete[] z2;
    }

    uint32_t h1(uint32_t x)
    {
//        uint32_t res = 0;
//        for (int i = 0; i < c; i++)
//        {
//            res ^= z1[i * 256 + (x & 0xFF)];
//            x = x >> 8;
//        }
//        return res;
          return z1[x & 0xFF] ^ z1[256 + ((x >> 8) & 0xFF)] ^ z1[512 + ((x >> 16) & 0xFF)]
                    ^ z1[768 + (x >> 24)];
    }
    
    uint32_t h2(uint32_t x)
    {
//        uint32_t res = 0;
//        for (int i = 0; i < c; i++)
//        {
//            res ^= z2[i * 256 + (x & 0xFF)];
//            x = x >> 8;
//        }
//        return res;
        return z2[x & 0xFF] ^ z2[256 + ((x >> 8) & 0xFF)] ^ z2[512 + ((x >> 16) & 0xFF)]
                    ^ z2[768 + (x >> 24)];
    }
        
    std::string getDescription()
        {
            return "simp-tab-8";
        }

    private:
        uint32_t* z1;
        uint32_t* z2;

};

class SimpleTab16: public HashFunction
{
    public:
    SimpleTab16()
    {
        z1 = new uint32_t[1<<17];
        z2 = new uint32_t[1<<17];
        
        boost::uniform_int<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
        boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint32_t> > rand (g_gen, dis);

        for (uint32_t i = 0; i < 1<<17; i++)
        {
            z1[i] = rand();
        }
        for (uint32_t i = 0; i < 1<<17; i++)
        {
            z2[i] = rand();
        }
    }

    virtual ~SimpleTab16()
    {
        delete[] z1;
        delete[] z2;
    }

    uint32_t h1(uint32_t x)
    {
//        uint32_t res = 0;
//        for (int i = 0; i < c; i++)
//        {
//            res ^= z1[i * 256 + (x & 0xFF)];
//            x = x >> 8;
//        }
//        return res;
          return z1[x & 0xFFFF] ^ z1[(1 << 16) + (x >> 16)];
    }
    
    uint32_t h2(uint32_t x)
    {
//        uint32_t res = 0;
//        for (int i = 0; i < c; i++)
//        {
//            res ^= z2[i * 256 + (x & 0xFF)];
//            x = x >> 8;
//        }
//        return res;
          return z2[x & 0xFFFF] ^ z2[(1 << 16) + (x >> 16)];
    }
        
    std::string getDescription()
        {
            return "simp-tab-16";
        }

    private:
        uint32_t* z1;
        uint32_t* z2;

};


class Murmur3: public HashFunction
{
    private:
        uint32_t h1_seed, h2_seed;

        // adapted from
        // https://raw.githubusercontent.com/PeterScott/murmur3/master/murmur3.c

        uint32_t MurmurHash3_x86_32 ( const void * key, int len,
                                uint32_t seed)
        {
        const uint8_t * data = (const uint8_t*)key;
        const int nblocks = len / 4;
        int i;

        uint32_t h1 = seed;

        uint32_t c1 = 0xcc9e2d51;
        uint32_t c2 = 0x1b873593;

        //----------
        // body

        const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

        for(i = -nblocks; i; i++)
        {
            uint32_t k1 = blocks[i];

            k1 *= c1;
            k1 = ROTL32(k1,15);
            k1 *= c2;
            
            h1 ^= k1;
            h1 = ROTL32(h1,13); 
            h1 = h1*5+0xe6546b64;
        }

        //----------
        // tail

        const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

        uint32_t k1 = 0;

        switch(len & 3)
        {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
        };

        //----------
        // finalization

        h1 ^= len;

        h1 ^= h1 >> 16;
        h1 *= 0x85ebca6b;
        h1 ^= h1 >> 13;
        h1 *= 0xc2b2ae35;
        h1 ^= h1 >> 16;

        return h1;
        } 
    public:
        Murmur3()
        {
            boost::mt19937_64 gen;
            gen.seed(time(NULL));
            boost::uniform_int<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
            boost::variate_generator<boost::mt19937_64&,boost::uniform_int<uint32_t> > rand (gen, dis);

            h1_seed = rand();
            h2_seed = rand();
        }

	virtual ~Murmur3()
	{
	}


        uint32_t h1(uint32_t x)
        {
            return MurmurHash3_x86_32 (&x, 4, h1_seed);
        }

        uint32_t h2(uint32_t x)
        {
            return MurmurHash3_x86_32 (&x, 4, h2_seed);
        }

        std::string getDescription()
        {
            return "Murmur3";
        }

};
