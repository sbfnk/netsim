/******************************************************************/
// RandomGenerator.hh
// implements random generators to be used by the simulation
/******************************************************************/

#ifndef RANDOMGENERATOR_HH
#define RANDOMGENERATOR_HH

#include <boost/random.hpp>
#include <time.h>

class RandomGenerator {

   public:

      RandomGenerator() {;}
      virtual ~RandomGenerator() {;}

      virtual double operator()() = 0;

};

class mt19937Uniform
   : virtual public RandomGenerator {

      typedef boost::mt19937 rng_t;
      typedef boost::uniform_real<double> dist_t;
      typedef boost::variate_generator<rng_t&, dist_t> generator_t;
      
      rng_t rng;
      dist_t dist;
      generator_t generator;
      
   public:

      mt19937Uniform()
         : rng(time(0)),
           dist(),
           generator(rng, dist)
      {}

      double operator()()
      {
         return generator();
      }
};

#endif
