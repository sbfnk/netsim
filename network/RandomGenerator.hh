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
      typedef boost::uniform_01<rng_t, double> generator_t;
      
      rng_t rng;
      generator_t generator;
      
   public:

      mt19937Uniform()
         : rng(time(0)),
           generator(rng)
      {}

      double operator()()
      {
         return generator();
      }
};

#endif
