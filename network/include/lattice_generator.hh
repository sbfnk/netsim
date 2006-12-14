/******************************************************************/
// lattice_generator.hh
// implements a lattice iterator for the boost graph library
/******************************************************************/
#ifndef LATTICE_GENERATOR_HH
#define LATTICE_GENERATOR_HH

#include <boost/graph/graph_traits.hpp>
#include <math.h>
#include <vector>

#include <iostream>

namespace boost {

   /******************************************************************/
   // generateLattice
   // creates the lattice graph
   /******************************************************************/
   template <typename Graph>
   class lattice_iterator
   {
         typedef typename graph_traits<Graph>::vertices_size_type
         vertices_size_type;
      public:
         typedef std::pair<vertices_size_type, vertices_size_type> value_type;
         typedef const value_type& reference;
         typedef const value_type* pointer;

         lattice_iterator(): past_the_end(true) {}
         lattice_iterator(vertices_size_type sl, unsigned int dim = 2,
                          bool pb = true)
            : past_the_end(false),periodic_boundary(pb),
              sideLength(sl), dimensions(dim),
              dimCount(0), source(0),target(1),current(0,1)
         {
            // calculate needed powers of sideLength to save computing time
            for (unsigned int i = 0; i <= dimensions; ++i) {
               long p = static_cast<unsigned int>(pow(sideLength, i));
               powers.push_back(p);
            }
         }
         
         reference operator*() const { return current; }
         pointer operator->() const { return &current; }

         lattice_iterator& operator++()
         {
            ++dimCount;
            if (dimCount == dimensions) {
               ++source;
               dimCount=0;
               if ((source + 1) % sideLength) {
                  target = source + 1;
               } else if (periodic_boundary) {
                  target = source + 1 - sideLength;
               } else {
                  (*this)++;
               }
            } else {
               target = source + powers[dimCount];
               if (target >= powers[dimCount+1]) {
                  if (periodic_boundary) {
                     target = target - powers[dimCount+1];
                  } else {
                     (*this)++;
                  }
               }
            }
            current.first = source;
            current.second = target;
            return *this;
         }
         
         lattice_iterator operator++(int)
         {
            lattice_iterator temp(*this);
            ++(*this);
            return temp;
         }

         bool operator==(const lattice_iterator& rhs) const
         {
            if (past_the_end && !rhs.past_the_end) {
               return rhs == *this;
            } else if (!past_the_end && rhs.past_the_end) {
               return source == powers[dimensions];
            } else if (past_the_end && rhs.past_the_end) {
               return true;
            }
            return source == rhs.source && target == rhs.target;
         }
         
         bool operator!=(const lattice_iterator& rhs) const
         { return !(*this == rhs); }

         bool past_the_end;
         bool periodic_boundary;
         vertices_size_type sideLength;
         unsigned int dimensions;
         unsigned int dimCount;
         vertices_size_type source;
         vertices_size_type target;
         value_type current;
         std::vector<unsigned long> powers;
   };
   
   /******************************************************************/
   // print_lattice
   // prints the lattice on the screen
   /******************************************************************/
   template <typename Graph>
   void print_lattice
   (Graph& g,
    unsigned int sideLength)
   {
      typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
      vertex_descriptor v;
      for (unsigned int i = 0; i < sideLength; i++) {
         for (unsigned int j = 0; j < sideLength; j++) {
            v = vertex(i*sideLength+j, g);
            std::cout << g[v].state << " ";
         }
         std::cout << std::endl;
      }
   }

}

#endif
