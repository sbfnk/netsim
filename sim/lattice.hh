/******************************************************************/
// lattice.hh
// implements a lattice graph for the boost graph library
/******************************************************************/
#include <boost/graph/graph_traits.hpp>
#include <math.h>

namespace boost {

   /******************************************************************/
   // generateLattice
   // creates the lattice graph
   /******************************************************************/
   template <typename MutableGraph>
   void generate_lattice
   (MutableGraph& g, 
    unsigned int dimensions,
    unsigned int sideLength)
   {
      typedef graph_traits<MutableGraph> Traits;
      typedef typename Traits::vertices_size_type v_size_t;
      typedef typename Traits::edges_size_type e_size_t;
      typedef typename Traits::vertex_descriptor vertex_descriptor;

      // calculate needed powers of sideLength to save computing time
      std::vector<unsigned long> powers;
      for (unsigned int i = 0; i <= dimensions; ++i) {
         long p = static_cast<unsigned int>(pow(sideLength, i));
         powers.push_back(p);
      }

      // add vertices
      for (v_size_t j = 0; j < powers[dimensions]; ++j) {
         add_vertex(g);
      }

      // add edges
      for (v_size_t k = 0; k < powers[dimensions]; ++k)
      {
         for (e_size_t l = 1; l <= dimensions; ++l) {
            if (((k%powers[l])+powers[l-1])<powers[l]) {
               add_edge(k, k+powers[l-1], g);
            } else {
               add_edge(k, k+powers[l-1]-powers[l], g);
            }
         }
      }
   }
   

   /******************************************************************/
   // print_lattice
   // displays the lattice on the screen
   /******************************************************************/
   template <typename Graph>
   void print_lattice
   (Graph& g,
    unsigned int sideLength)
   {
      typedef typename graph_traits<Graph>::vertex_descriptor
         vertex_descriptor;
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

