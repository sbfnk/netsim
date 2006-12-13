#include <boost/graph/graph_traits.hpp>
#include <math.h>

namespace boost {

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
      std::vector<unsigned long> powers;
      for (unsigned int i = 0; i <= dimensions; ++i) {
         long p = static_cast<unsigned int>(pow(sideLength, i));
         powers.push_back(p);
      }

      for (v_size_t j = 0; j < powers[dimensions]; ++j) {
         add_vertex(g);
      }
   
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
   

   template <typename MutableGraph>
   void print_lattice
   (MutableGraph& g,
    unsigned int sideLength)
   {
      std::cout << "print_lattice called with sideLength " << sideLength
                << std::endl;
      typedef typename graph_traits<MutableGraph>::vertex_descriptor vertex_descriptor;
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

