/******************************************************************/
// graph_structure.hh
// contains routines to add vertices and edges with a given structure to an
// existing graph 
/******************************************************************/
#ifndef GENERATE_GRAPH_HH
#define GENERATE_GRAPH_HH

#include <boost/graph/graph_traits.hpp>

namespace boost {
   template <typename Graph>
   void add_vertices(Graph& g, 
                     typename graph_traits<Graph>::vertices_size_type n)
      {
      for (unsigned int i = 0; i < n; i++) {
         add_vertex(g);
      }
   }

   template <typename Graph>
   void add_vertices(Graph& g, 
                     typename graph_traits<Graph>::vertices_size_type n,
                     const typename vertex_property_type<Graph>::type v)
   {
      for (unsigned int i = 0; i < n; i++) {
         add_vertex(v, g);
      }
   }

   template <typename Graph, typename EdgeIterator>
   void add_edge_structure(Graph& g, EdgeIterator ei, EdgeIterator ei_end)
   {
      while (ei != ei_end) {
         add_edge((*ei).first, (*ei).second, g);
         ++ei;
      }
   }

   template <typename Graph, typename EdgeIterator>
   void add_edge_structure(Graph& g, EdgeIterator ei, EdgeIterator ei_end,
                           const typename edge_property_type<Graph>::type e)
   {
      while (ei != ei_end) {
         add_edge((*ei).first, (*ei).second, e, g);
         ++ei;
      }
   }
}

#endif
