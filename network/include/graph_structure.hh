/******************************************************************/
// graph_structure.hh
// contains routines to add vertices and edges with a given structure to an
// existing graph 
/******************************************************************/
#ifndef GENERATE_GRAPH_HH
#define GENERATE_GRAPH_HH

#include <boost/graph/graph_traits.hpp>

//----------------------------------------------------------

namespace boost {
   template <typename Graph>
   void add_vertices(Graph& g, 
                     typename graph_traits<Graph>::vertices_size_type n)
   {
      for (unsigned int i = 0; i < n; i++) {
         add_vertex(g);
      }
   }
   
   //----------------------------------------------------------
   
   template <typename Graph>
   void add_vertices(Graph& g, 
                     typename graph_traits<Graph>::vertices_size_type n,
                     const typename vertex_property_type<Graph>::type v)
   {
      for (unsigned int i = 0; i < n; i++) {
         add_vertex(v, g);
      }
   }
   
   //----------------------------------------------------------
   
   template <typename Graph, typename EdgeIterator>
   void add_edge_structure(Graph& g, EdgeIterator ei, EdgeIterator ei_end)
   {
      while (ei != ei_end) {
         add_edge((*ei).first, (*ei).second, g);
         ++ei;
      }
   }

   //----------------------------------------------------------
   
   template <typename Graph, typename EdgeIterator>
   void add_edge_structure(Graph& g, EdgeIterator ei, EdgeIterator ei_end,
                           const typename edge_property_type<Graph>::type e)
   {
      while (ei != ei_end) {
         add_edge((*ei).first, (*ei).second, e, g);
         ++ei;
      }
   }

   //----------------------------------------------------------
   
   template <typename Graph>
   unsigned int mark_parallel_edges(Graph& g)
   {
      typename graph_traits<Graph>::edge_iterator ei, ei_end;
      typename graph_traits<Graph>::out_edge_iterator oi, oi_end;
      typename Graph::vertex_descriptor s, t;

      // counter
      unsigned int parallel_edges = 0;      
      
      // loop over all edges
      for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
         // setting source and target vertices
         s = source(*ei, g);
         t = target(*ei, g);
         
         // check if more than one edge is going to the same target
         // from the edge's source
       
         unsigned int parallel_outgoing = 0;
         
         for (tie(oi, oi_end) = out_edges(s, g); oi != oi_end; oi++) {
            if (target(*oi, g) == t) parallel_outgoing++;            
         }
         
         if (parallel_outgoing > 1) {
            
            // count parallel edges
            parallel_edges++;
            
            // mark both parallel edges
            g[*ei].parallel = true;
         }      
      }
      
      // since we looped over all edges, we have counted parallels twice
      parallel_edges = parallel_edges / 2;
      
      return parallel_edges;
      
   }
}

//----------------------------------------------------------
  
#endif
