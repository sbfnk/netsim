/******************************************************************/
// graph_structure.hh
// contains routines to add vertices and edges with a given structure to an
// existing graph 
/******************************************************************/
#ifndef GENERATE_GRAPH_HH
#define GENERATE_GRAPH_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <vector>
#include <iostream>
#include <sys/time.h>

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
  // checks if unseen edge among stubs
  // needed for copy_graph and random_regular_graph

  template <typename Graph>
  bool suitable(Graph& g,
                std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& stubs,
                std::vector<std::vector<bool> >& seen_edges)
  {
    typedef typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::iterator iter;
    
    for (iter s = stubs.begin(); s != stubs.end(); s++) 
      for (iter t = s + 1; t != stubs.end(); t++) 
        if (!seen_edges[*s][*t]) return 1; // success
    
    return 0; //failure
  }
  
  //----------------------------------------------------------
   
  template <typename Graph1, typename Graph2, typename RandomGenerator,
            typename EdgeType>
  int copy_graph(Graph1& source_graph, Graph2& target_graph,
                  RandomGenerator& r,
                  typename edge_property_type<Graph2>::type et =
                  edge_property_type<Graph2>::type(),
                  double rewireFraction = 0., double removeFraction = 0.,
                  double addFraction = 0.)
  {
    typedef typename boost::graph_traits<Graph1>::edge_iterator
      edge_iterator;
    
    typedef typename boost::graph_traits<Graph2>::edge_descriptor
      edge_descriptor;
    typedef typename boost::graph_traits<Graph2>::vertex_descriptor
      vertex_descriptor;

    std::vector< std::vector<bool> >
      seen_edges(num_vertices(source_graph),
                    std::vector<bool>(num_vertices(source_graph), false));
    std::vector<vertex_descriptor> stubs;

    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(source_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, source_graph), target(*ei, source_graph),
               et, target_graph);
      seen_edges[source(*ei, source_graph)][target(*ei, source_graph)] = true;
      seen_edges[target(*ei, source_graph)][source(*ei, source_graph)] = true;
    }

    et.m_value.rewired = true;

    unsigned int N = num_edges(target_graph);
    
    if (rewireFraction > 0.) {
      if (rewireFraction <= 1.) {
        boost::uniform_01<boost::mt19937, double> uni_gen(r);
        unsigned int num_rewire =
          static_cast<unsigned int>(rewireFraction * N);
        // remove num_rewire random edges
        
        for (unsigned int i = 0; i < num_rewire; i++) {
          edge_descriptor e = random_edge(target_graph, r);
          stubs.push_back(source(e, target_graph));
          stubs.push_back(target(e, target_graph));
          boost::remove_edge(e, target_graph);
        }
        // rewire removed edges
        while (num_rewire > 0) {
          // rewire two of the stubs
          unsigned int src =
            static_cast<unsigned int>(uni_gen() * stubs.size());
          unsigned int trg = 
            static_cast<unsigned int>(uni_gen() * stubs.size());
          // check if suitable pair
          if ((src != trg) && (!seen_edges[stubs[src]][stubs[trg]])) {
            boost::add_edge(stubs[src], stubs[trg], et, target_graph);
            // remove stubs
            stubs.erase(stubs.begin() + src);
            if (src > trg) {
              stubs.erase(stubs.begin() + trg);
            } else {
              stubs.erase(stubs.begin() + trg - 1);
            }
            --num_rewire;
          } else {
            if (!suitable(target_graph,stubs, seen_edges))
              return -1; // failure - no more suitable pairs
          }
        }
      } else {
        std::cerr << "ERROR: rewire fraction must be between 0 and 1" << std::endl;
        std::cerr << "no rewiring performed" << std::endl;
      }
    }
    if (removeFraction > 0.) {
      if (removeFraction < 1.) {
        unsigned int num_remove =
          static_cast<unsigned int>(removeFraction * N);
        
        while (num_remove > 0) {
          // select random edge for removal
          edge_descriptor e = random_edge(target_graph, r);

          if (!seen_edges[source(e, target_graph)][target(e, target_graph)]) {
            // remove edge
            boost::remove_edge(e, target_graph);
            --num_remove;
          }
        }
      } else {
        std::cerr << "ERROR: remove fraction must be between 0 and 1" << std::endl;
        std::cerr << "no removing performed" << std::endl;
      }
    } 
    if (addFraction > 0.) {
      unsigned int num_add =
        static_cast<unsigned int>(addFraction * N);

      
      while (num_add > 0) {
        // select random edge for adding
        bool found = false;
        vertex_descriptor v1;
        vertex_descriptor v2;
        while (!found) {
          v1 = random_vertex(target_graph, r);
          do {
            v2 = random_vertex(target_graph, r);
          } while (v2 == v1);
          found = !(edge(v1,v2,target_graph).second);
        }
        add_edge(v1,v2,et,target_graph);
        --num_add;
      }
    }
    return 0;
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

  //----------------------------------------------------------
  
  template <typename Graph, typename DistributionType>
  int random_regular_graph(Graph& g,
                           std::vector<std::pair<unsigned int, unsigned int> >& rrg_edges,
                           const unsigned int d,
                           const unsigned int N,
                           DistributionType& uni_gen)
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;

     
    // define seen_edges NxN zero matrix
    std::vector<std::vector<bool> >
      seen_edges(N, std::vector<bool>(N, false));
    
    // define stubs
    std::vector<vertex_descriptor> stubs;
    
    // init stubs
    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < d; ++j)
        stubs.push_back(i);
    
    // construct random regular graph
    while (stubs.size()) {
      
      // select source and target from stubs
      unsigned int src = static_cast<vertex_descriptor>(uni_gen() * stubs.size());
      unsigned int trg = static_cast<vertex_descriptor>(uni_gen() * stubs.size());
      unsigned int source = stubs[src];
      unsigned int target = stubs[trg];
      
      // check if pair is suitable
      if ( (source != target) && !(seen_edges[source][target]) ) {
        
        // removing source and target from stubs
        stubs.erase(stubs.begin() + src);
        if (src > trg) {
          stubs.erase(stubs.begin() + trg);
        } else {
          stubs.erase(stubs.begin() + trg - 1);
        }
        
        // update seen_edges
        seen_edges[source][target] = true;
        seen_edges[target][source] = true;         
        
        // add edge to rrg_edges
        rrg_edges.push_back(std::make_pair(source, target));
        
      } else { // check if suitable stubs left
        
        if (!suitable(g, stubs, seen_edges)) return 0; // failure - no more suitable pairs
      }
    }   
    
    return 1; // success
  }
  
} // namespace boost

//----------------------------------------------------------
  
#endif
