/*! \file community_structure.hh
  \brief Functions to determine the community structure of the networkk
*/

#ifndef COMMUNITY_STRUCTURE_HH
#define COMMUNITY_STRUCTURE_HH

#include <boost/graph/bc_clustering.hpp>
#include "Edge.hh"

//! \addtogroup graph_statistics Graph statistics
//! \addtogroup helper_functions Helper functions

//----------------------------------------------------------
namespace boost {


  template<typename MutableGraph, typename Done, typename EdgeCentralityMap,
           typename VertexIndexMap, typename WeightedMap>
  void 
  weighted_betweenness_centrality_clustering(MutableGraph& g, Done done,
                                             EdgeCentralityMap edge_centrality,
                                             VertexIndexMap vertex_index,
                                             WeightedMap edge_weights)
  {
    typedef typename property_traits<EdgeCentralityMap>::value_type
      centrality_type;
    typedef typename graph_traits<MutableGraph>::edge_iterator edge_iterator;
    typedef typename graph_traits<MutableGraph>::edge_descriptor edge_descriptor;
    typedef typename graph_traits<MutableGraph>::vertices_size_type
      vertices_size_type;

    if (has_no_edges(g)) return;

    // Function object that compares the centrality of edges
    indirect_cmp<EdgeCentralityMap, std::less<centrality_type> > 
       cmp(edge_centrality);

    bool is_done;
    do {
      brandes_betweenness_centrality(g, 
                                     edge_centrality_map(edge_centrality)
                                     .vertex_index_map(vertex_index)
                                     .weight_map(edge_weights));
      std::pair<edge_iterator, edge_iterator> edges_iters = edges(g);
      edge_descriptor e = *max_element(edges_iters.first, edges_iters.second, cmp);
      is_done = done(get(edge_centrality, e), e, g);
      if (!is_done) remove_edge(e, g);
    } while (!is_done && !has_no_edges(g));
  }

  /**
   * \overload
   */ 
  template<typename MutableGraph, typename Done, typename EdgeCentralityMap,
           typename WeightedMap>
  void 
  weighted_betweenness_centrality_clustering(MutableGraph& g, Done done,
                                    EdgeCentralityMap edge_centrality,
                                    WeightedMap edge_weights)
  {
    weighted_betweenness_centrality_clustering(g, done, edge_centrality,
                                      get(vertex_index, g), edge_weights);
  }
  
  template<typename Graph>
  class clustering_threshold
    : public bc_clustering_threshold<double>
  {
    typedef bc_clustering_threshold<double> inherited;
    
  public:
    clustering_threshold(double threshold, const Graph& g,
                         bool normalize, unsigned int v = 0,
                         bool pm = false,
                         const Graph* og = 0, Graph* bmg = 0,
                         double* bm = 0)
      : inherited(threshold, g, normalize), iter(1), verbose(v),
        printModularity(pm), original_graph(og), best_mod_graph(bmg),
        bestMod(0.), bestModPtr(bm)
    { }
    
    template<typename Edge>
    bool operator()(double max_centrality, Edge e, const Graph& g)
    {
      if (verbose >= 2) {
        std::cout << "Iteration: " << iter << ", max centrality: " 
                  << (max_centrality / dividend);
      }
      if ((printModularity || best_mod_graph) && original_graph) {
        double mod = graph_modularity(*original_graph, g);
        if (verbose >=2) {
          std::cout << ", modularity: " << mod;
        }
        if (mod > bestMod) {
          bestMod = mod;
          if (bestModPtr) {
            *bestModPtr = bestMod;
          }
          if (best_mod_graph) {
            *best_mod_graph = g;
          }
        }
      }
      if (verbose >=2) {
        std::cout << std::endl;
      }
      ++iter;
      return inherited::operator()(max_centrality, e, g);
    }
    
  private:
    unsigned int iter;
    unsigned int verbose;
    bool printModularity;
    const Graph* original_graph;
    Graph* best_mod_graph;
    double bestMod;
    double* bestModPtr;
  };
  
  template <class Graph>
  double community_structure(const Graph& og, Graph& g,
                             double threshold,
                             unsigned int verbose = 0,
                             bool saveGraph = false,
                             bool calcModularity = false)
  {
    typedef typename clustering_threshold<Graph>::centrality_type
      centrality_type;
    std::vector<centrality_type> edge_centrality(num_edges(g));
    typename boost::property_map<Graph, double Edge::*>::type 
      weight_pmap = get(&Edge::weight, g);

    Graph temp_graph;

    Graph* saveBest = 0;
    Graph* workGraph = &g;
    if (saveGraph) {
      saveBest = &g;
      temp_graph = g;
      workGraph = &temp_graph;
    }

    
    double mod = 0.;
    weighted_betweenness_centrality_clustering
      (*workGraph,
       clustering_threshold<Graph>(threshold, *workGraph, false,
                                   verbose, calcModularity, &og,
                                   saveBest, &mod),
       make_iterator_property_map(edge_centrality.begin(),
                                  get(&Edge::index, g)),
       weight_pmap
       );
    return mod; 
  }
  
} // namespace boost

//----------------------------------------------------------
  

#endif
