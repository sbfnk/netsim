/*! \file community_structure.hh
  \brief Functions to determine the community structure of the networkk
*/

#ifndef COMMUNITY_STRUCTURE_HH
#define COMMUNITY_STRUCTURE_HH

#include <boost/graph/bc_clustering.hpp>

//! \addtogroup graph_statistics Graph statistics
//! \addtogroup helper_functions Helper functions

//----------------------------------------------------------
namespace boost {

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
      if (printModularity && original_graph) {
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

    Graph temp_graph;

    Graph* saveBest = 0;
    Graph* workGraph = &g;
    if (saveGraph) {
      saveBest = &g;
      temp_graph = g;
      workGraph = &temp_graph;
    }

    double mod = 0.;
    betweenness_centrality_clustering
      (*workGraph,
       clustering_threshold<Graph>(threshold, *workGraph, false,
                                   verbose, calcModularity, &og,
                                   saveBest, &mod),
       make_iterator_property_map(edge_centrality.begin(),
                                  get(&Edge::index, g)));
    return mod; 
  }
  
} // namespace boost

//----------------------------------------------------------
  

#endif
