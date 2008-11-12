/*! \file community_structure.hh
  \brief Functions to determine the community structure of the networkk
*/

#ifndef COMMUNITY_STRUCTURE_HH
#define COMMUNITY_STRUCTURE_HH

#include <boost/graph/bc_clustering.hpp>
#include <boost/foreach.hpp>

//! \addtogroup graph_statistics Graph statistics
//! \addtogroup helper_functions Helper functions

//----------------------------------------------------------
namespace boost {

  class clustering_threshold
    : public bc_clustering_threshold<double>
  {
    typedef bc_clustering_threshold<double> inherited;
    
  public:
    template<typename Graph>
    clustering_threshold(double threshold, const Graph& g,
                         bool normalize, unsigned int v = 0)
      : inherited(threshold, g, normalize), iter(1), verbose(v) { }
    
    template<typename Graph, typename Edge>
    bool operator()(double max_centrality, Edge e, const Graph& g)
    {
      std::cout << "Iteration: " << iter << ", max centrality: " 
                << (max_centrality / dividend) << std::endl;
      ++iter;
      return inherited::operator()(max_centrality, e, g);
    }
    
  private:
    unsigned int iter;
    unsigned int verbose;
  };
  
  template <class Graph>
  void community_structure(Graph& g, double threshold)
  {
    typedef typename clustering_threshold::centrality_type centrality_type;
    std::vector<centrality_type> edge_centrality(num_edges(g));
    
    betweenness_centrality_clustering
      (g, clustering_threshold(threshold, g, false),
       make_iterator_property_map(edge_centrality.begin(), get(&Edge::index, g)));
  }
  
} // namespace boost

//----------------------------------------------------------
  

#endif
