/*! \file community_structure.hh
  \brief Functions to determine the community structure of the networkk
*/

#ifndef COMMUNITY_STRUCTURE_HH
#define COMMUNITY_STRUCTURE_HH

#include <set>
#include <boost/graph/bc_clustering.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
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
    double mod = 0.;

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
    
  template <class Graph, class RandomGenerator>
  double community_structure_rw(const Graph& og, Graph& g,
                                RandomGenerator& r,
                                unsigned int verbose = 0,
                                bool saveGraph = false,
                                bool calcModularity = false,
                                unsigned int numSteps = 0)
  {

    double mod = 0.;

    uniform_real<> uni_dist(0, 1);
    variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);

    if (numSteps == 0) {
      numSteps = num_vertices(og);
    }
    
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    
    boost::multi_array<int, 2> 
      similarity(boost::extents[num_vertices(og)][num_vertices(og)]);
    std::fill(similarity.origin(), similarity.origin()+similarity.size(), 0);
    
    vertex_iterator vi, vi_end;
    vertex_descriptor current_vertex;
    out_edge_iterator oi, oi_end;

    // assign similarities
    if (verbose >= 2) {
      std::cout << "Assigning similarities " << std::endl;
    }
    
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      std::set<unsigned int> C;
      current_vertex = *vi;
      for (unsigned int i = 0; i < numSteps; ++i) {
        double weightSum = 0.;
        for (tie(oi, oi_end) = out_edges(current_vertex, og);
             oi != oi_end; oi++) {
          weightSum += g[*oi].weight;
        }
        double chosenEdge = uni_gen() * weightSum;
        weightSum = 0.;
        bool found = false;
        for (tie(oi, oi_end) = out_edges(current_vertex, og);
             oi != oi_end && !found; oi++) {
          weightSum += g[*oi].weight;
          if (weightSum >= chosenEdge) {
            found = true;
          }
        }
        oi--;
        current_vertex = target(*oi, g);
        C.insert(current_vertex);
      }

      for (std::set<unsigned int>::iterator it = C.begin();
           it != C.end(); it++) {
        for (std::set<unsigned int>::iterator it2 = C.begin();
             it2 != C.end(); it2++) {
          if ((*it2) > (*it)) {
            similarity[*it][*it2]++;
          }
        }
      }
    }

    // identify communities

    if (verbose >= 2) {
      std::cout << "Identifying communities" << std::endl;
    }
    std::vector<std::set<unsigned int> > communities;
    for (unsigned int i = 0; i < num_vertices(og); ++i) {
      std::set<unsigned int> temp;
      temp.insert(i);
      communities.push_back(temp);
    }

    int max = 0;
    int max_i = -1, max_j = -1;
    for (unsigned int i = 0; i < num_vertices(og); ++i) {
      if (similarity[i][0] > -1) {
        for (unsigned int j = i; j < num_vertices(og); ++j) {
          if (similarity[i][j] > max) {
            max = similarity[i][j];
            max_i = i;
            max_j = j;
          }
        }
      }
    }

    while (max > 0) {
      
      if (verbose >= 2) {
        std::cout << "max = " << max << " " << max_i << " " << max_j << std::endl;
      }

      for (unsigned int i = 0; i < num_vertices(og); ++i) {
        if (i != static_cast<unsigned int>(max_i) && 
            similarity[max_i][i] > -1) {
          similarity[max_i][i] =
            (similarity[max_i][i] + similarity[max_j][i])/2;
          similarity[max_j][i] = -1;
          similarity[i][max_j] = -1;
        }
      }
      similarity[max_i][max_j] = -1;

      for (std::set<unsigned int>::iterator it = communities[max_j].begin();
             it != communities[max_j].end(); it++) {
        communities[max_i].insert(*it);
      }

      communities[max_j].clear();

      max = 0;
      for (unsigned int i = 0; i < num_vertices(og); ++i) {
        if (similarity[i][0] > -1) {
          for (unsigned int j = i; j < num_vertices(og); ++j) {
            if (similarity[i][j] > max) {
              max = similarity[i][j];
              max_i = i;
              max_j = j;
            }
          }
        }
      }
    }

    if (verbose >= 1) {
      unsigned int num_communities = 0;
      for (unsigned int i = 0; i < num_vertices(og); ++i) {
        if (communities[i].size() > 0) num_communities++;
      }
      std::cout << "Number of communities: " << num_communities << std::endl;
    }

    for (unsigned int i = 0; i < num_vertices(og); ++i) {
      for (std::set<unsigned int>::iterator it = communities[i].begin();
             it != communities[i].end(); it++) {
        std::cerr << *it << " " << i << std::endl;
      }
    }
    return mod; 
  }
    
} // namespace boost
  
//----------------------------------------------------------
  

#endif

