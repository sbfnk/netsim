/*! \file graph_statistics.hh
  \brief Functions for collecting graph statistics.
*/

#ifndef GRAPH_STATISTICS_HH
#define GRAPH_STATISTICS_HH

#include <math.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/multi_array.hpp>

//! \addtogroup graph_statistics Graph statistics

namespace boost {

  //----------------------------------------------------------
  /*! \brief Checks for target and type of an edge.
    \ingroup helper_functions
  */
  template <class Graph, class EdgeType>
  struct target_and_type {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    target_and_type(Graph& g, Vertex v, EdgeType et) :
      m_graph(g), m_target(v), m_edgetype(et) { }
    template <class StoredEdge>
    bool operator()(const StoredEdge& e) const {
      return (target(e, m_graph) == m_target) &&
        (m_graph[e].type == m_edgetype);
    }
    Graph& m_graph;
    Vertex m_target;
    EdgeType m_edgetype;
  };

  //----------------------------------------------------------
  /*! \brief Check if there is an edge of give type between two vertices.
    \param[in] g The graph to consider.
    \param[in] u The first vertex to consider.
    \param[in] v The second vertex to consider.
    \param[in] et The edge type to consider.
    \ingroup helper_functions
  */
  template <typename Graph, typename EdgeType>
  std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool>
  edge(typename boost::graph_traits<Graph>::vertex_descriptor u,
       typename boost::graph_traits<Graph>::vertex_descriptor v,
       Graph& g, EdgeType et)
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;

    out_edge_iterator oi, oi_end;
    tie(oi, oi_end) = out_edges(u, g);
    
    bool found;
    out_edge_iterator
      i = std::find_if(oi, oi_end,
                       target_and_type<Graph, EdgeType>(g, v, et));
    found = (i != oi_end);
    if (found)
      return std::make_pair(edge_descriptor(u, v, (*i).get_property()),
                            true);
    else
      return std::make_pair(edge_descriptor(u, v, 0), false);
  }
  
  //----------------------------------------------------------
  /*! \brief Calculate out_degree of a given edge type for a given vertex.
    \param[in] g The graph to consider.
    \param[in] v The vertex to calculate the degree of.
    \param[in] et The edge type to consider.
    \ingroup helper_functions
  */
  template <typename Graph, typename VertexDescriptor, typename EdgeType>
  unsigned int out_degree_type(VertexDescriptor v, Graph& g, EdgeType et)
  {
    unsigned int deg = 0;
    
    typename graph_traits<Graph>::out_edge_iterator oi, oi_end;
    for (tie(oi, oi_end) = out_edges(v, g); oi != oi_end; oi++) {
      if (g[*oi].type == et) ++deg;
    }
    return deg;
  }
  
  //----------------------------------------------------------
  /*! \brief Counts the number of triangles a vertex takes part in.

  Calculates the number of triangles a vertex takes part in, as well as the
  number of triples with the vertex in the middle. This information can be used
  to calculate clustering coefficients.
  
  \param[in] g The graph containing the vertices.
  \param[in] nEdgeTypes The number of edge types in the graph.
  \param[in] v The vertex under consideration.
  \return A pair of multi_arrays, the first one containing the number of triples
  for a given combination of edge types, the second the number of triangles.
  \ingroup graph_statistics
  */
  template <typename Graph>
  std::pair<boost::multi_array<unsigned int, 2>,
            boost::multi_array<unsigned int, 3> >
  count_triangles_vertex
  (Graph& g, unsigned int nEdgeTypes,
   typename boost::graph_traits<Graph>::vertex_descriptor v)
  {

    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;

    typedef boost::multi_array<unsigned int, 2> triple_type;
    typedef boost::multi_array<unsigned int, 3> triangle_type;
    
    triple_type triples(boost::extents
                       [nEdgeTypes]
                       [nEdgeTypes]);
    
    triangle_type triangles(boost::extents
                         [nEdgeTypes]
                         [nEdgeTypes]
                         [nEdgeTypes]);
    
    out_edge_iterator oi, oi_end, oi2, oi3, oi3_end;

    for (tie(oi, oi_end) = out_edges(v, g); oi != oi_end; oi++) {
      for(oi2 = oi + 1; oi2 != oi_end; oi2++) {
        // make sure we do not have two edges going to the same vertex
        if (target(*oi, g) != target(*oi2, g)) {
          std::vector<unsigned int> edges(3);
          edges[0] = g[*oi].type;
          edges[1] = g[*oi2].type;
          std::sort(edges.begin(), edges.begin()+1);
          ++triples
            [edges[0]]
            [edges[1]];
          for (tie(oi3, oi3_end) = out_edges(target(*oi, g), g);
               oi3 != oi3_end; oi3++) {
            if (target(*oi3, g) == target(*oi2, g)) {
              edges[2] = g[*oi3].type;
              ++triangles
                [edges[0]]
                [edges[1]]
                [edges[2]];
            }
          }
        }
      }
    }

    return std::make_pair(triples, triangles);
  }

  //----------------------------------------------------------
  /*! \brief Calculate local clustering coefficients.

  Calculates the clustering coefficients for a network using the averaged
  clustering coefficients of the vertices, as defined by Watts and Strogatz
  (1998). The clustering coefficients are printed to stdout.
  
  \param[in] g The graph containing the vertices.
  \param[in] nEdgeTypes The number of edge types in the graph.
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<double, 3>
  local_cluster_coeffs
  (Graph& g, unsigned int nEdgeTypes)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    typedef boost::multi_array<unsigned int, 2> triple_type;
    typedef boost::multi_array<unsigned int, 3> triangle_type;
    typedef boost::multi_array<double, 3> coeff_type;
    coeff_type cluster_coeffs(boost::extents
                              [nEdgeTypes]
                              [nEdgeTypes]
                              [nEdgeTypes]);

    for (unsigned int i = 0; i< nEdgeTypes; ++i) {
      for (unsigned int j = i; j< nEdgeTypes; ++j) {
        for (unsigned int k = 0; k< nEdgeTypes; ++k) {
          cluster_coeffs[i][j][k] = 0.;
        }
      }
    }
    
    unsigned int vertex_count = 0;
    
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      if (out_degree(*vi, g) > 1) {
        ++vertex_count;
        std::pair<triple_type, triangle_type> vertex_triangles =
          count_triangles_vertex(g, nEdgeTypes, *vi);
        for (unsigned int i = 0; i< nEdgeTypes; ++i) {
          for (unsigned int j = i; j< nEdgeTypes; ++j) {
            if (vertex_triangles.first[i][j] > 0) {
              for (unsigned int k = 0; k< nEdgeTypes; ++k) {
                cluster_coeffs[i][j][k] +=
                  static_cast<double>(vertex_triangles.second[i][j][k]) /
                  static_cast<double>(vertex_triangles.first[i][j]);
              }
            }
          }
        }
      }
    }

    if (vertex_count > 0) {
      for (unsigned int i = 0; i< nEdgeTypes; ++i) {
        for (unsigned int j = i; j< nEdgeTypes; ++j) {
          for (unsigned int k = 0; k< nEdgeTypes; ++k) {
            cluster_coeffs[i][j][k] /= static_cast<double>(vertex_count);
          }
        }
      }
    }

    return cluster_coeffs;
  }
  
  //----------------------------------------------------------
  /*! \brief Calculate global clustering coefficients.

  Calculates the clustering coefficients for a network using the ratio of
  triangles to triples in the graph. The clustering coefficients are printed to
  stdout. 
  
  \param[in] g The graph containing the vertices.
  \param[in] nEdgeTypes The number of edge types in the graph.
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<double, 3>
  global_cluster_coeffs
  (Graph& g, unsigned int nEdgeTypes)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;

    typedef boost::multi_array<unsigned int, 2> triple_type;
    typedef boost::multi_array<unsigned int, 3> triangle_type;
    typedef boost::multi_array<double, 3> coeff_type;
    
    triple_type triples(boost::extents
                       [nEdgeTypes]
                       [nEdgeTypes]);
    
    triangle_type triangles(boost::extents
                         [nEdgeTypes]
                         [nEdgeTypes]
                         [nEdgeTypes]);
    
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      if (out_degree(*vi, g) > 1) {
        std::pair<triple_type, triangle_type> vertex_triangles =
          count_triangles_vertex(g, nEdgeTypes, *vi);
        for (unsigned int i = 0; i< nEdgeTypes; ++i) {
          for (unsigned int j = i; j< nEdgeTypes; ++j) {
            triples[i][j] += vertex_triangles.first[i][j];
            for (unsigned int k = 0; k< nEdgeTypes; ++k) {
              triangles[i][j][k] += vertex_triangles.second[i][j][k];
            }
          }
        }
      }
    }
    
    coeff_type cluster_coeffs(boost::extents
                              [nEdgeTypes]
                              [nEdgeTypes]
                              [nEdgeTypes]);

    for (unsigned int i = 0; i< nEdgeTypes; ++i) {
      for (unsigned int j = i; j< nEdgeTypes; ++j) {
        for (unsigned int k = 0; k< nEdgeTypes; ++k) {
          if (triples[i][j] > 0) {
            cluster_coeffs[i][j][k] =
              static_cast<double>(triangles[i][j][k]) /
              static_cast<double>(triples[i][j]);
          } else {
            cluster_coeffs[i][j][k] = 0.;
          }
        }
      }
    }

    return cluster_coeffs;
  }
  
  //----------------------------------------------------------
  /*! \brief Count triangles.
  
  Counts the number of triangles in the graph.

  \param[in] g The graph containing the vertices
  \param[in] et1 The first edge type to consider
  \param[in] et2 The second edge type to consider
  \param[in] et3 The third edge type to consider
  \return The number of triangles found.
  \ingroup graph_statistics
  */
  template <typename Graph, typename EdgeType>
  unsigned int count_triangles(Graph& g, EdgeType et1, EdgeType et2,
                               EdgeType et3)
  {
    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;

    unsigned int count = 0;
    
    edge_iterator ei, ei_end;
    out_edge_iterator oi, oi_end;

    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      if (g[*ei].type == et1) {
        for (tie(oi, oi_end) = out_edges(target(*ei, g), g);
             oi != oi_end; oi++) {
          if (g[*oi].type == et2) {
            if (edge(source(*ei, g), target(*oi, g), et3, g).second) {
              ++count;
            }
          }
        }
      }
    }
    
    return count;
  }


  //----------------------------------------------------------
  /*! \brief Count triples centered on a given vertex.
  \param[in] g The graph containing the vertices.
  \param[in] nEdgeTypes The number of edge types in the graph.
  \param[in] v The vertex under consideration.
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<unsigned int, 2>
  count_triples_vertex
  (Graph& g, unsigned int nEdgeTypes,
   typename boost::graph_traits<Graph>::vertex_descriptor v)
  {
    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    typedef boost::multi_array<unsigned int, 2> array_type;
    
    array_type counts(boost::extents
                      [nEdgeTypes]
                      [nEdgeTypes]);
    
    out_edge_iterator oi, oi_end, oi2;

    if (out_degree(v, g) > 1) {
      for (tie(oi, oi_end) = out_edges(v, g); oi != oi_end; oi++) {
        for(oi2 = oi + 1; oi2 != oi_end; oi2++) {
          if (target(*oi, g) != target(*oi2, g)) {
            
            std::vector<unsigned int> edges(2);
            edges[0] = g[*oi].type;
            edges[1] = g[*oi2].type;
            std::sort(edges.begin(), edges.end());
            
            ++counts
              [edges[0]]
              [edges[1]];
          }
        }
      }
    }

    return counts;
  }
  
  //----------------------------------------------------------
  /*! \brief Count triples in a graph.
  \param[in] g The graph containing the vertices.
  \param[in] nEdgeTypes The number of edge types in the graph.
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<unsigned int, 2>
  count_triples(Graph& g, unsigned int nEdgeTypes)
  {
    typedef typename graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    typedef boost::multi_array<unsigned int, 2> array_type;
    
    array_type counts(boost::extents
                      [nEdgeTypes]
                      [nEdgeTypes]);

    vertex_iterator vi, vi_end;
    
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      array_type vertex_counts =
        count_triples_vertex(g, nEdgeTypes, *vi);
      for (unsigned int i = 0; i< nEdgeTypes; ++i) {
        for (unsigned int j = i; j < nEdgeTypes; ++j) {
          counts[i][j] += vertex_counts[i][j];
        }
      }
    }

    return counts;
  }

  //----------------------------------------------------------
  /*! \brief Count pairs.
  
  Counts the number of pairs in a graph.

  \param[in] g The graph containing the vertices and edges
  \param[in] nEdgeTypes The number of edge types in the graph.
  \return A vector of counts of the number of pairs of a given edge type
  \ingroup graph_statistics
  */
  template <typename Graph>
  std::vector<unsigned int>
  count_pairs(Graph& g, unsigned int nEdgeTypes)
  {
    typedef typename boost::graph_traits<Graph>::edge_iterator
      edge_iterator;

    std::vector<unsigned int> counts(nEdgeTypes, 0);
    
    // count all edges of type et, including parallel
    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      ++counts[g[*ei].type];
    }

    return counts;
  }


  //----------------------------------------------------------
  /*! \brief Print degrees

  Prints the degree of each vertex

  \param[in] g The graph containing the vertices and edges
  \param[in] nEdgeTypes The number of edge types in the graph.
  \ingroup graph_statistics
  */
  template <typename Graph>
  void print_degrees(Graph& g, unsigned int nEdgeTypes,
                     bool print = true, std::string baseFileName = "")
  {

    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::vector<std::string> symbols;
    symbols.push_back("*");
    symbols.push_back("%");
    symbols.push_back("@");
    symbols.push_back("&");

    std::stringstream s;
    s << num_vertices(g);
    
    unsigned int vLength = s.str().length();

    std::vector< std::vector<unsigned int> > degrees;
    degrees.resize(num_vertices(g));

    unsigned int maxDegree = 0;
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        unsigned int deg = out_degree_type(*vi, g, i);
        degrees[*vi].push_back(deg);
        if (deg > maxDegree) {
	  maxDegree = deg;
	}
      }
    }
    s.str("");
    s << maxDegree;
    unsigned int mdLength = s.str().length();

    s.str("");
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      bool firstLine = true;
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        if (firstLine) {
          s << std::setw(vLength) << *vi;
          firstLine = false;
        } else {
          for (unsigned int j = 0; j < vLength; ++j) s << " ";
        }
        
        s << " : ";
	s << std::setw(mdLength) << degrees[*vi][i] << " ";

        for (unsigned j = 0; j < degrees[*vi][i]; ++j) {
          s << symbols[i];
        }
        s << std::endl;
      }
    }

    if (print) {
      std::cout << s.str();
    }
    if (baseFileName.length() > 0) {
      std::string vDegName = baseFileName + ".stat.vdegree";
      std::ofstream vDegFile(vDegName.c_str(), std::ios::out);
      vDegFile << s.str();
      vDegFile.close();
    }    
  }

  //----------------------------------------------------------
  /*! \brief Write the degree distribution to a file
  
  Loops over all vertices and records both total degree distribution and degree
  distribution for each edge type. 

  \param[in] g The graph to calculate the degree distribution for.
  \param[in] nEdgeTypes The number of edge types.
  \param[out] degrees 2D array of edgetypes and corresponding degree 
  distributions, the first (nEdgeTypes) entries for the degree distribution of 
  edges of that type only (non-parallel), the next (nEdgeTypes) entries for all
  edges of that type (including parallel), the next for parallel and the last 
  for all edges
  \return the maximal degree
  \ingroup graph_statistics
  */
  template <typename Graph>
  unsigned int degree_dist(const Graph& g, unsigned int nEdgeTypes,
                           boost::multi_array<unsigned int, 2>& degrees)
  {
    typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
    typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;   

    unsigned int max_degree = 0;
    
    // counters
    degrees.resize(boost::extents[nEdgeTypes*2+2][max_degree+1]);
    for (unsigned int i = 0; i < nEdgeTypes*2+2; ++i) {
      degrees[i][0] = 0;
    }

    // loop over all vertices
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      
      // tmp sums
      std::vector<unsigned int> vertex_deg(nEdgeTypes*2+2, 0);
      
      // loop over the vertex out edges
      for (tie(ei, ei_end) = out_edges(*vi, g); ei != ei_end; ei++) {
         
        // count all edges including parallel
        if (g[*ei].parallel) {
          // parallel edges only
          ++vertex_deg[nEdgeTypes*2];
        } else {
          // all non-parallel edges of the type
          ++vertex_deg[g[*ei].type];
        }
        // all edges of the type including parallel
        ++vertex_deg[nEdgeTypes+g[*ei].type];
        // all edges
        ++vertex_deg[nEdgeTypes*2 + 1];
      }

      // count parallel edges only once
      vertex_deg[nEdgeTypes*2] /= 2;
      vertex_deg[nEdgeTypes*2+1] -= vertex_deg[nEdgeTypes*2];
                 
      // update max_degree
      if (vertex_deg[nEdgeTypes*2 + 1] > max_degree) {
        max_degree = vertex_deg[nEdgeTypes*2 + 1];
        degrees.resize(boost::extents[nEdgeTypes*2+2][max_degree+1]);
      }
      // update global degree
      for (unsigned int i = 0; i < nEdgeTypes*2+2; i++) {
        ++degrees[i][vertex_deg[i]];
      }
    }
           
    return max_degree;
  }
    
  //----------------------------------------------------------
  /*! \brief Calculate the graph modularity
  
  Loops over all vertices and calculates the graph modularity
  given a graph and a graph indicating communities

  \param[in] g The graph to calculate the degree distribution for.
  \param[in] cg The graph indicating the communities
  \return the graph modularity
  \ingroup graph_statistics
  */

  template <typename Graph>
  double graph_modularity(const Graph& g, const Graph& cg)
  {
    typedef typename graph_traits<Graph>::edge_descriptor
      edge_descriptor;

    std::vector<int> component(num_vertices(cg));
    int num = boost::connected_components(cg, &component[0]);
    std::size_t nEdges = num_edges(g);

    boost::multi_array<std::size_t, 2> community_matrix
      (boost::extents[num][num]);

    for (int i = 0; i < num; ++i) {
      for (int j = 0; j < num; ++j) {
        community_matrix[i][j] = 0;
      }
    }

    BOOST_FOREACH(edge_descriptor e, edges(g)) {
      ++community_matrix
        [component[source(e, g)]]
        [component[target(e, g)]];
      if (component[source(e, g)] == component[target(e, g)]) {
        std::cout << "Inner link: " << source(e, g) << "--"
                  << target(e, g) << std::endl;
      }
    }

    unsigned int trace = 0;
    unsigned int rowSqSum = 0;
      
    for (int i = 0; i < num; ++i) {
      std::cout << "Inner links, group " << i << ": " << community_matrix[i][i]
                << std::endl;
      trace += community_matrix[i][i];
      std::size_t tempSum = 0;
      for (int j = 0; j < num; ++j) {
        tempSum += community_matrix[i][j];
      }
      std::cout << "All links, group " << i << ": " << tempSum << std::endl;
      rowSqSum += (community_matrix[i][i] + tempSum) *
        (community_matrix[i][i] + tempSum);
    }
    
    return
      (trace / static_cast<double>(nEdges)) -
      (rowSqSum / static_cast<double>(4* nEdges * nEdges));
  }
} // namespace boost

//----------------------------------------------------------

#endif
