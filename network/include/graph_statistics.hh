/*! \file graph_statistics.hh
  \brief Functions for collecting graph statistics.
*/

#ifndef GRAPH_STATISTICS_HH
#define GRAPH_STATISTICS_HH

#include <boost/graph/adjacency_list.hpp>
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
  void print_degrees(Graph& g, unsigned int nEdgeTypes)
  {

    std::cout << "print degrees" << std::endl;
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

    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      bool firstLine = true;
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        if (firstLine) {
          std::cout << std::setw(vLength) << *vi;
          firstLine = false;
        } else {
          for (unsigned int j = 0; j < vLength; ++j) std::cout << " ";
        }
        
        std::cout << " : ";
        for (unsigned j = 0; j < out_degree_type(*vi, g, i); ++j) {
          std::cout << symbols[i];
        }
        std::cout << std::endl;
      }
    }
  }

  //----------------------------------------------------------
  /*! \brief Write the degree distribution to a file
  
  Loops over all vertices and records both total degree distribution and degree
  distribution for each edge type. 

  \param[in] g The graph to calculate the degree distribution for.
  \param[in] nEdgeTypes The number of edge types.
  \param[in] degreeFileName The name of the file to write the degree
  distribution to.
  \return 0 if successful, <0 if failed
  \ingroup graph_statistics
  */
  template <typename Graph>
  unsigned int write_degree(const Graph& g, unsigned int nEdgeTypes,
                            const std::string degreeFileName)
  {
    typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
    typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;   

    // counters
    typedef std::map<unsigned int, unsigned int> degree_map;
    std::vector<degree_map*> degree(nEdgeTypes+1); // including parallel
   
    for (unsigned int i = 0; i < degree.size(); i++)
      degree[i] = new degree_map;
   
    // open file
    std::ofstream file;
    try {      
      file.open(degreeFileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open degree file: " << e.what()
                << std::endl;
      return -1;
    }
   
    // loop over all vertices
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      
      // tmp sums
      unsigned int vertex_deg[nEdgeTypes+1];
      for (unsigned int i = 0; i < nEdgeTypes+1; i++) 
        vertex_deg[i] = 0;
      
      // loop over the vertex out edges
      for (tie(ei, ei_end) = out_edges(*vi, g); ei != ei_end; ei++) {
         
        // count all edges including parallel
        vertex_deg[nEdgeTypes] += 1;
         
        // count all other edges (whether they are parallel or not)
        vertex_deg[g[*ei].type] += 1;
      }

      // update global degree
      for (unsigned int i = 0; i < degree.size(); i++)
        (*degree[i])[vertex_deg[i]] += 1;
    }

    // printing
    degree_map::const_iterator it;      
   
    // setting max_degree
    unsigned int max_degree = 0;
    for (unsigned int i =0; i < degree.size(); ++i) {
      it = (*degree[i]).end();
      --it;
      if ((*it).first > max_degree)
        max_degree = (*it).first;
    }

    file << "# max_degree = " << max_degree << std::endl;
    file << "# degree | d-edges | i-edges | total-edges\n";
    file << "# ------ | ------- | ------- | -----------\n";
   
    // enumerate over all degrees and all maps
    for (unsigned int i = 0; i <= max_degree; i++) {
      file << i;
      for (unsigned int j = 0; j < degree.size(); j++)
        file << " " << (*degree[j])[i];
      file << std::endl;
    }    
   
    // close file
    file.close();
   
    // delete
    for (unsigned int i = 0; i <= nEdgeTypes; i++) 
      delete degree[i];
   
    return 0;
   
  }
} // namespace boost

//----------------------------------------------------------

#endif
