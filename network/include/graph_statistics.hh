/*! \file graph_statistics.hh
  \brief Functions for collecting graph statistics.
*/

#ifndef GRAPH_STATISTICS_HH
#define GRAPH_STATISTICS_HH

#include <fstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/multi_array.hpp>

#include "Model.hh"

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
  edge(Graph& g,
       typename boost::graph_traits<Graph>::vertex_descriptor u,
       typename boost::graph_traits<Graph>::vertex_descriptor v,
       EdgeType et)
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
  unsigned int out_degree_type(Graph& g, VertexDescriptor v, EdgeType et)
  {
    unsigned int deg = 0;
    
    typename graph_traits<Graph>::out_edge_iterator oi, oi_end;
    for (tie(oi, oi_end) = out_edges(v, g); oi != oi_end; oi++) {
      if (g[*oi].type == et) ++deg;
    }
    return deg;
  }
  
  //----------------------------------------------------------
  /*! \brief Count triangles.
  
  Counts the number of triangles in the network combining three given edge
  types. 

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
            if (edge(g, source(*ei, g), target(*oi, g), et3).second) {
              ++count;
            }
          }
        }
      }
    }
    
    return count;
  }

  //----------------------------------------------------------
  /*! \brief Count vertices.
  
  Counts the number of vertices of a given state in a graph

  \param[in] g The graph containing the vertices
  \param[in] m The model providing the vertex states
  \return A vector of state counts.
  \ingroup graph_statistics
  */
  template <typename Graph>
  std::vector<unsigned int>
  count_vertices(Graph& g, const Model& m)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::vector<unsigned int> counts(m.getVertexStates().size(),0);

    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      ++counts[g[*vi].state];
    }
   
    return counts;
  }

  //----------------------------------------------------------
  /*! \brief Count pairs.
  
  Counts the number of pairs in a given state in a graph.

  \param[in] g The graph containing the vertices and edges
  \param[in] m The model providing the vertex states and edge types
  \return A 3d array of pairs, the first index of which denounces edge type and
  the other two vertex states
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<unsigned int, 3>
  count_pairs(Graph& g, const Model& m)
  {
    typedef typename boost::graph_traits<Graph>::edge_iterator
      edge_iterator;
    typedef boost::multi_array<unsigned int, 3> array_type;

    array_type counts(boost::extents
                      [m.getEdgeTypes().size()]
                      [m.getVertexStates().size()]
                      [m.getVertexStates().size()]);
    
    // count all edges of type et, including parallel
    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      std::vector<unsigned int> states(2);
      states[0] = g[source(*ei, g)].state;
      states[1] = g[target(*ei, g)].state;
      std::sort(states.begin(), states.end());
      ++counts[g[*ei].type][states[0]][states[1]];
    }

    return counts;
  }

  //----------------------------------------------------------
  /*! \brief Count parallel edges.
  
  Loops over all edges and counts for which edges there is another edge
  connecting the same two vertices. 

  \param[in] g The graph containing the vertices and edges
  \param[in] m The model providing the vertex states and edge types
  \return A 2d array of parallel pairs the indices denoting the two vertex states
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<unsigned int, 2>
  count_parallel_edges(Graph& g, const Model& m)
  {
    typedef typename boost::graph_traits<Graph>::edge_iterator
      edge_iterator;

    // count only PARALLEL edges
    typedef boost::multi_array<unsigned int, 2> array_type;

    array_type counts(boost::extents
                      [m.getVertexStates().size()]
                      [m.getVertexStates().size()]);

    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      if ((g[*ei].type == 0) && g[*ei].parallel) {
        std::vector<unsigned int> states(2);
        states[0] = g[source(*ei, g)].state;
        states[1] = g[target(*ei, g)].state;
        std::sort(states.begin(), states.end());
        ++counts[states[0]][states[1]];
      }
    }

    // since we looped over all edges of one type, no need to divide by 2
    return counts;
  }

  //----------------------------------------------------------
  /*! \brief Count triples.
  
  Counts the number of triples in the network.

  \param[in] g The graph containing the vertices
  \param[in] m The model providing the vertex states
  \ingroup graph_statistics
  */
  template <typename Graph>
  boost::multi_array<unsigned int, 5>
  count_triples(Graph& g, const Model& m)
  {
    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;
    typedef typename graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    typedef boost::multi_array<unsigned int, 5> array_type;
    
    array_type counts(boost::extents
                      [m.getEdgeTypes().size()]
                      [m.getEdgeTypes().size()]
                      [m.getVertexStates().size()]
                      [m.getVertexStates().size()]
                      [m.getVertexStates().size()]);
    
    vertex_iterator vi, vi_end;
    out_edge_iterator oi, oi_end, oi2;

    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      if (out_degree(*vi, g) > 0) {
        for (tie(oi, oi_end) = out_edges(*vi, g); oi != oi_end; oi++) {
          for(oi2 = oi + 1; oi2 != oi_end; oi2++) {
            if (target(*oi, g) != target(*oi2, g)) {

              std::vector<unsigned int> edges(2);
              std::vector<unsigned int> states(2);
              edges[0] = g[*oi].type;
              edges[1] = g[*oi2].type;
              states[0] = g[target(*oi, g)].state;
              states[1] = g[target(*oi2, g)].state;
              std::sort(edges.begin(), edges.end());
              std::sort(states.begin(), states.end());
            
              ++counts
                [edges[0]]
                [edges[1]]
                [g[*vi].state]
                [states[0]]
                [states[1]];

            }
          }
        }
      }
    }

    return counts;
    
  }

  //----------------------------------------------------------
  /*! \brief Print graph statistics.
  
  Prints the number of vertices in all possible states and, if desired, the
  number of edges of all types connecting vertices of all possible states, as
  well as the number of parallel edges connecting such vertices.

  \param[in] g The graph containing the vertices and edges
  \param[in] m The model to be used to find all possible states
  \param[in] pairs Whether pairs should be counted as well
  \param[in] triples Whether triples should be counted as well
  \ingroup graph_statistics
  */
  template <typename Graph>
  void print_graph_statistics(Graph& g, Model& m, bool pairs = true,
                              bool triples = false)
  {
    std::cout << std::endl;
    std::cout << "Vertex count: " << std::endl;

    // count states
    std::vector<unsigned int> vertexCount = count_vertices(g, m);
    for (unsigned int i=0; i < m.getVertexStates().size(); i++) {
      if (vertexCount[i] > 0) {
        std::cout << m.getVertexStates()[i] << ": " << vertexCount[i]
                  << std::endl;
      }
    }
   
    std::cout << std::endl;
    if (pairs) {
    
      std::cout << "Pair count: " << std::endl;
   
      // count state pairs
      boost::multi_array<unsigned int, 3> pairCount = count_pairs(g, m);
      for (unsigned int i=0; i < m.getEdgeTypes().size(); i++) {
        std::cout << m.getEdgeTypes()[i] << "-type: " << std::endl;
        for (unsigned int j=0; j < m.getVertexStates().size(); j++) {
          for (unsigned int k=j; k < m.getVertexStates().size(); k++) {
            if (pairCount[i][j][k] > 0) {
              std::cout << m.getVertexStates()[j] << m.getVertexStates()[k]
                        << ": " << pairCount[i][j][k] << std::endl;
            }
          }
        }
      }
   
      // count parallel pairs
      std::cout << "parallel:" << std::endl;

      boost::multi_array<unsigned int, 2> parallelCount =
        count_parallel_edges(g, m);

      for (unsigned int j=0; j < m.getVertexStates().size(); j++) {
        for (unsigned int k=j; k < m.getVertexStates().size(); k++) {
          if (parallelCount[j][k] > 0) {
            std::cout << m.getVertexStates()[j] << m.getVertexStates()[k]
                      << ": " << parallelCount[j][k] << std::endl;
          }
        }
      }
      std::cout << std::endl;
    }

    if (triples) {
      std::cout << "Triple count: " << std::endl;

      boost::multi_array<unsigned int, 5> tripleCount =
        count_triples(g, m);

      for (unsigned int i=0; i < m.getEdgeTypes().size(); i++) {
        for (unsigned int j=i; j < m.getEdgeTypes().size(); j++) {
          std::cout << m.getEdgeTypes()[i] << m.getEdgeTypes()[j]
                    << "-triples: " << std::endl;
          for (unsigned int k=0; k < m.getVertexStates().size(); k++) {
            for (unsigned int l=0; l < m.getVertexStates().size(); l++) {
              for (unsigned int n=l; n < m.getVertexStates().size(); n++) {
                std::cout << m.getVertexStates()[l] << m.getVertexStates()[k]
                          << m.getVertexStates()[n] << ": "
                          << tripleCount[i][j][k][l][n] << std::endl;
              }
            }
          }
        }
      }
      std::cout << std::endl;
    }
  }

  //----------------------------------------------------------
  /*! \brief Write graph statistics to a file
  
  Writes a line containing the current time (passed as a parameter) to a file,
  as well as the number of vertices in all possible states and, if desired, the 
  number of edges of all types connecting vertices of all possible states, as
  well as the number of parallel edges connecting such vertices.

  \param[in] g The graph containing the vertices and edges
  \param[in] m The model to be used to find all possible states
  \param[in] t The current time to be written to the first column
  \param[in] ofile A stream to the output file
  \param[in] pairs Whether pairs should be counted as well (intensive on
  computation) or not 
  \ingroup graph_statistics
  */
  template<typename Graph>
  std::string write_graph_data(const Graph& g, const Model& m,
                               double t, std::ofstream& ofile, bool pairs = true)
  {
    std::stringstream line("");

    // first in line is current time
    ofile << t << '\t';

    // count vertices' states
    std::vector<unsigned int> vertexCount = count_vertices(g, m);
    for (unsigned int i = 0; i < m.getVertexStates().size(); i++) {
      line << vertexCount[i] << '\t';
    }

    if (pairs) {
      boost::multi_array<unsigned int, 3> pairCount = count_pairs(g, m);
      // count pairs
      for (unsigned int i = 0; i < m.getEdgeTypes().size(); i++) {
        for (unsigned int j = 0; j < m.getVertexStates().size(); j++) {
          for (unsigned int k = j; k < m.getVertexStates().size(); k++) {
            line << pairCount[j][k][i] << '\t';
          }
        }
      }

      // count parallel pairs
      boost::multi_array<unsigned int, 2> parallelCount =
        count_parallel_edges(g, m);

      for (unsigned int j = 0; j < m.getVertexStates().size(); j++) {
        for (unsigned int k = j; k < m.getVertexStates().size(); k++) {
          line << parallelCount[j][k] << '\t';
        }
      }
    }
  
    line << std::endl;

    ofile << line.str();
    return line.str();
  }

  //----------------------------------------------------------
  /*! \brief Print degrees

  Prints the degree of each vertex

  \param[in] g The graph containing the vertices and edges
  \param[in] m The model to be used to find all possible edge types
  \ingroup graph_statistics
  */
  template <typename Graph>
  void print_degrees(Graph& g, Model& m)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::vector<std::string> symbols;
    symbols.push_back("*");
    symbols.push_back("-");

    std::stringstream s;
    s << num_vertices(g);
    
    unsigned int vLength = s.str().length();

    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      bool firstLine = true;
      for (unsigned int i = 0; i < m.getEdgeTypes().size(); ++i) {
        if (firstLine) {
          std::cout << std::setw(vLength) << *vi;
          firstLine = false;
        } else {
          for (unsigned int j = 0; j < vLength; ++j) std::cout << " ";
        }
        
        std::cout << " : ";
        for (unsigned j = 0; j < out_degree_type(g, *vi, i); ++j) {
          std::cout << symbols[i];
        }
        std::cout << std::endl;
      }
    }
  }

} // namespace boost

//----------------------------------------------------------

#endif
