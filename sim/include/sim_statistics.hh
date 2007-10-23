/*! \file sim_statistics.hh
  \brief Functions for simulation statistics.
*/

#ifndef SIM_STATISTICS_HH
#define SIM_STATISTICS_HH

#include <fstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/multi_array.hpp>

#include "Model.hh"

//! \addtogroup sim_statistics Simulation statistics

  //----------------------------------------------------------
  /*! \brief Count vertices.
  
  Counts the number of vertices of a given state in a graph

  \param[in] g The graph containing the vertices
  \param[in] nVertexStates The number of vertex states.
  \return A vector of state counts.
  \ingroup sim_statistics
  */
  template <typename Graph>
  std::vector<unsigned int>
  count_vertices(Graph& g, unsigned int nVertexStates)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::vector<unsigned int> counts(nVertexStates,0);

    vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      ++counts[g[*vi].state.base];
    }
   
    return counts;
  }

  //----------------------------------------------------------
  /*! \brief Count state pairs.
  
  Counts the number of pairs of states in a graph.

  \param[in] g The graph containing the vertices and edges
  \param[in] nVertexStates The number of vertex states.
  \param[in] nEdgeTypes The number of edge types.
  \return A 3d array of pairs, the first index of which denounces edge type and
  the other two vertex states
  \ingroup sim_statistics
  */
template <typename Graph>
boost::multi_array<unsigned int, 3>
count_state_pairs(Graph& g, unsigned int nVertexStates,
                  unsigned int nEdgeTypes)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
  typedef boost::multi_array<unsigned int, 3> array_type;

  array_type counts(boost::extents
                    [nEdgeTypes]
                    [nVertexStates]
                    [nVertexStates]);
    
  // count all edges of type et, including parallel
  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    std::vector<unsigned int> states(2);
    states[0] = g[source(*ei, g)].state.base;
    states[1] = g[target(*ei, g)].state.base;
    std::sort(states.begin(), states.end());
    ++counts[g[*ei].type][states[0]][states[1]];
  }

  return counts;
}

//----------------------------------------------------------
/*! \brief Count parallel edges.
  
Loops over all edges and counts for which edges there is another edge
connecting the same two vertices. 

\param[in] g The graph containing the vertices and edges.
\param[in] nVertexStates The number of vertex states.
\return A 2d array of parallel pairs the indices denoting the two vertex states.
\ingroup sim_statistics
*/
template <typename Graph>
boost::multi_array<unsigned int, 2>
count_parallel_edges(Graph& g, unsigned int nVertexStates)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;

  // count only PARALLEL edges
  typedef boost::multi_array<unsigned int, 2> array_type;

  array_type counts(boost::extents
                    [nVertexStates]
                    [nVertexStates]);

  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    if ((g[*ei].type == 0) && g[*ei].parallel) {
      std::vector<unsigned int> states(2);
      states[0] = g[source(*ei, g)].state.base;
      states[1] = g[target(*ei, g)].state.base;
      std::sort(states.begin(), states.end());
      ++counts[states[0]][states[1]];
    }
  }

  // since we looped over all edges of one type, no need to divide by 2
  return counts;
}

//----------------------------------------------------------
/*! \brief Count state triples.
  
Counts the number of triples of states in the network.

\param[in] g The graph containing the vertices.
\param[in] nVertexStates The number of vertex states.
\param[in] nEdgeTypes The number of edge types.
\ingroup sim_statistics
*/
template <typename Graph>
boost::multi_array<unsigned int, 5>
count_state_triples(Graph& g, unsigned int nVertexStates,
                    unsigned int nEdgeTypes)
{
  typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
  typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;

  typedef boost::multi_array<unsigned int, 5> array_type;
    
  array_type counts(boost::extents
                    [nEdgeTypes]
                    [nEdgeTypes]
                    [nVertexStates]
                    [nVertexStates]
                    [nVertexStates]);
    
  vertex_iterator vi, vi_end;
  out_edge_iterator oi, oi_end, oi2;

  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
    if (out_degree(*vi, g) > 0) {
      for (boost::tie(oi, oi_end) = out_edges(*vi, g); oi != oi_end; oi++) {
        for(oi2 = oi + 1; oi2 != oi_end; oi2++) {
          if (target(*oi, g) != target(*oi2, g)) {

            std::vector<unsigned int> edges(2);
            std::vector<unsigned int> states(2);
            edges[0] = g[*oi].type;
            edges[1] = g[*oi2].type;
            states[0] = g[target(*oi, g)].state.base;
            states[1] = g[target(*oi2, g)].state.base;
            std::sort(edges.begin(), edges.end());
            std::sort(states.begin(), states.end());
            
            ++counts
              [edges[0]]
              [edges[1]]
              [g[*vi].state.base]
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
/*! \brief Print simulation status.
  
Prints the number of vertices in all possible states and, if desired, the
number of edges of all types connecting vertices of all possible states, as
well as the number of parallel edges connecting such vertices.

\param[in] g The graph containing the vertices and edges
\param[in] m The model to be used to find all possible states
\param[in] pairs Whether pairs should be counted as well
\param[in] triples Whether triples should be counted as well
\ingroup sim_statistics
*/
template <typename Graph>
void print_sim_status(Graph& g, const Model& m,
                      bool pairs = false, bool triples = false)
{

  unsigned int nVertexStates = m.getVertexStates().size();
  unsigned int nEdgeTypes = m.getEdgeTypes().size();
    
  std::cout << std::endl;
  std::cout << "Vertex count: " << std::endl;

  // count states
  std::vector<unsigned int> vertexCount = count_vertices(g, nVertexStates);
  for (unsigned int i=0; i < nVertexStates; i++) {
    if (vertexCount[i] > 0) {
      std::cout << m.getVertexStates()[i] << ": " << vertexCount[i]
                << std::endl;
    }
  }
   
  std::cout << std::endl;
  if (pairs) {
    
    std::cout << "Pair count: " << std::endl;
   
    // count state pairs
    boost::multi_array<unsigned int, 3> pairCount =
      count_state_pairs(g, nVertexStates, nEdgeTypes);
    for (unsigned int i=0; i < nEdgeTypes; i++) {
      std::cout << m.getEdgeTypes()[i] << "-type: " << std::endl;
      for (unsigned int j=0; j < nVertexStates; j++) {
        for (unsigned int k=j; k < nVertexStates; k++) {
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
      count_parallel_edges(g, nVertexStates);

    for (unsigned int j=0; j < nVertexStates; j++) {
      for (unsigned int k=j; k < nVertexStates; k++) {
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
      count_state_triples(g, nVertexStates, nEdgeTypes);

    for (unsigned int i=0; i < nEdgeTypes; i++) {
      for (unsigned int j=i; j < nEdgeTypes; j++) {
        std::cout << m.getEdgeTypes()[i] << m.getEdgeTypes()[j] << "-triples: "
                  << std::endl;
        for (unsigned int k=0; k < nVertexStates; k++) {
          for (unsigned int l=0; l < nVertexStates; l++) {
            for (unsigned int n=l; n < nVertexStates; n++) {
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
/*! \brief Write simulations data statistics to a file
  
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
\param[in] triples Whether triples should be counted as well (intensive on
computation) or not 
\ingroup sim_statistics
*/
template<typename Graph>
std::string write_sim_data(const Graph& g, const Model& m,
                           double t, std::ofstream& ofile,
                           bool pairs = true, bool triples = true)
{
  unsigned int nVertexStates = m.getVertexStates().size();
  unsigned int nEdgeTypes = m.getEdgeTypes().size();
    
  std::stringstream line("");

  // first in line is current time
  ofile << t << '\t';

  // count vertices' states
  std::vector<unsigned int> vertexCount = count_vertices(g, nVertexStates);
  for (unsigned int i = 0; i < nVertexStates; i++) {
    line << vertexCount[i] << '\t';
  }

  if (pairs) {
    boost::multi_array<unsigned int, 3> pairCount =
      count_state_pairs(g, nVertexStates, nEdgeTypes);
    // count pairs
    for (unsigned int i = 0; i < nEdgeTypes; i++) {
      for (unsigned int j = 0; j < nVertexStates; j++) {
        for (unsigned int k = j; k < nVertexStates; k++) {
          line << pairCount[i][j][k] << '\t';
        }
      }
    }

    // count parallel pairs
    boost::multi_array<unsigned int, 2> parallelCount =
      count_parallel_edges(g, nVertexStates);

    for (unsigned int j = 0; j < nVertexStates; j++) {
      for (unsigned int k = j; k < nVertexStates; k++) {
        line << parallelCount[j][k] << '\t';
      }
    }
  }

  if (triples) {
    boost::multi_array<unsigned int, 5> tripleCount =
      count_state_triples(g, nVertexStates, nEdgeTypes);

    for (unsigned int i=0; i < nEdgeTypes; i++) {
      for (unsigned int j=i; j < nEdgeTypes; j++) {
        for (unsigned int k=0; k < nVertexStates; k++) {
          for (unsigned int l=0; l < nVertexStates; l++) {
            for (unsigned int n=l; n < nVertexStates; n++) {
              line << tripleCount[i][j][k][l][n] << '\t';
            }
          }
        }
      }
    }
  }
  
  line << std::endl;

  ofile << line.str();
  return line.str();
}

//----------------------------------------------------------

#endif
