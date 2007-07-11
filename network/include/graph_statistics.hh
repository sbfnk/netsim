/*! \file graph_statistics.hh
  \brief Functions for collecting graph statistics.
*/

#ifndef GRAPH_STATISTICS_HH
#define GRAPH_STATISTICS_HH

#include <fstream>

#include <boost/graph/adjacency_list.hpp>

#include "Model.hh"

//! \addtogroup graph_statistics Graph statistics

//----------------------------------------------------------
/*! \brief Count vertices.
  
Counts the number of vertices in a graph being in a given state 

\param[in] g The graph containing the vertices
\param[in] vs The vertex state to be counted
\return The number of vertices found in the given state
\ingroup graph_statistics
*/
template <typename Graph>
unsigned int count_vertices(Graph& g, unsigned int vs)
{
  typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;

  unsigned int count = 0;
  vertex_iterator vi, vi_end;
  for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
    if (g[*vi].state == vs) {
      count++;
    }
  }
   
  return count;
}

//----------------------------------------------------------
/*! \brief Count edges.
  
Counts the number of edges of a given type connecting two vertices in a
given state in a graph.

\param[in] g The graph containing the vertices and edges
\param[in] vs1 The state the vertex on one end of the edge is supposed to be in
\param[in] vs2 The state the vertex on the other end of the edge is supposed
to be in
\param[in] et The type the edge is supposed to be of
\return The number of edges fulfilling the given criteria
\ingroup graph_statistics
*/
template <typename Graph>
unsigned int count_edges(Graph& g,
                         unsigned int vs1, unsigned int vs2, unsigned int et)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;

  // count all edges of type et, including parallel
  unsigned int count = 0;
  edge_iterator ei, ei_end;
  for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    if ((g[*ei].type == et) && //!g[*ei].parallel &&
        ((g[source(*ei, g)].state == vs1 && g[target(*ei, g)].state == vs2) ||
         (g[source(*ei, g)].state == vs2 && g[target(*ei, g)].state == vs1))) {
      count++;
    }
  }
   
  return count;
}

//----------------------------------------------------------
/*! \brief Count parallel edges.
  
Loops over all edges of a given type, connecting two vertices in given states,
and counts for which edges there is another edge connecting the same two
vertices. 

\param[in] g The graph containing the vertices and edges
\param[in] vs1 The state the vertex on one end of the edge is supposed to be in
\param[in] vs2 The state the vertex on the other end of the edge is supposed
to be in
\param[in] et The type the edge is supposed to be of
\return The number of edges fulfilling the given criteria and having a
parallel edge
\ingroup graph_statistics
*/
template <typename Graph>
unsigned int count_parallel_edges(Graph& g, unsigned int vs1,
                                  unsigned int vs2, unsigned int et)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;

  // count only PARALLEL edges of type et
  unsigned int count = 0;
  edge_iterator ei, ei_end;
  for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    if ((g[*ei].type == et) && g[*ei].parallel &&
        ((g[source(*ei, g)].state == vs1 && g[target(*ei, g)].state == vs2) ||
         (g[source(*ei, g)].state == vs2 && g[target(*ei, g)].state == vs1))) {
      count++;
    }
  }

  // since we looped over all edges of one type, no need to divide by 2
  return count;
}

//----------------------------------------------------------
/*! \brief Print graph statistics.
  
Prints the number of vertices in all possible states and, if desired, the
number of edges of all types connecting vertices of all possible states, as
well as the number of parallel edges connecting such vertices.

\param[in] g The graph containing the vertices and edges
\param[in] m The model to be used to find all possible states
\param[in] pairs Whether pairs should be counted as well (intensive on
computation) or not 
\ingroup graph_statistics
*/
template <typename Graph>
void print_graph_statistics(Graph& g, Model& m, bool pairs = true)
{
  std::cout << std::endl;
  std::cout << "Vertex count: " << std::endl;

  // count states
  for (unsigned int i=0; i < m.getVertexStates().size(); i++) {
    unsigned int count = count_vertices(g, i);
    if (count > 0) {
      std::cout << m.getVertexStates()[i] << ": " << count << std::endl;
    }
  }
   
  std::cout << std::endl;
  if (pairs) {
    
    std::cout << "Edge count: " << std::endl;
   
    // count state pairs
    for (unsigned int i=0; i < m.getEdgeTypes().size(); i++) {
      std::cout << m.getEdgeTypes()[i] << "-type: " << std::endl;
      for (unsigned int j=0; j < m.getVertexStates().size(); j++) {
        for (unsigned int k=j; k < m.getVertexStates().size(); k++) {
          unsigned int count = count_edges(g, j, k, i);
          if (count > 0) {
            std::cout << m.getVertexStates()[j] << m.getVertexStates()[k]
                      << ": " << count << std::endl;
          }
        }
      }
    }
   
    // count parallel pairs
    std::cout << "parallel:" << std::endl;
    for (unsigned int j=0; j < m.getVertexStates().size(); j++) {
      for (unsigned int k=j; k < m.getVertexStates().size(); k++) {
        unsigned int i = 0; // counting one type is enough
        unsigned int count = count_parallel_edges(g, j, k, i);
        if (count > 0) {
          std::cout << m.getVertexStates()[j] << m.getVertexStates()[k]
                    << ": " << count << std::endl;
        }
      }
    }
  }
  
  std::cout << std::endl;
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
  for (unsigned int i = 0; i < m.getVertexStates().size(); i++) {
    line << count_vertices(g, i) << '\t';
  }

  if (pairs) {
    // count pairs
    for (unsigned int i = 0; i < m.getEdgeTypes().size(); i++) {
      for (unsigned int j = 0; j < m.getVertexStates().size(); j++) {
        for (unsigned int k = j; k < m.getVertexStates().size(); k++) {
          line << count_edges(g, j, k, i) << '\t';
        }
      }
    }

    // count parallel pairs
    for (unsigned int j = 0; j < m.getVertexStates().size(); j++) {
      for (unsigned int k = j; k < m.getVertexStates().size(); k++) {
        unsigned int i = 0; // counting one type is enough
        line << count_parallel_edges(g, j, k, i) << '\t';
      }
    }
  }
  
  line << std::endl;

  ofile << line.str();
  return line.str();
}

//----------------------------------------------------------

#endif
