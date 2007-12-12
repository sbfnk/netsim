/*! \file path_length.hh
  \brief Functions for network properties based on path length.
*/

#ifndef PATH_LENGTH_HH
#define PATH_LENGTH_HH

#include <boost/multi_array.hpp>

#include "breadth_first_search.hpp"

//! \addtogroup graph_statistics Graph statistics
//! \addtogroup helper_functions Helper functions

//----------------------------------------------------------
namespace boost {

  /*! \brief Visitor for getting the distance to the nearest infected

  This differs from the standard boost::distance_recorder by searching only for
  the nearest infected. Once the an infected has been found, an
  exception is thrown which can be caught in the main routine.
  
  \ingroup graph_statistics
  */
  template <class Model, class DistanceMap, class Tag>
  struct infected_distance_recorder
    : public base_visitor<distance_recorder<DistanceMap, Tag> >
  {
    typedef Tag event_filter;
    
    infected_distance_recorder(DistanceMap pa_all, const Model& model)
      : m_distance(pa_all), m(model)
    {}
    template <class Edge, class Graph>
    void operator()(Edge e, const Graph& g) {
      typename graph_traits<Graph>::vertex_descriptor
        u = source(e, g), v = target(e, g);
      put(m_distance, v, get(m_distance, u) + 1);
      
      if (m.isInfected(g[v].state)) throw get(m_distance, v);
    }
    DistanceMap m_distance;
    const Model& m;
  };
  template <class Model, class DistanceMap, class Tag>
  infected_distance_recorder<Model, DistanceMap, Tag>
  record_infected_distance(DistanceMap pa_all, const Model& m, Tag) {
    return infected_distance_recorder<Model, DistanceMap, Tag>
      (pa_all, m);
  }

  /*! \brief Calculate the distance to the nearest infected for all vertices in
    a graph..

  Calculates the average path length to the nearest infected over a given edge type
  
  \param[in] g The graph to examine.
  \param[in] m The model providing information about infected
  \param[in] et The edge type to pick the neighbours from.
  */

  template <class Graph, class EdgeType, class Model>
  std::vector<unsigned int>
  nearest_infected(Graph& g, EdgeType et, const Model& m)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::vertices_size_type
      vertices_size_type;
    typedef typename graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::vector<unsigned int> infected_distances;

    vertices_size_type d[num_vertices(g)];

    vertex_iterator vi, vi_end;

    // loop over all vertices
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {

      std::fill_n(d, num_vertices(g), 0);

      unsigned int distance = 0;
      
      // breadth_first_search throws an exception as soon as all the shortest
      // path lengths for all entries in et1_neighbours have been found
      try {
        breadth_first_search(g, *vi, et,
            visitor(make_bfs_visitor(record_infected_distance(d, m,
                                                              on_tree_edge()))));
      } catch (unsigned int e) {
        distance = e;
      }
      
      infected_distances.push_back(distance);
      
    }

    return infected_distances;
    
  }

} // namespace boost

//----------------------------------------------------------
  

#endif
