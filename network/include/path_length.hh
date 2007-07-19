/*! \file path_length.hh
  \brief Functions for network properties based on path length.
*/

#ifndef PATH_LENGTH_HH
#define PATH_LENGTH_HH

#include "breadth_first_search.hpp"
#include "graph_statistics.hh"

//! \addtogroup graph_statistics Graph statistics
//! \addtogroup helper_functions Helper functions

//----------------------------------------------------------
namespace boost {

  /*! \brief Calculate the average shortest path length.

  Calculates the average shortest path length for edges of a given type by
  looping over all vertices and then using a breadth_first_search to record the
  distances to all other edges.
  
  \param[in] g The graph to examine.
  \param[in] et The edge type to consider.
  \ingroup graph_statistics
  */

  template <class Graph, class EdgeType>
  double avg_shortest_path_length(Graph& g, EdgeType et)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::vertices_size_type
      vertices_size_type;
    typedef typename graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    
    vertices_size_type d[num_vertices(g)];

    double distance_sum = 0.;
    vertex_iterator vi, vi_end;

    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      std::fill_n(d, num_vertices(g), 0);

      breadth_first_search(g, *vi, et,
        visitor(make_bfs_visitor(record_distances(d, on_tree_edge()))));

      for (unsigned int i = 0; i < num_vertices(g); i++) distance_sum += d[i];
    }

    return (distance_sum/(num_vertices(g)*num_vertices(g)));
  }

  /*! \brief Visitor for getting the distance to a number of vertices

  \param[in] et The edge type to consider.
  \ingroup graph_statistics
  */
  template <class VertexDescriptor, class DistanceMap, class Tag>
  struct vertices_distance_recorder
    : public base_visitor<distance_recorder<DistanceMap, Tag> >
  {
    typedef Tag event_filter;
    typedef typename std::vector<VertexDescriptor> vertex_vector;
    
    vertices_distance_recorder(DistanceMap pa_all, DistanceMap pa_found,
                      vertex_vector& vec)
      : search_list(vec), m_distance(pa_all), m_distance_found(pa_found),
        found(0) {}
    template <class Edge, class Graph>
    void operator()(Edge e, const Graph& g) {
      typename graph_traits<Graph>::vertex_descriptor
        u = source(e, g), v = target(e, g);
      put(m_distance, v, get(m_distance, u) + 1);
      
      if (std::find(search_list.begin(), search_list.end(), v) !=
          search_list.end()) {
        put(m_distance_found, v, get(m_distance, v));
        ++found;
        if (found == search_list.size()) throw 14;
      }
    }
    vertex_vector& search_list;
    DistanceMap m_distance;
    DistanceMap m_distance_found;
    unsigned int found;
  };
  template <class VertexDescriptor, class DistanceMap, class Tag>
  vertices_distance_recorder<VertexDescriptor, DistanceMap, Tag>
  record_vertex_distances(DistanceMap pa_all, DistanceMap pa_found,
                          typename std::vector<VertexDescriptor> vec, Tag) {
    return vertices_distance_recorder<VertexDescriptor, DistanceMap, Tag>
      (pa_all, pa_found, vec);
  }

  /*! \brief Calculate the average shortest path length.

  Calculates the average shortest path length for edges of a given type by
  looping over all vertices and then using a breadth_first_search to record the
  distances to all other edges.
  
  \param[in] g The graph to examine.
  \param[in] et The edge type to consider.
  \ingroup graph_statistics
  */

  template <class Graph, class EdgeType>
  double avg_nb_shortest_path_length(Graph& g, EdgeType et1, EdgeType et2)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::vertices_size_type
      vertices_size_type;
    typedef typename graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    if (et1.type == et2.type) {
      return 1.;
    }
    
    vertices_size_type d[num_vertices(g)];
    vertices_size_type d_nb[num_vertices(g)];

    double distance_sum = 0.;
    unsigned int distance_count = 0;

    vertex_iterator vi, vi_end;

    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      std::fill_n(d, num_vertices(g), 0);
      std::fill_n(d_nb, num_vertices(g), 0);

      std::vector<vertex_descriptor> et1_neighbours;
      typename graph_traits<Graph>::out_edge_iterator oi, oi_end;
      for (tie(oi, oi_end) = out_edges(*vi, g); oi != oi_end; oi++) {
        if (g[*oi].type == et1.type) et1_neighbours.push_back(target(*oi, g));
      }
      
      if (et1_neighbours.size() > 0) {
        try {
          breadth_first_search(g, *vi, et2,
            visitor(make_bfs_visitor(record_vertex_distances(d, d_nb,
                                                             et1_neighbours,
                                                             on_tree_edge()))));
        } catch (...) {}

        unsigned int nb_count = 0;
        for (unsigned int i = 0; nb_count < et1_neighbours.size(); ++i) {
          if (d_nb[i] > 0) {
            ++nb_count;
            distance_sum += d_nb[i];
          }
        }
        distance_count += nb_count;
      }
    }

    if (distance_count > 0) return (distance_sum/distance_count);
    else return 0.;
  }

} // namespace boost

//----------------------------------------------------------
  

#endif
