/*! \file breadth_first_search.hh
  \brief Additional functions for boost::breadth_first_search to allow for
  searches along edges of only one type.
*/
#ifndef BREADTH_FIRST_SEARCH_HPP
#define BREADTH_FIRST_SEARCH_HPP

#include <boost/graph/breadth_first_search.hpp>

namespace boost {

  template <class IncidenceGraph, class Buffer, class BFSVisitor,
            class ColorMap, class EdgeType>
  void breadth_first_visit
    (const IncidenceGraph& g,
     typename graph_traits<IncidenceGraph>::vertex_descriptor s,
     EdgeType e,
     Buffer& Q, BFSVisitor vis, ColorMap color)
  {
    function_requires< IncidenceGraphConcept<IncidenceGraph> >();
    typedef graph_traits<IncidenceGraph> GTraits;
    typedef typename GTraits::vertex_descriptor Vertex;
    typedef typename GTraits::edge_descriptor Edge;
    function_requires< BFSVisitorConcept<BFSVisitor, IncidenceGraph> >();
    function_requires< ReadWritePropertyMapConcept<ColorMap, Vertex> >();
    typedef typename property_traits<ColorMap>::value_type ColorValue;
    typedef color_traits<ColorValue> Color;
    typename GTraits::out_edge_iterator ei, ei_end;

    put(color, s, Color::gray());             vis.discover_vertex(s, g);
    Q.push(s);
    while (! Q.empty()) {
      Vertex u = Q.top(); Q.pop();            vis.examine_vertex(u, g);
      for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
        if (g[*ei].type == e.type) { // check if edge has desired type
          Vertex v = target(*ei, g);            vis.examine_edge(*ei, g);
          ColorValue v_color = get(color, v);
          if (v_color == Color::white()) {      vis.tree_edge(*ei, g);
            put(color, v, Color::gray());       vis.discover_vertex(v, g);
            Q.push(v);
          } else {                              vis.non_tree_edge(*ei, g);
            if (v_color == Color::gray())       vis.gray_target(*ei, g);
            else                                vis.black_target(*ei, g);
          }
        } // end if
      } // end for
      put(color, u, Color::black());          vis.finish_vertex(u, g);
    } // end while
  } // breadth_first_visit


  template <class VertexListGraph, class Buffer, class BFSVisitor,
            class ColorMap, class EdgeType>
  void breadth_first_search
    (const VertexListGraph& g,
     typename graph_traits<VertexListGraph>::vertex_descriptor s,
     EdgeType e,
     Buffer& Q, BFSVisitor vis, ColorMap color)
  {
    // Initialization
    typedef typename property_traits<ColorMap>::value_type ColorValue;
    typedef color_traits<ColorValue> Color;
    typename boost::graph_traits<VertexListGraph>::vertex_iterator i, i_end;
    for (tie(i, i_end) = vertices(g); i != i_end; ++i) {
      vis.initialize_vertex(*i, g);
      put(color, *i, Color::white());
    }
    breadth_first_visit(g, s, e, Q, vis, color);
  }

  namespace detail {

    template <class VertexListGraph, class ColorMap, class BFSVisitor,
      class P, class T, class R, class EdgeType>
    void bfs_helper
      (VertexListGraph& g,
       typename graph_traits<VertexListGraph>::vertex_descriptor s,
       EdgeType e,
       ColorMap color,
       BFSVisitor vis,
       const bgl_named_params<P, T, R>& params)
    {
      typedef graph_traits<VertexListGraph> Traits;
      // Buffer default
      typedef typename Traits::vertex_descriptor Vertex;
      typedef boost::queue<Vertex> queue_t;
      queue_t Q;
      detail::wrap_ref<queue_t> Qref(Q);
      breadth_first_search
        (g, s, e,
         choose_param(get_param(params, buffer_param_t()), Qref).ref,
         vis, color);
    }

    //-------------------------------------------------------------------------
    // Choose between default color and color parameters. Using
    // function dispatching so that we don't require vertex index if
    // the color default is not being used.

    template <class ColorMap>
    struct bfs_dispatch_edgetype {
      template <class VertexListGraph, class P, class T, class R, class EdgeType>
      static void apply
      (VertexListGraph& g,
       typename graph_traits<VertexListGraph>::vertex_descriptor s,
       EdgeType e,
       const bgl_named_params<P, T, R>& params,
       ColorMap color)
      {
        bfs_helper
          (g, s, e, color,
           choose_param(get_param(params, graph_visitor),
                        make_bfs_visitor(null_visitor())),
           params);
      }
    };

    template <>
    struct bfs_dispatch_edgetype<detail::error_property_not_found> {
      template <class VertexListGraph, class P, class T, class R, class EdgeType>
      static void apply
      (VertexListGraph& g,
       typename graph_traits<VertexListGraph>::vertex_descriptor s,
       EdgeType e,
       const bgl_named_params<P, T, R>& params,
       detail::error_property_not_found)
      {
        std::vector<default_color_type> color_vec(num_vertices(g));
        default_color_type c = white_color;
        null_visitor null_vis;

        bfs_helper
          (g, s, e,
           make_iterator_property_map
           (color_vec.begin(),
            choose_const_pmap(get_param(params, vertex_index),
                              g, vertex_index), c),
           choose_param(get_param(params, graph_visitor),
                        make_bfs_visitor(null_vis)),
           params);
      }
    };

  } // namespace detail


  // Named Parameter Variant
  template <class VertexListGraph, class P, class T, class R, class EdgeType>
  void breadth_first_search
    (const VertexListGraph& g,
     typename graph_traits<VertexListGraph>::vertex_descriptor s,
     EdgeType e,
     const bgl_named_params<P, T, R>& params)
  {
    // The graph is passed by *const* reference so that graph adaptors
    // (temporaries) can be passed into this function. However, the
    // graph is not really const since we may write to property maps
    // of the graph.
    VertexListGraph& ng = const_cast<VertexListGraph&>(g);
    typedef typename property_value< bgl_named_params<P,T,R>,
      vertex_color_t>::type C;
    detail::bfs_dispatch_edgetype<C>::apply(ng, s, e, params,
                                            get_param(params, vertex_color));
  }


  // This version does not initialize colors, user has to.

  template <class IncidenceGraph, class P, class T, class R, class EdgeType>
  void breadth_first_visit
    (const IncidenceGraph& g,
     typename graph_traits<IncidenceGraph>::vertex_descriptor s,
     EdgeType e,
     const bgl_named_params<P, T, R>& params)
  {
    // The graph is passed by *const* reference so that graph adaptors
    // (temporaries) can be passed into this function. However, the
    // graph is not really const since we may write to property maps
    // of the graph.
    IncidenceGraph& ng = const_cast<IncidenceGraph&>(g);

    typedef graph_traits<IncidenceGraph> Traits;
    // Buffer default
    typedef typename Traits::vertex_descriptor vertex_descriptor;
    typedef boost::queue<vertex_descriptor> queue_t;
    queue_t Q;
    detail::wrap_ref<queue_t> Qref(Q);

    breadth_first_visit
      (ng, s, e,
       choose_param(get_param(params, buffer_param_t()), Qref).ref,
       choose_param(get_param(params, graph_visitor),
                    make_bfs_visitor(null_visitor())),
       choose_pmap(get_param(params, vertex_color), ng, vertex_color)
       );
  }

} // namespace boost

#endif // BOOST_GRAPH_BREADTH_FIRST_SEARCH_HPP

