/*! \file community_generator.hh
  \brief Contains the iterator for weighted community graphs.
*/
#ifndef BOOST_GRAPH_COMMUNITY_GENERATOR_HH
#define BOOST_GRAPH_COMMUNITY_GENERATOR_HH

#include <set>

#include <iterator>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graph_traits.hpp>

namespace boost {

  //----------------------------------------------------------
  /*! \brief Iterator for community graphs.
    
  Provides and interface to generate graphs according to Kumpula et al
  (arXiv: 0708.0925v1)
  
  \ingroup graph_generators
  */
  template<typename RandomGenerator, typename Graph>
  class community_iterator
  {
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::vertices_size_type
    vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type
    edges_size_type;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef typename graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename boost::variate_generator
    <RandomGenerator&, boost::uniform_real<> > uniform_gen;

  public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef void difference_type;

    //! Constructor for past-the-end iterator.
    community_iterator() :
      gen(0), graph(0)
    {}

    /*! \brief Constructor.
      \param[in] gen Random generator to use.
      \param[in] n Number of vertices to create.
      \param[in] m Number of vertices to connect new nodes to at each step.
    */
    community_iterator(RandomGenerator& r, Graph* const g, 
                       double delta = 0.5, double pl = 1e-3,
                       double pr = 5e-4, double pd = 1e-3, 
                       unsigned int iterations = 25000)
      : gen(new uniform_gen(r, boost::uniform_real<> (0,1))),
        graph(g), n(num_vertices(*graph)), delta(delta), pl(pl), pr(pr), pd(pd),
        itCount(iterations)
         
    {
      if (n<2) {
        // directly become past_the_end iterator
        itCount = 0;
      } else {
        // create all edges to insert in first step -- all random
        double x;
        unsigned int s, t;
        for (s = 0; s < n; ++s) {
//         std::cout << "Choosing random link for " << s << std::endl;
          do {
            x = (*gen)();
            t = static_cast<unsigned int>(x*n);
            // choose as t only if it is not a t yet
          } while ((s == t) || (edge(s, t, *graph).second));
          edges.push_back(std::make_pair(s, t));
//           std::cout << s << "--" << t << " (random)" << std::endl;
        }
        current.first = edges[0].first;
        current.second = edges[0].second;
        edges.erase(edges.begin());
      }
    }

    reference operator*() const { return current; }
    pointer operator->() const { return &current; }

    community_iterator& operator++()
    {
      if (edges.size() == 0) {
        // we have added all edges, so we go through a deletion cycle
        double x;
        vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(*graph); vi != vi_end; vi++) {
          x = (*gen)();
          if (x < pd) {
            clear_vertex(*vi, *graph);
          }
        }

        --itCount;
        std::cout << itCount << " cycles left" << std::endl;
        if (itCount > 0) {
          // recreate edge list
          unsigned int s(0), t(0);
          for (s = 0; s < n; ++s) {
            x = (*gen)();
            if (x < pr || out_degree(s, *graph) == 0) {
              // choose new t at random
              do {
                x = (*gen)();
                t = static_cast<unsigned int>(x*n);
              } while ((s == t) || (edge(s, t, *graph).second));
              edges.push_back(std::make_pair(s, t));
            } else {
              // perform first local search
              out_edge_iterator oi, oi_end;
              double strength = 0;
              for (tie(oi, oi_end) = out_edges(s, *graph);
                   oi != oi_end; oi++) {
                strength += (*graph)[*oi].weight;
              }
              x = (*gen)();
              x *= strength;
              tie(oi, oi_end) = out_edges(s, *graph);
              while (x > (*graph)[*oi].weight) {
                x -= (*graph)[*oi].weight;
                oi++;
              }
              edge_descriptor first_edge = *oi;
              vertex_descriptor temp_vertex = target(first_edge, *graph);
              // perform second local search
              if (out_degree(temp_vertex, *graph) > 1) {
                strength = 0;
                for (tie(oi, oi_end) = out_edges(temp_vertex, *graph);
                     oi != oi_end; oi++) {
                  if (target(*oi, *graph) != s) {
                    strength += (*graph)[*oi].weight;
                  }
                }
                x = (*gen)();
                x *= strength;
                tie(oi, oi_end) = out_edges(temp_vertex, *graph);
                while (x > (*graph)[*oi].weight) {
                  if (target(*oi, *graph) != s) x -= (*graph)[*oi].weight;
                  oi++;
                }
                if (target(*oi, *graph) == s) oi++;
		edge_descriptor second_edge = *oi;
                std::pair<edge_descriptor, bool> chosen_edge =
                  edge(s, target(second_edge, *graph), *graph);
                if (chosen_edge.second) {
                  // edge exists, so weight will be increased
//                   edges.push_back(std::make_pair(s, target(second_edge, *graph)));
                  (*graph)[chosen_edge.first].weight += delta;
                } else {
                  x = (*gen)();
                  if (x < pl) {
                    edges.push_back(std::make_pair(s, target(second_edge, *graph)));
                  }
                }
                (*graph)[first_edge].weight += delta;
                (*graph)[second_edge].weight += delta;
//                   edges.push_back(std::make_pair(source(first_edge, *graph),
//                                                  target(first_edge, *graph)));
//                   edges.push_back(std::make_pair(source(second_edge, *graph),
//                                                  target(second_edge, *graph)));
              }
            }
          }
          if (edges.size() > 0) {
            current.first = edges[0].first;
            current.second = edges[0].second;
            edges.erase(edges.begin());
          } else {
            ++(*this);
          }
        } else {
	  return *this;
	}
      } else {
        current.first = edges[0].first;
        current.second = edges[0].second;
        edges.erase(edges.begin());
        // check if we have already added this one
        std::pair<edge_descriptor, bool> newEdge =
          edge(current.first, current.second, *graph);
        if (newEdge.second) {
          (*graph)[newEdge.first].weight += delta;
          ++(*this);
        }

      }

      return *this;
    }

    community_iterator operator++(int)
    {
      community_iterator temp(*this);
      ++(*this);
      return temp;
    }

    bool operator==(const community_iterator& rhs) const
    {
      if (!gen && rhs.gen) {
        return rhs == *this;
      } else if (gen && !rhs.gen) {
        return (itCount == 0);
      } else if (!gen && !rhs.gen) {
        return true;
      }
      return current.first == rhs.current.first &&
        current.second == rhs.current.second;
    }
    
    bool operator!=(const community_iterator& rhs) const
    { return !(*this == rhs); }
    
  private:
        
    uniform_gen* gen;
    Graph* graph;
    unsigned int n;
    double delta, pl, pr, pd;
    unsigned int itCount;
    std::vector<std::pair<vertices_size_type, vertices_size_type> > edges;
    value_type current;
  };

} // end namespace boost

#endif // BOOST_GRAPH_ERDOS_RENYI_GENERATOR_HPP
