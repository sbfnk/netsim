/*! \file tree_generator.hh
  \brief Contains an iterator for tree-like networks.
*/
#ifndef TREE_GENERATOR_HH
#define TREE_GENERATOR_HH

#include <boost/graph/graph_traits.hpp>
#include <math.h>
#include <vector>

namespace boost {

  //----------------------------------------------------------
  /*! \brief Iterator for a tree.
  \ingroup graph_generators
  */
  template <typename Graph>
  class tree_iterator
  {
    typedef typename graph_traits<Graph>::vertices_size_type
    vertices_size_type;
  public:
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;

    //! Constructor for past-the-end iterator.
    tree_iterator(): past_the_end(true) {}
    /*! \brief Constructor.
      \param[in] n The number of vertices to connect.
      \param[in] b The number of branches to create from each vertex.
    */
    tree_iterator(vertices_size_type n, unsigned int b)
      : past_the_end(false), vertices(n), branches(b), depth(0), depth_last(0),
        source(0), target(1), current(0,1), link_count(0)
    {}
         
    reference operator*() const { return current; }
    pointer operator->() const { return &current; }

    tree_iterator& operator++()
    {
      ++link_count;
      if (link_count == branches) {
        if (source == depth_last) {
          ++depth;
          depth_last += static_cast<vertices_size_type>(pow(branches, depth));
          target = depth_last;
        }
        ++source;
        link_count = 0;
      }
      ++target;
      current.first = source;
      current.second = target;
      return *this;
    }
         
    tree_iterator operator++(int)
    {
      tree_iterator temp(*this);
      ++(*this);
      return temp;
    }

    bool operator==(const tree_iterator& rhs) const
    {
      if (past_the_end && !rhs.past_the_end) {
        return rhs == *this;
      } else if (!past_the_end && rhs.past_the_end) {
        return (target >= vertices);
      } else if (past_the_end && rhs.past_the_end) {
        return true;
      }
      return source == rhs.source && target == rhs.target;
    }
         
    bool operator!=(const tree_iterator& rhs) const
    { return !(*this == rhs); }

    bool past_the_end;
    vertices_size_type vertices;
    unsigned int branches;
    unsigned int depth;
    vertices_size_type depth_last;
    vertices_size_type source;
    vertices_size_type target;
    value_type current;
    unsigned int link_count;
  };

}
#endif
