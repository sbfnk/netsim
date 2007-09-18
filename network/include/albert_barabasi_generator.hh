/*! \file albert_barabasi_generator.hh
  \brief Contains the iterator for Albert-Barabasi graphs.
*/
#ifndef BOOST_GRAPH_ALBERT_BARABASI_GENERATOR_HH
#define BOOST_GRAPH_ALBERT_BARABASI_GENERATOR_HH

#include <set>

#include <iterator>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graph_traits.hpp>

namespace boost {

  //----------------------------------------------------------
  /*! \brief Iterator for Albert-Barabasi graphs.
    
  Provides and interface to generate graphs according to Albert-Barabasi's
  preferential attachment rule.
  
  \ingroup graph_generators
  */
  template<typename RandomGenerator, typename Graph>
  class albert_barabasi_iterator
  {
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::vertices_size_type
    vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type
    edges_size_type;
    typedef typename boost::variate_generator
    <RandomGenerator&, boost::uniform_real<> > uniform_gen;

  public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef void difference_type;

    //! Constructor for past-the-end iterator.
    albert_barabasi_iterator() :
      gen(0), n(0), m(0), nodes(), targets(), source(0), target(0)
    {}

    /*! \brief Constructor.
      \param[in] gen Random generator to use.
      \param[in] n Number of vertices to create.
      \param[in] m Number of vertices to connect new nodes to at each step.
    */
    albert_barabasi_iterator(RandomGenerator& gen, vertices_size_type n, 
                          edges_size_type m = 1)
      : gen(new uniform_gen(gen, boost::uniform_real<> (0,1))),
        n(n), m(m), nodes(), targets(), source(m), target(0),
        current(source, target)
         
    {
      if (n<1) {
        // directly become past_the_end iterator
        source = n;
      }
      
      targets.insert(target);
    }

    reference operator*() const { return current; }
    pointer operator->() const { return &current; }

    albert_barabasi_iterator& operator++()
    {
      if (targets.size() == m) {
        // we have added m edges for this vertex, so we go to the next one
        for (std::set<unsigned int>::iterator it = targets.begin();
             it != targets.end(); it++) {
          nodes.push_back(*it);
          nodes.push_back(source);
        }
        targets.clear();
        source++;
      }

      if (source == m) {

        // we are still at first vertex, just add all links
        target++;
        targets.insert(target);

      } else {
        
        // choose random vertex
        double x;
        do {
          x = (*gen)();
          target = nodes[static_cast<unsigned int>(x*nodes.size())];
          // choose as target only if it is not a target yet
        } while (targets.insert(target).second == false);

      }

      current.first = source;
      current.second = target;
      return *this;
    }

    albert_barabasi_iterator operator++(int)
    {
      albert_barabasi_iterator temp(*this);
      ++(*this);
      return temp;
    }

    bool operator==(const albert_barabasi_iterator& rhs) const
    {
      if (!gen && rhs.gen) {
        return rhs == *this;
      } else if (gen && !rhs.gen) {
        return source == n;
      } else if (!gen && !rhs.gen) {
        return true;
      }
      return source == rhs.source && target == rhs.target;
    }
    
    bool operator!=(const albert_barabasi_iterator& rhs) const
    { return !(*this == rhs); }
    
  private:
        
    uniform_gen* gen;
    unsigned int n;
    edges_size_type m;
    std::vector<vertices_size_type> nodes;
    std::set<unsigned int> targets;
    vertices_size_type source;
    vertices_size_type target;
    value_type current;
  };

} // end namespace boost

#endif // BOOST_GRAPH_ERDOS_RENYI_GENERATOR_HPP
