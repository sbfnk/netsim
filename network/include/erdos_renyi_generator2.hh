#ifndef BOOST_GRAPH_ERDOS_RENYI_GENERATOR2_HH
#define BOOST_GRAPH_ERDOS_RENYI_GENERATOR2_HH

#include <iostream>

#include <iterator>
#include <utility>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost {

  template<typename RandomGenerator, typename Graph>
  class erdos_renyi_iterator2
  {
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::vertices_size_type
    vertices_size_type;

    BOOST_STATIC_CONSTANT
    (bool,
     is_undirected = (is_base_and_derived<undirected_tag,
                      directed_category>::value
                      || is_same<undirected_tag, directed_category>::value));

  public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef void difference_type;

    erdos_renyi_iterator2() : gen(0), n(0), allow_self_loops(false), prob(0.0),
                              source(0), target(0)
    {}
    erdos_renyi_iterator2(RandomGenerator& gen, vertices_size_type n, 
                          double prob = 0.0, bool allow_self_loops = false)
       : gen(&gen), n(n), allow_self_loops(allow_self_loops), prob(prob),
         source(0), target(allow_self_loops? 0 : 1),
         current(0, allow_self_loops? 0 : 1)
    { }

    reference operator*() const { return current; }
    pointer operator->() const { return &current; }
    
    erdos_renyi_iterator2& operator++()
    { 
      next();
      return *this;
    }

    erdos_renyi_iterator2 operator++(int)
    {
      erdos_renyi_iterator2 temp(*this);
      ++(*this);
      return temp;
    }

    bool operator==(const erdos_renyi_iterator2& rhs) const
    {
      if (!gen && rhs.gen) {
        return rhs == *this;
      } else if (gen && !rhs.gen) {
        return source == n - 1;
      } else if (!gen && !rhs.gen) {
        return true;
      }
      return source == rhs.source && target == rhs.target;
    }
    
    bool operator!=(const erdos_renyi_iterator2& rhs) const
    { return !(*this == rhs); }

  private:
    void next()
    {
      ++target;
      if (target == n) {
        ++source;
        if (is_undirected) {
          target = source+1;
        } else {
          target = 0;
        }
      }
      
      uniform_01<RandomGenerator, double> rand01(*gen);
      double x = rand01();
      *gen = rand01.base();
      
      if (x > prob) {
        ++(*this);
      }

      current.first = source;
      current.second = target;
    }        
        
    RandomGenerator* gen;
    vertices_size_type n;
    bool allow_self_loops;
    double prob;
    vertices_size_type source;
    vertices_size_type target;
    value_type current;
  };

} // end namespace boost

#endif // BOOST_GRAPH_ERDOS_RENYI_GENERATOR_HPP
