#ifndef BOOST_GRAPH_ALBERT_BARABASI_GENERATOR_HH
#define BOOST_GRAPH_ALBERT_BARABASI_GENERATOR_HH

#include <set>

#include <iterator>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/graph_traits.hpp>

namespace boost {

  template<typename RandomGenerator, typename Graph>
  class albert_barabasi_iterator
  {
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::vertices_size_type
    vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type
    edges_size_type;

  public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef void difference_type;

    albert_barabasi_iterator() :
      gen(0), n(0), m(0), nodes(), targets(), source(0), target(0)
    {}
    albert_barabasi_iterator(RandomGenerator& gen, vertices_size_type n, 
                          edges_size_type m = 1)
      : gen(&gen), n(n), m(m), nodes(), targets(),
        source(m), target(0), current(source, target)
         
    {
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
          uniform_01<RandomGenerator, double> rand01(*gen);
          x = rand01();
          target = nodes[static_cast<unsigned int>(x*nodes.size())];
          *gen = rand01.base();
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
        return source == n;      } else if (!gen && !rhs.gen) {
        return true;
      }
      return source == rhs.source && target == rhs.target;
    }
    
    bool operator!=(const albert_barabasi_iterator& rhs) const
    { return !(*this == rhs); }
    
  private:
        
    RandomGenerator* gen;
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
