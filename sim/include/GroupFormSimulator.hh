/*! \file GroupFormSimulator.hh
  \brief The Simulators::GroupFormSimulator class.
*/

#ifndef GROUPFORMSIMULATOR_HH
#define GROUPFORMSIMULATOR_HH

#include "Tree.hh"
#include "Simulator.hh"

#include <boost/graph/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//! \addtogroup gillespie_simulator Gillespie simulator

namespace Simulators {
  
  /*! \brief The Group formation simulation
    
  \ingroup gillespie_simulator
  */
  template <typename RandomGenerator, typename Graph>
  class GroupFormSimulator :
    virtual public Simulator
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
    edge_descriptor;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
    typedef typename boost::graph_traits<Graph>::adjacency_iterator
    adjacency_iterator;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type
    vertex_index_type;
    typedef typename boost::vertex_property_type<Graph>::type::value_type
    vertex_property_type;
    typedef typename boost::variate_generator
    <RandomGenerator&, boost::uniform_real<> > uniform_gen;

  public:

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    GroupFormSimulator(RandomGenerator& r, Graph& g, const Model& m, 
                       unsigned int numGroups = 1, double alpha=0.,
                       unsigned int v = 0) :
      Simulator(m, v), randGen(r, boost::uniform_real<> (0,1)), graph(g),
      groups(numGroups), alpha(alpha), active(), verbose(v) {;}
    
  
    ~GroupFormSimulator() {;}
      
    void initialise();
    bool updateState();
  
   void print();

  private:
  
    uniform_gen randGen; //!< The random generator to be used for choosing events
    Graph& graph; //!< The graph determining how vertices can effect another

    unsigned int groups; //!< Number of seed groups
    double alpha;
    std::vector<vertex_descriptor> active; //!< Active nodes

    unsigned int verbose;
  
  };

  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
  
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  void GroupFormSimulator<RandomGenerator, Graph>::initialise()
  {
    Simulator::initialise();

    // set trait for each vertex 
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      graph[*vi].state = State(0, (randGen)());
      if (verbose >=1) {
        std::cout << "Associating trait " << graph[*vi].state.detail
                  << " with vertex " << *vi << std::endl;
      }
    }

    // seed groups
    for (unsigned int i = 1; i <= groups; ++i) {
      vertex_descriptor v = random_vertex(graph, randGen);
      if (graph[v].state.base == 0) {
        graph[v].state.base = i;
        active.push_back(v);
        if (verbose >=1) {
          std::cout << "Seeding group " << i << " with vertex " << v << std::endl;
        }
      } else {
        --i;
      }
    }
    
  }

  //----------------------------------------------------------
  /*! \brief Perform an update and process one event.
  
  \return true if an event is processed, false if no events can happen or
  something goes wrong 
  */
  template <typename RandomGenerator, typename Graph>
  bool GroupFormSimulator<RandomGenerator, Graph>::updateState()
  {
    
    if (active.size() > 0) {

      if (verbose >=2) {
        std::cout << "Active nodes: " << std::endl;
        for (unsigned int i = 0; i < active.size(); ++i) {
          std::cout << active[i] << std::endl;
        }
      }

      // group formation stage
      if (verbose >=1) {
        std::cout << "Group formation" << std::endl;
      }

      // choose an active node at random
      unsigned int randActive =
        static_cast<unsigned int>((randGen)() * active.size());
      if (verbose >=2) {
        std::cout << "Active node " << active[randActive]
                  << " is sending out invitations" << std::endl;
      }
    
      // invite all neighbours to group
      out_edge_iterator oi, oi_end;;
   
      for (tie(oi, oi_end) = boost::out_edges(active[randActive], graph);
           oi != oi_end; ++oi) {
        if (graph[target(*oi, graph)].state.base == 0) {
          double randAccept = randGen();
          if (verbose >=2) {
            std::cout << "Invitation to " << target(*oi, graph) << " ("
                      << graph[active[randActive]].state.detail << ","
                      << graph[target(*oi, graph)].state.detail << ")";
          }
          if (randAccept < alpha*(1-fabs(graph[active[randActive]].state.detail -
                                         graph[target(*oi, graph)].state.detail))) {
            // accept invitation
            graph[target(*oi, graph)].state.base =
              graph[active[randActive]].state.base;
            // add to active nodes
            active.push_back(target(*oi, graph));
            if (verbose >=2) {
              std::cout << " accepted." << std::endl;
            }
          } else {
            if (verbose >=2) {
              std::cout << " refused." << std::endl;
            }
          }
        }
      }

      // remove randActive from active nodes
      active.erase(active.begin() + randActive);

      // update time if group formation phase is finished
      if (active.size() == 0) updateTime(1.);
      
    } else {
      // network update stage
      if (verbose >=1) {
        std::cout << "Network update stage, time=" << getTime() << std::endl;
      }

      // go through members group by group
      std::vector<vertex_descriptor> group;
      
      for (unsigned int i = 1; i <= groups; ++i) {
        if (verbose >=2) {
          std::cout << "Group " << i << std::endl;
        }
        group.clear();
        
        vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
          if (graph[*vi].state.base == i) group.push_back(*vi);
        }

        // loop over all member of group and update ties
        for (unsigned int i = 0; i < group.size(); ++i) {
          if (verbose >=2) {
            std::cout << "Group member " << i << ": " << group[i] << std::endl;
          }
          
          out_edge_iterator oi, oi_end;;

          // find least similar neighbour
          double biggest_distance = 0;
          edge_descriptor weakest_link;
          
          for (tie(oi, oi_end) = boost::out_edges(group[i], graph);
               oi != oi_end; ++oi) {
            double distance = fabs(graph[group[i]].state.detail -
                                   graph[target(*oi, graph)].state.detail);
            
            if (distance > biggest_distance) {
              weakest_link = *oi;
              biggest_distance = distance;
            }
          }
          // find most similar neighbour in group
          double lowest_distance = 1.;
          vertex_descriptor similar_neighbour;

          for (unsigned int j = 0; j < group.size(); ++j) {
            double distance = 1.;
            if (i != j) {
              distance = fabs(graph[group[i]].state.detail -
                              graph[group[j]].state.detail);
            }
            if (distance < lowest_distance &&
                edge(group[i], group[j], graph).second == false) {
              similar_neighbour = group[j];
              lowest_distance = distance;
            }
          }

          if (lowest_distance < 1 && biggest_distance > 0) {
            if (verbose >=2) {
              std::cout << "removing link to " << target(weakest_link, graph)
                        << std::endl;
            }
            // remove weakest link
            boost::remove_edge(weakest_link, graph);
            // link to group node
            if (verbose >=2) {
              std::cout << "adding link to " << similar_neighbour
                        << std::endl;
            }
            boost::add_edge(group[i], similar_neighbour, graph);
          }
        }
      }
      // seed groups again
      vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        graph[*vi].state.base = 0;
      }
      for (unsigned int i = 1; i <= groups; ++i) {
        vertex_descriptor v = random_vertex(graph, randGen);
        if (graph[v].state.base == 0) {
          graph[v].state.base = i;
          active.push_back(v);
          if (verbose >=2) {
            std::cout << "Seeding group " << i << " with vertex " << v << std::endl;
          }
        } else {
          --i;
        }
      }
    }

    return true;
  }

  //----------------------------------------------------------
  /*! \brief Print the state of the simulation. 
  
  Prints the state of all vertices, as well as the events that can happen to
  them and at which rate.
  */
  template <typename RandomGenerator, typename Graph>
  void GroupFormSimulator<RandomGenerator, Graph>::print()
  {
    std::cout << std::endl;
  }
  
}

#endif
