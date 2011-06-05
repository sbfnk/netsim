/*! \file NetEvoSimulator.hh
  \brief The Simulators::NetEvoSimulator class.
*/

#ifndef NETEVOSIMULATOR_HH
#define NETEVOSIMULATOR_HH

#include <math.h>
#include <errno.h>
#include <algorithm>

#include "GillespieSimulator.hh"
#include "Vertex.hh"
#include "EvoModel.hh"

#include "network/include/community_structure.hh"

#include <boost/graph/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//! \addtogroup simulator Simulator

namespace Simulators {

  /*! \brief The Group formation simulation
    
  \ingroup simulator
  */
  template <typename RandomGenerator, typename Graph>
  class NetEvoSimulator :
    virtual public GillespieSimulator<RandomGenerator, Graph>
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
    edge_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
    typedef typename boost::graph_traits<Graph>::adjacency_iterator
    adjacency_iterator;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type
    vertex_index_type;
    typedef typename boost::vertex_property_type<Graph>::type::value_type
    vertex_property_type;

  public:

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    NetEvoSimulator(RandomGenerator& r, Graph& g, unsigned int v = 0);
    ~NetEvoSimulator() {;}
      
    bool initialise();
    bool updateState();
  
    void print();
    
  };

  template <typename RandomGenerator, typename Graph>
  NetEvoSimulator<RandomGenerator, Graph>::
  NetEvoSimulator(RandomGenerator& r, Graph& g, unsigned int v) :
    GillespieSimulator<RandomGenerator, Graph>(r, g, v)
  {
    this->knownModels.push_back
      (std::make_pair("Evo",
                      new Models::EvoModel<Graph>
                      (this->getVerbose())));
  }
  

  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
    
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  bool NetEvoSimulator<RandomGenerator, Graph>::initialise()
  {
    bool result = true;

    result &= GillespieSimulator<RandomGenerator, Graph>::initialise();

    Graph& graph = this->getGraph();
    
    vertex_iterator vi, vi_end, vi2, vi_end2;
    edge_iterator ei, ei_end;

    // give everyone a random trait
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      EvoModelState* state = dynamic_cast<EvoModelState*>(graph[*vi].state);
      state->setTrait(this->getRandom());
    }

    // connect everyone on the base layer
    for (tie(vi, vi_end) = vertices(graph); (vi + 1) != vi_end; ++vi) {
      for (vi2 = (vi + 1); vi2 != vi_end;  ++vi2) {
        add_edge(*vi, *vi2, Edge(Models::EvoModel<Graph>::Base), graph);
      }
    }

    return result;
  }
  
  //----------------------------------------------------------
  /*! \brief Perform an update and process one event.
  
  \return true if an event is processed, false if no events can happen or
  something goes wrong 
  */
  template <typename RandomGenerator, typename Graph>
  bool NetEvoSimulator<RandomGenerator, Graph>::updateState()
  {
    bool res = GillespieSimulator<RandomGenerator, Graph>::updateState();
    
    Graph& graph = this->getGraph();
    
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      EvoModelState* state = dynamic_cast<EvoModelState*>(graph[*vi].state);
      if (state->getPending() == true) {
        add_edge(*vi, state->getNewNb(),
                 Edge(Models::EvoModel<Graph>::Interaction), graph);
        state->clearPending();
      }
    }

    return res;
  }

  //----------------------------------------------------------
  /*! \brief Print the state of the simulation. 
  
  Prints the state of all vertices, as well as the events that can happen to
  them and at which rate.
  */
  template <typename RandomGenerator, typename Graph>
  void NetEvoSimulator<RandomGenerator, Graph>::print()
  {
    std::cout << std::endl;
  }
  
}

#endif
