/*! \file GillespieSimulator.hh
  \brief The Simulators::GillespieSimulator class.
*/

#ifndef GILLESPIESIMULATOR_HH
#define GILLESPIESIMULATOR_HH

#include "Tree.hh"
#include "Simulator.hh"

//! \addtogroup gillespie_simulator Gillespie simulator

namespace Simulators {
  
  /*! \brief The Gillespie simulation
    
  This contains the Gillespie simulation, able to run a simulation on a graph
  using a given model, with a sequential update rule according to
  Gillespie (1977). The rates of the events are stored in a Tree for quick
  access, which is generated using the member function initialise. The actual
  updating of the graph is done in updateState.
  
  \ingroup gillespie_simulator
  */
  template <typename RandomGenerator, typename Graph>
  class GillespieSimulator :
    virtual public Simulator
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;
    typedef typename boost::graph_traits<Graph>::adjacency_iterator
    adjacency_iterator;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type
    vertex_index_type;
    typedef typename boost::vertex_property_type<Graph>::type::value_type
    vertex_property_type;
    typedef typename boost::uniform_01<RandomGenerator, double> uniform_gen;

  public:

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    GillespieSimulator(RandomGenerator& r, Graph& g, const Model& m,
                       unsigned int v = 0) :
      Simulator(m, v), randGen(r), graph(g) {;}
    ~GillespieSimulator() {;}
      
    void initialise();
    bool updateState();
  
    void print();

  private:
  
    uniform_gen randGen; //!< The random generator to be used for choosing events
    Graph& graph; //!< The graph determining how vertices can effect another
    Tree::Tree<unsigned int> tree; //!< The tree holding the event rates
  
  };

  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
  
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  void GillespieSimulator<RandomGenerator, Graph>::initialise()
  {
    // get simulation variables
    Graph& graph = this->getGraph();
    Model& model = this->getModel();
    unsigned int verbose = this->getVerbose();

    // generate event list for each vertex
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      generateEventList(graph, *vi, model, verbose);
    }

    // generate the tree holding rate sums corresponding to vertices
    generateTree(tree,graph,
                 get(&vertex_property_type::rateSum,graph),
                 get(boost::vertex_index,graph));
  }

  //----------------------------------------------------------
  /*! \brief Perform an update and process one event.
  
  Updates the state of the graph by advancing one time step
  and choosing a random event to process. After that, the list of events is
  updated for affected vertices.

  \return true if an event is processed, false if no events can happen or
  something goes wrong 
  */
  template <typename RandomGenerator, typename Graph>
  bool GillespieSimulator<RandomGenerator, Graph>::updateState()
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;

    // get simulation variables
    const Model& model = this->getModel();
    unsigned int verbose = this->getVerbose();

    // exit if nothing can happen -- exit also if the rate sum is very small
    // to prevent rounding problems
    if (tree.getTopBin()->getRateSum() < 1e-8) {
      return false;
    }
  
    if (verbose >= 2) {
      std::cout << "choose event, total sum of rates is "
                << tree.getTopBin()->getRateSum() << std::endl;
      print();
    }
         
    // draw a random number from [0,1) for the timestep advance
    double randNo = (randGen)();
    updateTime(-log(randNo)/tree.getTopBin()->getRateSum());
         
    // draw another random number from [0,1) for picking the event
    randNo = (randGen)();
    unsigned int* eventVertex = tree.pickRandomElement(randNo);
    if (eventVertex) {
      // process vertex event
      double tempSum = .0;
      vertex_descriptor v = vertex(*eventVertex, graph);

      // find event corresponding to randNo
      std::vector<event>::iterator it = graph[v].events.begin();
      while (it != graph[v].events.end() && tempSum < randNo) {
        tempSum += (*it).rate;
        it++;
      }
      if (tempSum < randNo) {
        // should not happen
        std::cerr << "Could not pick event" << std::endl;
        return false;
      }
      it--;
            
      if (verbose >= 2) {
        std::cout << "Vertex #" << v << " changes state: " 
                  << model.getVertexStates()[graph[v].state]
                  << "->" << model.getVertexStates()[(*it).newState]
                  << std::endl;
      }

      // collect statistics
      if (model.isInfection(graph[v].state, it->newState)) addInfection();
      else if (model.isInformation(graph[v].state, it->newState)) addInformation();
      else if (model.isRecovery(graph[v].state, it->newState)) addRecovery();
      else if (model.isForgetting(graph[v].state, it->newState)) addForgetting();

      // process the change of state
      graph[v].state = (*it).newState;
            
      // update vertex event list
      double rateDiff = generateEventList(graph, v, model, verbose);
            
      // update sum of rates for the vertex
      vertex_index_type index = get(boost::vertex_index, graph);
      tree.getLeaves()[index[v]]->updateRateSum(rateDiff);
            
      //update neighbours
      out_edge_iterator oi, oi_end;;
   
      for (tie(oi, oi_end) = boost::out_edges(v, graph);
           oi != oi_end; ++oi) {
        rateDiff = updateEventList(graph, *oi, model, verbose);
        tree.getLeaves()[index[target(*oi, graph)]]->updateRateSum(rateDiff);
      }
      return true;
    } else {
      // something went wrong in picking the vertex, should not happen
      return false;
    }
  }

  //----------------------------------------------------------
  /*! \brief Print the state of the simulation. 
  
  Prints the state of all vertices, as well as the events that can happen to
  them and at which rate.
  */
  template <typename RandomGenerator, typename Graph>
  void GillespieSimulator<RandomGenerator, Graph>::print()
  {
    // get simulation variable
    const Model& model = this->getModel();
    
    std::vector<Tree::Leaf<unsigned int>*>::iterator it;
    for (it = tree.getLeaves().begin(); it != tree.getLeaves().end(); it++) {
      if ((*it)->getRateSum() > 0) {
        std::cout << "Vertex #" << *((*it)->getItem()) << " ["
                  << model.getVertexState(graph[*((*it)->getItem())].state)
                  << "]:";
        eventList::iterator eit;
        for (eit = graph[*((*it)->getItem())].events.begin();
             eit != graph[*((*it)->getItem())].events.end(); eit++) {
          std::cout << " " << model.getVertexState((*eit).newState) << " ("
                    << (*eit).rate << ")";
        }
        std::cout << std::endl;
      }
    }
  }
  
}

#endif
