/*! \file ChrisSimulator.hh
  \brief The Simulators::ChrisSimulator class.
*/

#ifndef CHRISSIMULATOR_HH
#define CHRISSIMULATOR_HH

#include <queue>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "Tree.hh"
#include "Simulator.hh"

namespace Simulators {
  
  /*! \brief The Chris simulation (dysfunctional)
    
  This contains an unfinished attempt to design a Chris simulation, i.e. a
  simulation storing the events in a priority queue and allowing for arbitrary
  distributions of interevent times by immediately assigning each event a time
  at which it will happen according to that distriubtion, as proposed by
  Watkins(2006). 
  
  */
  template <typename RandomGenerator, typename Graph>
  class ChrisSimulator :
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
    
  public:
    
    /*! \brief Constructor
      
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    ChrisSimulator(RandomGenerator& r, Graph& g, Model& m,
                   unsigned int v = 0) :
      Simulator(m, v), randGen(r), graph(g) {;}
    ~ChrisSimulator() {;}
    
    void initialize();
    bool updateState();
    
    void print();
    
  private:

    //! \brief An event in the Chris simulator
    struct ChrisEvent
    {
      State newState; //! The new state after the event happens.
      vertex_descriptor vertex; //! The affected vertex.
      double eventTime; //! The time at which the event will happen.
      bool valid; //! Whether the event can still happen.
      
      /*! \brief The operator < for Chris events
      Stores elements in reverse order --> lowest time has highest priority.
      \param[in] rhs The right hand side of the comparison
      */
      bool operator<(const ChrisEvent& rhs) const
      { return eventTime>rhs.eventTime; }
      
    };
    
    //! The random generator to be used for choosing event times
    RandomGenerator& randGen;
    Graph& graph; //!< The graph determining how vertices can effect another
    
    std::priority_queue<ChrisEvent*> events; //!< The queue holding the events.
    // A map from vertices to the events that can happen to them
    std::multimap<vertex_descriptor, ChrisEvent*> eventPtrs;

    void generateEvents(vertex_descriptor v);
    
  };
  
  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
    
  Initialises the graph with the rates for the possible processes.
  */
  template <typename RandomGenerator, typename Graph>
  void ChrisSimulator<RandomGenerator, Graph>::initialize()
  {
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      generateEvents(*vi);
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Perform an update and process one event.

  Chooses the next possible event and processes that event. After that, the
  priority queue shoudl be updated
  */
  template <typename RandomGenerator, typename Graph>
  bool ChrisSimulator<RandomGenerator, Graph>::updateState()
  {
//    bool found = false;
//     while (!found && !events.empty()) {
//       ChrisEvent* ev;  = events.pop();
//       if (ev->valid) {
//         Simulator::updateTime(ev->eventTime);
        
//         found = true;
//       }
//     }
    return true;
  }
  
  //----------------------------------------------------------
  /*! \brief Print the state of the simulation. 
  
  Prints the priority queue.
  */
  template <typename RandomGenerator, typename Graph>
  void ChrisSimulator<RandomGenerator, Graph>::print()
  {
    const Model& model = this->getModel();

    ChrisEvent* ev = events.top();
    
    std::cout << "Next event: Vertex #" << ev->vertex << " ["
              << model.getVertexStates()[graph[ev->vertex].state.base]
              << " --> " << model.getVertexStates()[ev->newState.base]
              << "]" << std::endl;
  }
  
  //----------------------------------------------------------
  /*! \brief Generate events.
    
  Generates the queue of events.
  */
  template <typename RandomGenerator, typename Graph>
  void ChrisSimulator<RandomGenerator, Graph>::generateEvents(vertex_descriptor v)
  {
    const Model& model = this->getModel();

    generateEventList(graph, v, model);
    
    for (eventList::iterator it = graph[v].events.begin();
         it != graph[v].events.end(); it++) {
      // generate time
      boost::exponential_distribution<double> exp_dist(it->rate);
      boost::variate_generator<RandomGenerator&, boost::exponential_distribution<double> >
        exp_sampler(randGen, exp_dist);
      
      ChrisEvent* newEvent = new ChrisEvent;
      newEvent->newState = it->newState;
      newEvent->vertex = v;
      newEvent->eventTime = exp_sampler();
      newEvent->valid = true;
      
      events.push(newEvent);
      eventPtrs.insert(std::make_pair(v, newEvent));
    }
  }

}

#endif
