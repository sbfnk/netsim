/*! \file GillespieSimulator.hh
  \brief The Simulators::GillespieSimulator class.
*/

#ifndef GILLESPIESIMULATOR_HH
#define GILLESPIESIMULATOR_HH

#include <math.h>
#include <iomanip>

#include "Tree.hh"
#include "Simulator.hh"
#include "Vertex.hh"

#include "InfoSIRS.hh"
#include "DimInfoSIRS.hh"
#include "VaccinationSIRS.hh"
#include "ProtectiveSIRS.hh"
#include "SingleSIRS.hh"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//! \addtogroup epi_simulator Epidemic simulator
//! \addtogroup gillespie_simulator Gillespie simulator

namespace Simulators {
  
  /*! \brief The Gillespie simulation
    
  This contains the Gillespie simulation, able to run a simulation on a graph
  using a given model, with a sequential update rule according to
  Gillespie (1977). The rates of the events are stored in a Tree for quick
  access, which is generated using the member function initialise. The actual
  updating of the graph is done in updateState.
  
  \ingroup epi_simulator
  \ingroup gillespie_simulator
  */
  template <typename RandomGenerator, typename Graph>
  class GillespieSimulator :
    public Simulator<Graph>
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
    typedef typename boost::variate_generator
    <RandomGenerator&, boost::uniform_real<> > uniform_gen;

  public:

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    GillespieSimulator(RandomGenerator& r, Graph& g,
                       unsigned int v = 0);
    virtual ~GillespieSimulator() {;}
      
    void initialise();
    bool updateState();
    void changeVertexState(vertex_descriptor v, vertex_descriptor nb,
                           State* before, State* after);
    virtual void updateEventStats(State* before, State* after,
                                  vertex_descriptor v, vertex_descriptor nb) {;}
  
    virtual bool stopCondition() const
    { return Simulator<Graph>::stopCondition(); }

    void print();

  private:
  
    uniform_gen randGen; //!< The random generator to be used for choosing events
    Tree::Tree<unsigned int> tree; //!< The tree holding the event rates

  };

  template <typename RandomGenerator, typename Graph>
  GillespieSimulator<RandomGenerator, Graph>::
  GillespieSimulator(RandomGenerator& r, Graph& g, unsigned int v) :
    Simulator<Graph>(g, v), randGen(r, boost::uniform_real<> (0,1))
  {
    this->knownModels.push_back
      (std::make_pair("DimInfoSIRS", new Models::DimInfoSIRS<Graph>(this->getVerbose())));
    this->knownModels.push_back
      (std::make_pair("InfoSIRS", new Models::InfoSIRS<Graph>(this->getVerbose())));
    this->knownModels.push_back
      (std::make_pair("VaccinationSIRS", new Models::VaccinationSIRS<Graph>(this->getVerbose())));
    this->knownModels.push_back
      (std::make_pair("ProtectiveSIRS", new Models::ProtectiveSIRS<Graph>(this->getVerbose())));
    this->knownModels.push_back
      (std::make_pair("SingleSIRS", new Models::SingleSIRS<Graph>(this->getVerbose())));
  }
  
  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
  
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  void GillespieSimulator<RandomGenerator, Graph>::initialise()
  {
    Simulator<Graph>::initialise();

    // get simulation variables
    Graph& graph = this->getGraph();
    unsigned int verbose = this->getVerbose();
    const Model<Graph>* model = Simulator<Graph>::getModel();

    // generate event list for each vertex
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi) {
      generateEventList(graph, *vi, *model, verbose);
    }

    // generate the tree holding rate sums corresponding to vertices
    generateTree(tree,graph,
                 get(&vertex_property_type::rateSum,graph),
                 get(boost::vertex_index,graph));
  }

  
  template <typename RandomGenerator, typename Graph>
  void GillespieSimulator<RandomGenerator, Graph>::
  changeVertexState(vertex_descriptor v, vertex_descriptor nb,
                    State* before, State* after)
  {
    unsigned int verbose = this->getVerbose();
    const Model<Graph>* model = this->getModel();
    Graph& graph = this->getGraph();

    if (verbose >= 2) {
      std::cout << "Vertex #" << v << " changes state: (";
      std::cout << model->printState(before)
                << ")->("
                << model->printState(after)
                << ")";
      if (nb == v) {
        std::cout << std::endl;
      } else {
        std::cout << ", induced by " << nb << std::endl;
      }
    }

    // process the change of state
    graph[v].state = after->clone();

    // update stats
    updateEventStats(before, after, v, nb);
    // process the change of state
    delete before;

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
    unsigned int verbose = this->getVerbose();
    Graph& graph = this->getGraph();
    const Model<Graph>* model = this->getModel();

    // exit if nothing can happen
    if (tree.getTopBin()->getRateSum() < 1) {
      if (verbose >=2) {
        std::cout << "Nothing can happen. Stopping run. " << std::endl;
      }
      this->updateStats(true);
      return false;
    }
  
    if (verbose >= 2) {
      std::cout << "choose event, total sum of rates is "
                << tree.getTopBin()->getRateSum()/1e+4 << std::endl;
      print();
    }
         
    // draw a random number from [0,1) for the timestep advance
    double randTime = (randGen)();
    this->updateTime(-log(randTime)*1e+4/tree.getTopBin()->getRateSum());
         
    // draw another random number from [0,rateSum) for picking the event
    unsigned int randEvent =
      static_cast<unsigned int>((randGen)() * tree.getTopBin()->getRateSum())+1;
    unsigned int* eventVertex = tree.pickRandomElement(randEvent);
    if (eventVertex) {
      // process vertex event
      unsigned int tempSum = 0;
      vertex_descriptor v = vertex(*eventVertex, graph);

      // find event corresponding to randNo
      std::vector<Event>::iterator it = graph[v].events.begin();
      while (it != graph[v].events.end() && tempSum < randEvent) {
        tempSum += (*it).rate;
        it++;
      }
      if (tempSum < randEvent) {
        // should not happen
        std::cerr << "Could not pick event" << std::endl;
        return false;
      }
      it--;
            
      changeVertexState(v, it->nb, graph[v].state, it->newState);

      // update vertex event list
      unsigned int rateDiff = generateEventList(graph, v, *model, verbose);
            
      // update sum of rates for the vertex
      vertex_index_type index = get(boost::vertex_index, graph);
      tree.getLeaves()[index[v]]->updateRateSum(rateDiff);
            
      //update neighbours
      out_edge_iterator oi, oi_end;;
   
      for (tie(oi, oi_end) = boost::out_edges(v, graph);
           oi != oi_end; ++oi) {
        rateDiff = updateEventList(graph, *oi, *model, verbose);
        tree.getLeaves()[index[target(*oi, graph)]]->updateRateSum(rateDiff);
      }
      return true;
    } else {
      // something went wrong in picking the vertex, should not happen
      std::cout << "something went wrong in picking the vertex" << std::endl;
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
    Graph& graph = this->getGraph();
    const Model<Graph>* model = Simulator<Graph>::getModel();

    std::cout << std::endl;

    std::vector<Tree::Leaf<unsigned int>*>::iterator it;
    for (it = tree.getLeaves().begin(); it != tree.getLeaves().end(); it++) {
      if ((*it)->getRateSum() > 0) {
        std::cout << "Vertex #" << *((*it)->getItem()) << " ["
                  << model->printState(graph[*((*it)->getItem())].state)
                  << "]:";
        eventList::iterator eit;
        for (eit = graph[*((*it)->getItem())].events.begin();
             eit != graph[*((*it)->getItem())].events.end(); eit++) {
          std::cout << " " << model->printState(eit->newState) << " ("
                    << (*eit).rate/1e+4;
          if (eit->nb != *((*it)->getItem())) {
            std::cout << "," << eit->nb;
          }
          std::cout << ")";
        }
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  }
  
}

//----------------------------------------------------------
/*! \brief Generate event list for a vertex.

Generates an event list for a given vertex. This collects both events which can
happen to a node by itself, as well as events transmitted over edges.

\param[in] graph The graph variable defining the edge structure.
\param[in] v The vertex to consider.
\param[in] model The model to use to generate the events.
\param[in] verbose The level of verbosity.
\return The sum of rates for all events which have been generated.
\ingroup gillespie_simulator
*/
template <class Graph, class Model>
unsigned int generateEventList(Graph& graph,
                     typename boost::graph_traits<Graph>::vertex_descriptor v,
                               const Model& model,
                               unsigned int verbose = 0)
{
  // definitions of boost types for quick access
  typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
  typedef typename boost::graph_traits<Graph>::edge_descriptor
    edge_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;

  typename boost::property_map<Graph, boost::vertex_index_t>::type 
    id = get(boost::vertex_index, graph);
   
  // temporary sum for the new sum of rates of all events
  // that can affect the vertex v
  unsigned int tempSum = 0;

  if (verbose >= 2) {
    std::cout << "Generating events list for vertex #" << v << " ("
              << model.printState(graph[v].state) << ")"
              << std::endl;
  }

  // clear event list
  graph[v].events.clear();
   
  // get node events
  tempSum += model.getNodeEvents(graph[v].events, graph[v].state, id[v]);
   
  // get edge events
  out_edge_iterator oi, oi_end;;
   
  for (tie(oi, oi_end) = boost::out_edges(v, graph);
       oi != oi_end; ++oi) {
    edge_descriptor e = *oi;
    vertex_descriptor t =  target(e, graph);
    tempSum +=
      model.getEdgeEvents(graph[v].events, graph[v].state, graph[e].type,
                          graph[t].state, id[t]);
     
  }
   
  // calculate difference between new and old sum of rates
  unsigned int diff = tempSum - graph[v].rateSum;
  // set rateSum to new value
  graph[v].rateSum = tempSum;
   
  return diff;
}

//----------------------------------------------------------
/*! \brief Update event list for a vertex.

Updates an event list generated by the source of a edge for the target of the
same edge. The events happening over that edge which were previously in the list
are deleted and a new list of these events is generated

\param[in] graph The graph variable defining the edge structure.
\param[in] e The edge to consider.
\param[in] model The model to use to generate the events.
\param[in] verbose The level of verbosity.
\return The sum of rates for all events which have been generated.
\ingroup gillespie_simulator
*/
template <class Graph, class Model>
unsigned int updateEventList(Graph& graph,
                             typename boost::graph_traits<Graph>::edge_descriptor e,
                             const Model& model,
                             unsigned int verbose = 0)
{
  typename boost::graph_traits<Graph>::vertex_descriptor v = target(e, graph);
  typename boost::graph_traits<Graph>::vertex_descriptor n = source(e, graph);
  
  typename boost::property_map<Graph, boost::vertex_index_t>::type 
    id = get(boost::vertex_index, graph);
   
  // temporary sum for the new sum of rates of all events
  // that can affect the vertex v
  unsigned int tempSum = 0;

  if (verbose >= 2) {
    std::cout << "Updating events list for vertex #" << v << " ("
              << model.printState(graph[v].state);
    std::cout << "), as generated by neighbour #" << n << " ("
              << model.printState(graph[n].state);
    std::cout << ") along "
              << model.getEdgeTypes()[graph[e].type] << "-edge" << std::endl;
  }
  for (unsigned int i = 0; i < graph[v].events.size(); i++) {
    if (graph[v].events[i].nb == id[n] &&
        graph[v].events[i].et == graph[e].type) {
      tempSum -= graph[v].events[i].rate;
      graph[v].events.erase(graph[v].events.begin() + i);
      --i;
    }
  }

  tempSum +=
    model.getEdgeEvents(graph[v].events, graph[v].state,
                        graph[e].type, graph[n].state, id[n]);
  
  // update rateSum
  graph[v].rateSum += tempSum;
    
  return tempSum;
}

#endif
