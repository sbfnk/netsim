/*! \file Simulator.hh
  \brief The Simulator class.
*/
#ifndef SIMULATOR_HH
#define SIMULATOR_HH

#include "Model.hh"

//! \addtogroup gillespie_simulator Gillespie simulator

/*! \brief Base class for simulation classes
  
Classes derived from this class implement a given simulation algorithm, e.g. the
Gillespie algorithm or, possibly, other algorithms with different waiting times
etc. Most importantly, this class provides a declaration for the updateState and
getTime methods used by the main simulation code.

\ingroup gillespie_simulator
*/
class Simulator
{

public:

  /*! \brief Constructor.
  \param[in] m model initialiser
  \param[in] v verbose intialiser
  */
  Simulator(const Model& m, unsigned int v = 0) :
    model(m), verbose(v), time(0.),
    numInfections(0), numInformations(0), numRecoveries(0), numForgettings(0)
  {}

  //! Destructor.
  virtual ~Simulator() {;}

  //! Initialise the simulation. Implemented by derived classes.
  virtual void initialise()
  {
    numInfections = numInformations = numRecoveries = numForgettings = 0;
  }
  
  //! Perform an update. Implemented by derived classes.
  virtual bool updateState() = 0;

  //! Print the state of the simulation. Implemented by derived classes.
  virtual void print() {;}

  //! Accessor for the time variable.
  double getTime() const { return time; };
  /*! \brief Update the current time in the simulation

  \param[in] t The timestep to advance the time by.
  */
  void updateTime(double t) { time += t; };

  //! Accessor for the numInfections variable.
  unsigned int getNumInfections() const { return numInfections; }
  //! Accessor for the numRecoveries variable.
  unsigned int getNumRecoveries() const { return numRecoveries; }
  //! Accessor for the numInformations variable.
  unsigned int getNumInformations() const { return numInformations; }
  //! Accessor for the numForgettings variable.
  unsigned int getNumForgettings() const { return numForgettings; }

  //! Increase the numInfections variable by one.
  void addInfection() { ++numInfections; }
  //! Increase the numInformations variable by one.
  void addInformation() { ++numInformations; }
  //! Increase the numRecoveries variable by one.
  void addRecovery() { ++numRecoveries; }
  //! Increase the numInfections variable by one.
  void addForgetting() { ++numForgettings; }

  //! Accessor for the model variable.
  const Model& getModel() const { return model; }
  //! Accessor for the verbose variable.
  unsigned int getVerbose() const { return verbose; }
  
private:

  const Model& model; //!< The model to be used by the simulation.
  unsigned int verbose; //!< The verbosity level.

  double time; //!< The current time of the simulation.
  unsigned int numInfections; //!< A counter for the number of infections.
  unsigned int numInformations; //!< A counter for the number of informations.
  unsigned int numRecoveries; //!< A counter for the number of recoveries.
  unsigned int numForgettings; //!< A counter for the number of forgettings.
      
};

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
double generateEventList(Graph& graph,
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
  double tempSum = .0;

  if (verbose >= 2) {
    std::cout << "Generating events list for vertex #" << v << " ("
              << model.getVertexStates()[graph[v].state] << ")" << std::endl;
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
      model.getEdgeEvents(graph[v].events, graph[v].state,
                          graph[v].state_detail, graph[e].type, graph[t].state,
                          graph[t].state_detail, id[t]);
     
  }
   
  // calculate difference between new and old sum of rates
  double diff = tempSum - graph[v].rateSum;
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
double updateEventList(Graph& graph,
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
  double tempSum = .0;

  if (verbose >= 2) {
    std::cout << "Updating events list for vertex #" << v << " ("
              << model.getVertexStates()[graph[v].state];
    if (graph[v].state_detail > 0.) {
      std::cout << "," << std::setprecision(1) << std::fixed << graph[v].state_detail;
    }
    std::cout << "), as generated by neighbour #" << n << " ("
              << model.getVertexStates()[graph[n].state];
    if (graph[n].state_detail > 0.) {
      std::cout << "," << std::setprecision(1) << std::fixed << graph[n].state_detail;
    }
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
    model.getEdgeEvents(graph[v].events, graph[v].state, graph[v].state_detail,
                        graph[e].type,   graph[n].state, graph[n].state_detail,
                        id[n]);
  
  // update rateSum
  graph[v].rateSum += tempSum;
    
  return tempSum;
}

//----------------------------------------------------------
/*! \brief Simulation algorithms.
*/
namespace Simulators {} // just define the namespace

#endif
