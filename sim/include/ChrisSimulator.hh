#ifndef CHRISSIMULATOR_HH
#define CHRISSIMULATOR_HH

#include <queue>
#include <boost/random/exponential_distribution.hpp>

#include "Tree.hh"
#include "Simulator.hh"

template <typename RandomGenerator, typename Graph>
class ChrisSimulator :
  virtual public Simulator
{
  
  // type definitions for easy access of boost graph properties
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

  ChrisSimulator(RandomGenerator r, Graph& g, Model& m) :
    Simulator(m), randGen(r), graph(g) {;}
  ~ChrisSimulator() {;}
      
  void initialize();
  bool updateState();
  
  void print();

private:

  struct ChrisEvent
  {
    int newState;
    vertex_descriptor vertex;
    double eventTime;
    bool valid;

    // store elements in reverse order --> lowest time has highest priority
    bool operator<(const ChrisEvent& rhs) const
    { return eventTime>rhs.eventTime; }

  };

  RandomGenerator& randGen;
  Graph& graph;
  
  std::priority_queue<ChrisEvent*> events;
  std::multimap<vertex_descriptor, ChrisEvent*> eventPtrs;

  void generateEvents(vertex_descriptor v);
  
};

/******************************************************************/
// ChrisSimulator<RandomGenerator>::initialize
// initializes the graph with the rates for the possible processes
/******************************************************************/
template <typename RandomGenerator, typename Graph>
void ChrisSimulator<RandomGenerator, Graph>::initialize()
{
  vertex_iterator vi, vi_end;
  for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
    generateEvents(*vi);
  }
}

/******************************************************************/
// ChrisSimulator<RandomGenerator>::updateState
// updates the state of the graph by advancing one time step
// and choosing an event to process
/******************************************************************/
template <typename RandomGenerator, typename Graph>
bool ChrisSimulator<RandomGenerator, Graph>::updateState()
{
  bool found = false;
  while (!found && !events.empty()) {
    ChrisEvent* ev = events.pop();
    if (ev->valid) {
      Simulator::updateTime(ev->time);
      
      found = true;
    }
  }
  return true;
}

template <typename RandomGenerator, typename Graph>
void ChrisSimulator<RandomGenerator, Graph>::print()
{
  ChrisEvent* ev = events.top();
  
  std::cout << "Next event: Vertex #" << ev->vertex << " [" <<
            << model.getVertexStates()[graph[ev->vertex].state] << " --> "
            << model.getVertexStates()[ev->newState] << "]" << std::endl;
}

template <typename RandomGenerator, typename Graph>
void ChrisSimulator<RandomGenerator, Graph>::generateEvents(vertex_descriptor v)
{
  generateEventList(graph, v, model);

  for (eventList::iterator it = graph[v].events.begin();
       it != graph[v].events.end(); it++) {
    // generate time
    boost::exponential_distribution<double> exp_dist(it->rate);
    boost::variate_generator
      <RandomGenerator&, boost::exponential_distribution<double> >
      exp_sampler(randGen, exp_dist);

    ChrisEvent* newEvent = new ChrisEvent;
    newEvent->newState = it->newState;
    newEvent->vertex_descriptor = v;
    newEvent->time = exp_sampler();
    newEvent->valid = true;

    events.push(newEvent);
    eventPtrs.insert(make_pair(v, newEvent));
  }
}

  
#endif
