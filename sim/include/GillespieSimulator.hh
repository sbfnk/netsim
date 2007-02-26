#ifndef GILLESPIESIMULATOR_HH
#define GILLESPIESIMULATOR_HH

#include "Tree.hh"
#include "Simulator.hh"

template <typename RandomGenerator, typename Graph, typename Model>
class GillespieSimulator :
  virtual public Simulator<Model>
{
  
  // type definitions for easy access of boost graph properties
  typedef boost::graph_traits<dualtype_graph>::vertex_descriptor
  vertex_descriptor;
  typedef boost::graph_traits<dualtype_graph>::vertex_iterator
  vertex_iterator;
  typedef boost::graph_traits<dualtype_graph>::adjacency_iterator
  adjacency_iterator;
  typedef boost::property_map<dualtype_graph, boost::vertex_index_t>::type
  vertex_index_type;
  typedef boost::uniform_01<RandomGenerator, double> uniform_gen;

public:

  GillespieSimulator(RandomGenerator r, Graph& g) : randGen(r), graph(g) {;}
  ~GillespieSimulator() {;}
      
  void initialize(const Model& model);
  bool updateState(const Model& model);
  
  void print();

private:
  
  uniform_gen randGen; // random generator
  Graph graph;
  Tree<unsigned int> tree; // contains the tree structure

};

/******************************************************************/
// GillespieSimulator<RandomGenerator>::initialize
// initializes the graph with the rates for the possible processes
/******************************************************************/
template <typename RandomGenerator, typename Graph, typename Model>
void GillespieSimulator<RandomGenerator, Graph, Model>::initialize
(const Model& model)
{
  vertex_iterator vi, vi_end;
  for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
    generateEventList(graph, *vi, model);
  }
}

/******************************************************************/
// GillespieSimulator<RandomGenerator>::updateState
// updates the state of the graph by advancing one time step
// and choosing an event to process
/******************************************************************/
template <typename RandomGenerator, typename Graph, typename Model>
bool GillespieSimulator<RandomGenerator, Graph, Model>::updateState
(const Model& model)
{
  // exit if nothing can happen -- exit also if the rate sum is very small
  // to prevent rounding problems
  if (tree.getTopBin()->getRateSum() < 1e-10) {
    return false;
  }
         
  // draw a random number from [0,1) for the timestep advance
  double randNo = (randGen)();
  Simulator<Model>::updateTime(-log(randNo)/tree.getTopBin()->getRateSum());
         
  // draw another random number from [0,1) for picking the event
  randNo = (randGen)();
  unsigned int* eventVertex = tree.pickRandomElement(randNo);
  if (eventVertex) {
    // process vertex event
    double tempSum = .0;
    vertex_descriptor v = vertex(*eventVertex, graph);
            
    std::list<event>::iterator it = graph[v].events.begin();
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
            
    // process the change of state
    graph[v].state = (*it).newState;
            
    // update vertex event list
    double rateDiff = generateEventList(graph, v, model);
            
    // update tree
    vertex_index_type index = get(boost::vertex_index, graph);
    tree.getLeaves()[index[v]]->updateRateSum(rateDiff);
            
    //update neighbours
    adjacency_iterator ai, ai_end;
    for (tie(ai, ai_end) = adjacent_vertices(v, graph);
         ai != ai_end; ai++) {
      rateDiff =
        generateEventList(graph, *ai, model);
               
      tree.getLeaves()[index[*ai]]->updateRateSum(rateDiff);
    }
    return true;
  } else {
    return false;
  }
}

template <typename RandomGenerator, typename Graph, typename Model>
void GillespieSimulator<RandomGenerator, Graph, Model>::print
()
{
  std::vector<Leaf<unsigned int>*>::iterator it;
  for (it = tree.getLeaves().begin(); it != tree.getLeaves().end(); it++) {
    if ((*it)->getRateSum() > 0) {
      std::cout << "Vertex #" << *((*it)->getItem()) << " ["
                << graph[*((*it)->getItem())].state << "]:";
      eventList::iterator eit;
      for (eit = graph[*((*it)->getItem())].events.begin();
           eit != graph[*((*it)->getItem())].events.end(); eit++) {
        std::cout << " " << (*eit).newState << " ("
                  << (*eit).rate << ")";
      }
      std::cout << std::endl;
    }
  }
}

#endif
