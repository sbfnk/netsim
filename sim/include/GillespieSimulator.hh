#ifndef GILLESPIESIMULATOR_HH
#define GILLESPIESIMULATOR_HH

#include "Tree.hh"
#include "Simulator.hh"

template <typename RandomGenerator, typename Graph>
class GillespieSimulator :
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

  GillespieSimulator(RandomGenerator r, Graph& g, Model& m, unsigned int v = 0) :
    Simulator(m, v), randGen(r), graph(g) {;}
  ~GillespieSimulator() {;}
      
  void initialize();
  bool updateState();
  
  void print();

private:
  
  uniform_gen randGen; // random generator
  Graph& graph;
  Tree<unsigned int> tree; // contains the tree structure
  
};

/******************************************************************/
// GillespieSimulator<RandomGenerator>::initialize
// initializes the graph with the rates for the possible processes
/******************************************************************/
template <typename RandomGenerator, typename Graph>
void GillespieSimulator<RandomGenerator, Graph>::initialize()
{
  vertex_iterator vi, vi_end;
  for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
    generateEventList(graph, *vi, model, verbose);
  }

  generateTree(tree,graph,
               get(&vertex_property_type::rateSum,graph),
               get(boost::vertex_index,graph));
}

/******************************************************************/
// GillespieSimulator<RandomGenerator>::updateState
// updates the state of the graph by advancing one time step
// and choosing an event to process
/******************************************************************/
template <typename RandomGenerator, typename Graph>
bool GillespieSimulator<RandomGenerator, Graph>::updateState()
{
  typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;

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
            
    // process the change of state
    if (verbose >= 2) {
      std::cout << "Vertex #" << v << " changes state: " 
                << model.getVertexStates()[graph[v].state]
                << "->" << model.getVertexStates()[(*it).newState] << std::endl;
    }

    if (model.isInfection(graph[v].state, it->newState)) addInfection();
    else if (model.isInformation(graph[v].state, it->newState)) addInformation();
    else if (model.isRecovery(graph[v].state, it->newState)) addRecovery();
    else if (model.isForgetting(graph[v].state, it->newState)) addForgetting();
    
    graph[v].state = (*it).newState;
            
    // update vertex event list
    double rateDiff = generateEventList(graph, v, model, verbose);
            
    // update tree
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
    return false;
  }
}

template <typename RandomGenerator, typename Graph>
void GillespieSimulator<RandomGenerator, Graph>::print()
{
  std::vector<Leaf<unsigned int>*>::iterator it;
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

#endif
