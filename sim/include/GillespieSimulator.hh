#ifndef GILLESPIESIMULATOR_HH
#define GILLESPIESIMULATOR_HH

#include "Tree.hh"
#include "Model.hh"
#include "Vertex.hh"

// define gillespie_graph derived from boost adjacency list
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> gillespie_graph;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> onetype_graph;

template <typename RandomGenerator>
class GillespieSimulator
{

  // type definitions for easy access of boost graph properties
  typedef boost::graph_traits<gillespie_graph>::vertex_descriptor
  vertex_descriptor;
  typedef boost::graph_traits<gillespie_graph>::vertex_iterator
  vertex_iterator;
  typedef boost::graph_traits<gillespie_graph>::adjacency_iterator
  adjacency_iterator;
  typedef boost::property_map<gillespie_graph, boost::vertex_index_t>::type
  vertex_index_type;
  typedef boost::uniform_01<RandomGenerator, double> uniform_gen;

public:
      
  GillespieSimulator(RandomGenerator& rg)
    : randGen(uniform_gen(rg)),time(0.) {}
  ~GillespieSimulator();
      
  /******************************************************************/
  // GillespieSimulator<RandomGenerator>::initialize
  // initializes the graph with the rates for the possible processes
  /******************************************************************/
  void initialize(const Model& model)
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
  bool updateState(const Model& model)
  {
    // exit if nothing can happen -- exit also if the rate sum is very small
    // to prevent rounding problems
    if (tree.getTopBin()->getRateSum() < 1e-10) {
      return false;
    }
         
    // draw a random number from [0,1) for the timestep advance
    double randNo = (randGen)();
    time += - log(randNo)/tree.getTopBin()->getRateSum();
         
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
      tree.leaves[index[v]]->updateRateSum(rateDiff);
            
      //update neighbours
      adjacency_iterator ai, ai_end;
      for (tie(ai, ai_end) = adjacent_vertices(v, graph);
           ai != ai_end; ai++) {
        rateDiff =
          generateEventList(graph, *ai, model);
               
        tree.leaves[index[*ai]]->updateRateSum(rateDiff);
      }
      return true;
    } else {
      return false;
    }
  }
      
  double getTime() { return time; };

  gillespie_graph graph; // contains the adjacency structure
  Tree<unsigned int> tree; // contains the tree structure

private:

  uniform_gen randGen; // random generator
  double time;
      
};

#endif
