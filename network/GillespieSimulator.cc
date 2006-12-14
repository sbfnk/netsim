/*******************************************************************/
//
// Class GillespieSimulator
// --------------------
//
// This is a graph implementation using the Gillespie algorithm
// for event selection. It contains a tree structure for quick access
// of randomly chosen events, and a boost graph for the adjacency
// structure
//
/******************************************************************/

#include <iostream>

#include <boost/random.hpp>

#include "RandomGenerator.hh"
#include "GillespieSimulator.hh"

/******************************************************************/
// GillespieSimulator constructor
// creates the random generator for event selection
/******************************************************************/
GillespieSimulator::GillespieSimulator()
   : randGen(new mt19937Uniform()),
     time(0.)
{
}

/******************************************************************/
// GillespieSimulator destructor
/******************************************************************/
GillespieSimulator::~GillespieSimulator()
{
}

/******************************************************************/
// GillespieSimulator::initialize
// initializes the graph with the rates for the possible processes
/******************************************************************/
void GillespieSimulator::initialize(const Model& model)
{
   vertex_iterator vi, vi_end;
   for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      generateEventList(graph, *vi, model);
   }
}

/******************************************************************/
// GillespieSimulator::updateState
// updates the state of the graph by advancing one time step
// and choosing an event to process
/******************************************************************/
bool GillespieSimulator::updateState(const Model& model)
{
   // exit if nothing can happen
   if (tree.getTopBin()->getRateSum() ==0) {
      std::cout << "Nothing can happen" << std::endl;
      return false;
   }

   // draw a random number from [0,1) for the timestep advance
   double randNo = (*randGen)();
   randNo *= tree.getTopBin()->getRateSum();
   time += (1/randNo);

   // draw another random number from [0,1) for picking the event
   randNo = (*randGen)();
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
         std::cout << "Could not pick event" << std::endl;
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
      
   }
   return true;
}

