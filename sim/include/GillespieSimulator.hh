#ifndef GILLESPIESIMULATOR_HH
#define GILLESPIESIMULATOR_HH

#include "Tree.hh"
#include "Model.hh"
#include "Vertex.hh"

class RandomGenerator;

// define gillespie_graph derived from boost adjacency list
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> gillespie_graph;

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

   public:
      
      GillespieSimulator();
      ~GillespieSimulator();
      
      void initialize(const Model& model);
      bool updateState(const Model& model);

      double getTime() { return time; };

      gillespie_graph graph; // contains the adjacency structure
      Tree<unsigned int> tree; // contains the tree structure

   private:

      RandomGenerator* randGen; // random generator
      double time;
      
};

#endif
