#ifndef GILLESPIEGRAPH_HH
#define GILLESPIEGRAPH_HH

#include "Tree.hh"
#include "Model.hh"
#include "Vertex.hh"

class RandomGenerator;

// define gillespie_graph derived from boost adjacency list
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> gillespie_graph;

class GillespieGraph
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
      
      GillespieGraph();
      ~GillespieGraph();
      
      void initialize(const Model& model);
      void updateState(const Model& model);

      gillespie_graph graph; // contains the adjacency structure
      Tree<unsigned int> tree; // contains the tree structure

   private:

      RandomGenerator* randGen; // random generator
      
};

#endif
