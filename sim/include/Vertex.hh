/******************************************************************/
// Vertex.hh
// contains the Edge and Vertex classes, to be used by the different
// graph implementations
/******************************************************************/

#ifndef VERTEX_HH
#define VERTEX_HH

#include <list>

#include "Model.hh"

class Edge {
   public:
      Edge(): type(Disease) {;}
      Edge(EdgeType eType): type(eType) {;}
      EdgeType type;
};

class Vertex {

   public:

      Vertex() : state(VertexState(Susceptible, Uninformed)) {;}
      Vertex(int diseaseState, int infoState)
         : state(VertexState(diseaseState, infoState)) {;}

      Vertex(VertexState vs)
         : state(vs) {;}
      
      virtual ~Vertex() {;}

      VertexState state; // present (disease and information) state
      double rateSum; // sum of the rates of all events that can
                      // affect the vertex (change its state)

      eventList events; // list of events associated with vertex

};


/******************************************************************/
// generateEventList function
// generates an event list for a given vertex in a given graph,
// using a given model. The vertices of the event must be of type
// Vertex, otherwise this will throw a compile error
/******************************************************************/
template <class Graph, class VertexClass>
double generateEventList(Graph& graph, VertexClass v,
                         const Model& model)
{
   // definitions of boost types for quick access
   typedef typename boost::graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
   typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge_descriptor;
   typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;

   // temporary sum for the new sum of rates of all events
   // that can affect the vertex v
   double tempSum = .0;

   // clear event list
   graph[v].events.clear();
   // get node events
   tempSum += model.getNodeEvents(graph[v].events, graph[v].state);

   // get edge events
   out_edge_iterator oi, oi_end;;

//    unsigned int out_degrees;
   
//    for (tie(oi, oi_end) = boost::out_edges(v, graph);
//         oi != oi_end; ++oi) {
//    }
   
   for (tie(oi, oi_end) = boost::out_edges(v, graph);
        oi != oi_end; ++oi) {
      edge_descriptor e = *oi;
      vertex_descriptor t =  target(e, graph);
      tempSum +=
         model.getEdgeEvents(graph[v].events, graph[v].state,
                             graph[e].type, graph[t].state,
                             out_degree(v, graph)/2);
   }
   
   // calculate difference between new and old sum of rates
   double diff = tempSum - graph[v].rateSum;
   // set rateSum to new value
   graph[v].rateSum = tempSum;

   return diff;
}

#endif
