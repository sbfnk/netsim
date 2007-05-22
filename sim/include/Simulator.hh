#ifndef SIMULATOR_HH
#define SIMULATOR_HH

#include "Model.hh"

class Simulator
{

public:
      
  Simulator(const Model& m, unsigned int v = 0) : model(m), verbose(v), time(0.) {}
  virtual ~Simulator() {;}
      
  virtual void initialize() {;}
  
  virtual bool updateState() = 0;

  virtual void print() {;}

  double getTime() const { return time; };
  void updateTime(double t) { time += t; };

protected:
  
  const Model& model;
  unsigned int verbose;

private:

  double time;
      
};


/******************************************************************/
// generateEventList function
// generates an event list for a given vertex in a given graph,
// using a given model. The vertices of the event must be of type
// Vertex, otherwise this will throw a compile error
/******************************************************************/
template <class Graph, class Model>
double generateEventList(Graph& graph,
                         typename boost::graph_traits<Graph>::vertex_descriptor v,
                         const Model& model, unsigned int verbose = 0)
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

   if (verbose >= 2) {
     std::cout << "Generating events list for vertex #" << v << " ("
               << model.getVertexStates()[graph[v].state] << ")" << std::endl;
   }

   // get node events
   tempSum += model.getNodeEvents(graph[v].events, graph[v].state);

   // get edge events
   out_edge_iterator oi, oi_end;;

   for (tie(oi, oi_end) = boost::out_edges(v, graph);
        oi != oi_end; ++oi) {
      edge_descriptor e = *oi;
      vertex_descriptor t =  target(e, graph);
      tempSum +=
         model.getEdgeEvents(graph[v].events, graph[v].state,
                             graph[e].type, graph[t].state);
   }
   
   // calculate difference between new and old sum of rates
   double diff = tempSum - graph[v].rateSum;
   // set rateSum to new value
   graph[v].rateSum = tempSum;


   return diff;
}

#endif
