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

#endif
