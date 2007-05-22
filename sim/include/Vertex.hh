/******************************************************************/
// Vertex.hh
// contains the Edge and Vertex classes, to be used by the different
// graph implementations
/******************************************************************/

#ifndef VERTEX_HH
#define VERTEX_HH

#include <list>

class Edge
{      
public:
      
  Edge(): type(0), parallel(false), rewired(false) {;}
  Edge(int eType): type(eType), parallel(false), rewired(false) {;}
  Edge(int eType, bool p): type(eType), parallel(false), rewired(p) {;}
      
  unsigned int type;
  bool parallel;
  bool rewired;
      
};

class Vertex
{      
public:
      
  Vertex() : state(0) {;}
  Vertex(int s)
    : state(s) {;}
      
  virtual ~Vertex() {;}
      
  unsigned int state; // present (disease and information) state
  double rateSum; // sum of the rates of all events that can
  // affect the vertex (change its state)
      
  eventList events; // list of events associated with vertex
      
};

#endif
