/*! \file Vertex.hh
  \brief Contains the Vertex class
*/

#ifndef VERTEX_HH
#define VERTEX_HH

#include <vector>
#include <iostream>

//----------------------------------------------------------
/*! \brief The state a vertex can assume
  
\ingroup models
\ingroup gillespie_simulator
*/
//! Type denoting the state of an individual
class State {
  
public:
  /*! \brief Constructor.
      \param[in] b base state.
  */
  State(int b = 0) : state(b) {;}
  virtual ~State() {;}

  virtual State* clone() { return new State(*this); }

  //! Print a coloured letter representing the vertex state.
  virtual void print(std::ostream &os) const
  { os << state; }

  unsigned int getState() const { return state; }
  void setState(int s) { state = s; }

private:

  int state;
  
};

std::ostream& operator<<(std::ostream& os, const State& s);

//----------------------------------------------------------
/*! \brief An event which can happen to a vertex.
  
\ingroup models
\ingroup gillespie_simulator
*/
class Event
{
public:
  /*! \brief Constructor.
  \param[in] r rate initialiser
  \param[in] s newState initialiser
  \param[in] n nb initialiser
  \param[in] e et initialiser
  */
  Event(unsigned int r = 0, State* s = 0, unsigned int n=0, unsigned int e=0) :
    rate(r), newState(s), nb(n), et(e)
  {;}
  
  ~Event() { delete newState; }
  Event(const Event& e)
  {
    newState = e.newState->clone();
    rate = e.rate; nb = e.nb; et = e.et;
  }

  Event& operator=(const Event& e)
  {
    if (this == &e) return *this;
    newState = e.newState->clone();
    rate = e.rate; nb = e.nb; et = e.et;
    return *this;
  }

  //! The rate at which an event occurs (depends on model parameters),
  //  multiplied by 10^5 for integer representation
  unsigned int rate; 
  //! The state an event will change a vertex to
  State* newState; 
  //! The neighbour "responsible" for the event
  unsigned int nb; 
  //! The edge type over which event is transmitted (if applicable)
  unsigned int et;

};

typedef std::vector<Event> eventList;

/*! \brief Container for vertex properties.
\ingroup gillespie_simulator
*/
class Vertex
{      
public:

  //! Constructor.
  Vertex()
    : state(new State(0)) {;}
  
  /*! \brief Constructor.
    \param[in] s state initialiser.
  */
  Vertex(State* s)
    : state(s) {;}

  Vertex(const Vertex& v)
  { state = v.state->clone(); }

  Vertex& operator=(const Vertex& v)
  { if (&v != this) {this->state = v.state->clone();} return *this; }
    
  //! Destructor
  virtual ~Vertex() { if (state) delete state;}

  State* state;
  
  //! The sum of rates of all events that can change the state of the vertex.
  unsigned int rateSum;
      
  eventList events; //!< A list of events associated with vertex.
};

#endif
