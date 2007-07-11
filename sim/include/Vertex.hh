/*! \file Vertex.hh
  \brief Contains the Edge and Vertex classes
*/

#ifndef VERTEX_HH
#define VERTEX_HH

#include <list>

/*! \brief Container for edge properties.
\ingroup gillespie_simulator
*/
class Edge
{      
public:

  //! Constructor.
  Edge(): type(0), parallel(false) {;}

  /*! \brief Constructor.
    \param[in] eType type initialiser.
  */
  Edge(int eType): type(eType), parallel(false) {;}

  /*! \brief Constructor.
    \param[in] eType type initialiser.
    \param[in] p parallel initialiser.
  */
  Edge(int eType, bool p): type(eType), parallel(false) {;}
      
  //! A number corresponding to an edge type, to be defined by the used Model.
  unsigned int type; 
  bool parallel; //!< Whether the edge has a parallel counterpart or not
      
};

/*! \brief Container for vertex properties.
\ingroup gillespie_simulator
*/
class Vertex
{      
public:
      
  //! Constructor.
  Vertex() : state(0) {;}
  
  /*! \brief Constructor.
    \param[in] s state initialiser.
  */
  Vertex(int s)
    : state(s) {;}

  //! Destructor
  virtual ~Vertex() {;}
      
  /*! A number corresponding to the current state of the vertex, to be defined
    by the used Model.
  */
  unsigned int state; 
  //! The sum of rates of all events that can change the state of the vertex.
  double rateSum;
      
  eventList events; //!< A list of events associated with vertex.
      
};

#endif
