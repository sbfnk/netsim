/*! \file Vertex.hh
  \brief Contains the Vertex class
*/

#ifndef VERTEX_HH
#define VERTEX_HH

#include <list>

/*! \brief Container for vertex properties.
\ingroup gillespie_simulator
*/
class Vertex
{      
public:
      
  //! Constructor.
  Vertex() : state(0), infected(0) {;}
  
  /*! \brief Constructor.
    \param[in] s state initialiser.
  */
  Vertex(int s)
    : state(s), infected(0) {;}

  //! Destructor
  virtual ~Vertex() {;}
      
  /*! A number corresponding to the current state of the vertex, to be defined
    by the used Model.
  */
  unsigned int state;
  //! The sum of rates of all events that can change the state of the vertex.
  double rateSum;
      
  eventList events; //!< A list of events associated with vertex.

  bool infected; //!< Whether the vertex got infected in its lifetime
      
};

#endif
