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
  Vertex()
    : state(0), state_detail(0.) {;}
  
  /*! \brief Constructor.
    \param[in] s state initialiser.
  */
  Vertex(unsigned int s)
    : state(s), state_detail(0.) {;}

  //! Destructor
  virtual ~Vertex() {;}
      
  /*! Integer variable corresponding to the current state of the vertex, to be defined
    by the used Model.
  */
  unsigned int state;
  /*! Real-valued variable refining the current state of the vertex, to be defined
    by the used Model.
  */
  double state_detail;
  
  //! The sum of rates of all events that can change the state of the vertex.
  double rateSum;
      
  eventList events; //!< A list of events associated with vertex.
};

#endif
