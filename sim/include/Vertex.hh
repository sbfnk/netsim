/*! \file Vertex.hh
  \brief Contains the Vertex class
*/

#ifndef VERTEX_HH
#define VERTEX_HH

#include <list>

#include <Model.hh>

/*! \brief Container for vertex properties.
\ingroup gillespie_simulator
*/
class Vertex
{      
public:

  //! Constructor.
  Vertex()
    : state() {;}
  
  /*! \brief Constructor.
    \param[in] s state initialiser.
  */
  Vertex(unsigned int s)
    : state(s) {;}

  //! Destructor
  virtual ~Vertex() {;}

  State state;
  
  //! The sum of rates of all events that can change the state of the vertex.
  double rateSum;
      
  eventList events; //!< A list of events associated with vertex.
};

#endif
