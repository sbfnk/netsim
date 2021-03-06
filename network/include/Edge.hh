/*! \file Edge.hh
  \brief Contains the Edge class
*/

#ifndef EDGE_HH
#define EDGE_HH

/*! \brief Container for edge properties.
\ingroup gillespie_simulator
*/
class Edge
{      
public:

  //! Constructor.
  Edge(): type(0), parallel(false), weight(1.) {;}

  /*! \brief Constructor.
    \param[in] eType type initialiser.
  */
  Edge(int eType): type(eType), parallel(false), weight(1.) {;}

  /*! \brief Constructor.
    \param[in] eType type initialiser.
    \param[in] p parallel initialiser.
  */
  Edge(int eType, bool p): type(eType), parallel(p), weight(1.) {;}
  Edge(int eType, double w): type(eType), parallel(false), weight(w) {;}
      
  //! A number corresponding to an edge type, to be defined by the used Model.
  unsigned int type; 
  bool parallel; //!< Whether the edge has a parallel counterpart or not

  double weight;

  std::size_t index;

  bool operator==(const Edge& rhs)
  { return type == rhs.type; }
      
};

#endif
