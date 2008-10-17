/*! \file Model.cc
  \brief Implementation of the Model class
*/
#include "Model.hh"

/*! \brief Stream operator for Label.

Streams a letter associated VertexState/EdgeType, useful for easy printout
\param[in, out] os The stream to write the Label to
\param[in] l The Label to stream

\return A reference to the stream written to
*/
std::ostream& operator<<(std::ostream& os, const Label& l)
{ l.print(os); return os; }

