#include "Vertex.hh"

std::ostream& operator<<(std::ostream& os, const State& s)
{ s.print(os); return os; }
