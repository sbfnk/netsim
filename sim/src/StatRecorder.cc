#include <iomanip>

#include "StatRecorder.hh"

std::string generateFileName(std::string nameBase, unsigned int id)
{
  std::stringstream s;
  s << nameBase << std::setw(6) << std::setfill('0') << id;
  return s.str();
}

