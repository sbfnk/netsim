#include <iomanip>

#include "StatRecorder.hh"

std::string generateFileName(std::string nameBase, unsigned int id,
                             std::string ext)
{
  std::stringstream s;
  s << nameBase << std::setw(6) << std::setfill('0') << id;
  if (ext.size() > 0) s << "." << ext;
  return s.str();
}

