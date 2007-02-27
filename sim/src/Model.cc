#include "Model.hh"

double Model::assignParam(po::variables_map& vm, std::string p)
{
  if (vm.count(p)) {
    return vm[p].as<double>();
  } else {
    std::cerr << "WARNING: no " << p << " given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    return 0;
  }
}

void Model::Init(po::variables_map& vm)
{
  for (std::map<std::string, double*>::iterator it = params.begin();
       it != params.end(); it++) {
    *((*it).second) = assignParam(vm, (*it).first);
  }
}

std::ostream& operator<<(std::ostream& os, const Label& l)
{ l.print(os); return os; }
