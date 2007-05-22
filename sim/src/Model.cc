#include "Model.hh"

void Model::Init(po::variables_map& vm)
{
  for (std::map<std::string, double*>::iterator it = params.begin();
       it != params.end(); it++) {
    if (vm.count(it->first)) {
      *(it->second) = vm[it->first].as<double>();
    } else {
      std::cerr << "WARNING: no " << it->first << " given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      *(it->second) = 0;
    }
  }
}

void Model::Print()
{
  std::cout << "Model parameters:" << std::endl;
  std::cout << "=================" << std::endl;
  for (std::map<std::string, double*>::iterator it = params.begin();
       it != params.end(); it++) {
    std::cout << it->first << ": " << *(it->second) << std::endl; 
  }
}

std::ostream& operator<<(std::ostream& os, const Label& l)
{ l.print(os); return os; }
