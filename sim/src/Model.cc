/*! \file Model.cc
  \brief Implementation of the Model class
*/
#include "Model.hh"

/*! \brief Initialise model paramters.
   
This should be called after the command line parameters of the model haven been
assigned. It initialises the model parameter variables with the values found in
the command line parameters.

\param[in] vm The map of command line parameters
*/

void Model::Init(po::variables_map& vm)
{
  // loop over all model parameters
  for (std::map<std::string, double*>::iterator it = params.begin();
       it != params.end(); it++) {
    if (vm.count(it->first)) {
      // command line parameter has been specified, assign to model variable
      *(it->second) = vm[it->first].as<double>();
    } else {
      std::cerr << "WARNING: no " << it->first << " given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      *(it->second) = 0;
    }
  }
}

//! Print the model parameters.
void Model::Print()
{
  std::cout << "Model parameters:" << std::endl;
  std::cout << "=================" << std::endl;
  for (std::map<std::string, double*>::iterator it = params.begin();
       it != params.end(); it++) {
    std::cout << it->first << ": " << *(it->second) << std::endl; 
  }
}

/*! \brief Stream operator for Label.

Streams a letter associated VertexState/EdgeType, useful for easy printout
\param[in, out] os The stream to write the Label to
\param[in] l The Label to stream

\return A reference to the stream written to
*/
std::ostream& operator<<(std::ostream& os, const Label& l)
{ l.print(os); return os; }
