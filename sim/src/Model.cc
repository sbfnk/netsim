/*******************************************************************/
//
// Class Model
// --------------------
//
// This is the class implementing the used model, i.e. the parameters
// and the dependance of the rates of various processes on node states
// and edge properties
//
/******************************************************************/

#include <iostream>
#include <fstream>

#include "Model.hh"

double read_dbl_val(std::ifstream& ifile);

/******************************************************************/
// Model constructor
/******************************************************************/
Model::Model()
{
}

/******************************************************************/
// Vertex destructor
/******************************************************************/
Model::~Model()
{
}

/******************************************************************/
// Model::InitDefaultParams
// initialize parameters to some very unimaginative default values
/******************************************************************/
void Model::InitDefaultParams()
{
   gamma[0] = 1.;
   gamma[1] = 1.;
   delta[0] = 1.;
   delta[1] = 1.;
   beta[0][0] = 1.;
   beta[0][1] = 1.;
   beta[1][0] = 1.;
   beta[1][1] = 1.;
   alpha = 1.;
   nu = 1.;
   lambda = 1.;
}

/******************************************************************/
// Model::InitFromFile
// initialize parameters from file
/******************************************************************/
int Model::InitFromFile(std::string fileName)
{
   std::string line;
   std::ifstream ifile;
   
   ifile.open (fileName.c_str(), std::ios::in); 
   if(ifile.fail())
   {
      return 1;
   }

   // reading Model content of init.dat   
   std::cout << "Reading model parameters"; 
   
   int i=-1;
   while ((i<0) && !(ifile.eof()))
   {
      getline(ifile,line); // reading a line
      i=line.find("Model parameters"); // looking for "Model ... %%%"
   }
   
   // beta_++
   beta[0][0]=read_dbl_val(ifile);
   
   // gamma_d 
   gamma[0]=read_dbl_val(ifile);
  
   // delta_d 
   delta[0]=read_dbl_val(ifile);
      
   // beta_--
   beta[1][1]=read_dbl_val(ifile);
   
   // gamma_i 
   gamma[1]=read_dbl_val(ifile);
   
   // delta_i 
   delta[1]=read_dbl_val(ifile);
   
   // beta_m1 
   beta[0][1]=read_dbl_val(ifile);
   
   // beta_m2 
   beta[1][0]=read_dbl_val(ifile);
   
   // alpha 
   alpha=read_dbl_val(ifile);
   
   // nu 
   nu=read_dbl_val(ifile);
   
   // lambda 
   lambda=read_dbl_val(ifile);
   
   // Qd 
   //Qd=read_dbl_val(ifile);
   
   // Qi 
   //Qi=read_dbl_val(ifile);
   
   // N
//    N=read_dbl_val(ifile);
   
   // njac 
//    njac=obj.GetNvars();

   ifile.close(); 
   if(!ifile.is_open())
      std::cout << " ... done\n";
   return 0;
}

/******************************************************************/
// Model::getNodeEvents
// get the events that can happen for a given state of a node. Stores the
// events in the events list and returns the sum of their rates
/******************************************************************/
double Model::getNodeEvents(eventList& events,
                            VertexState state) const
{
   double rateSum(.0);

   // loss of immunity
   if (state.getDisease() == Recovered) {
      event immunityLoss;
      immunityLoss.rate = delta[state.getInfo()];
      immunityLoss.newState.setDisease(Susceptible);
      immunityLoss.newState.setInfo(state.getInfo());
      events.push_back(immunityLoss);
      rateSum += immunityLoss.rate;
   } else
   // recovery
   if (state.getDisease() == Infected) {
      event recovery;
      recovery.rate = gamma[state.getInfo()];
      recovery.newState.setDisease(Recovered);
      recovery.newState.setInfo(state.getInfo());
      events.push_back(recovery);
      rateSum += recovery.rate;
   }
   // information loss
   if (state.getInfo() == Informed) {
      event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState.setDisease(state.getDisease());
      infoLoss.newState.setInfo(Uninformed);
      events.push_back(infoLoss);
      rateSum += infoLoss.rate;
   }

   return rateSum;
}

/******************************************************************/
// Model::getNodeEvents
// get the events that can happen for a given edge type and the states of the
// neighbouring nodes. Stores the events in the events list and returns the sum
// of their rates 
/******************************************************************/
double Model::getEdgeEvents(eventList& events,
                            VertexState state, int edge,
                            VertexState nbState) const
{
   double rateSum(.0);
   if (edge == Disease) {
      // infection
      if (state.getDisease() == Susceptible &&
          nbState.getDisease() == Infected) {
         event infection;
         infection.rate = beta[state.getInfo()][nbState.getInfo()];
         infection.newState.setInfo(state.getInfo());
         infection.newState.setDisease(Infected);
         events.push_back(infection);
         rateSum += infection.rate;
      }
   } else if (edge == Information) {
      // information transmission
      if (state.getInfo() == Uninformed && nbState.getInfo() == Informed) {
         event infoTransmission;
         infoTransmission.rate = alpha;
         infoTransmission.newState.setDisease(state.getDisease());
         infoTransmission.newState.setInfo(Informed);
         events.push_back(infoTransmission);
         rateSum += infoTransmission.rate;
      }
      // information generation
      if (state.getInfo() == Uninformed && nbState.getDisease() == Infected) {
         event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState.setDisease(state.getDisease());
         infoGeneration.newState.setInfo(Informed);
         events.push_back(infoGeneration);
         rateSum += infoGeneration.rate;
      }
   }

   return rateSum;
}

/******************************************************************/
// VertexState::print
// for easier printout of the vertex state
/******************************************************************/
void VertexState::print(std::ostream& os) const
{
   if (info == Informed) {
      switch (disease) {
       case Susceptible:
        os << "\033[01;34m" << "s" << "\033[0m";
        break;
       case Infected:
        os << "\033[01;31m" << "i" << "\033[0m";
        break;
       case Recovered:
        os << "\033[01;32m" << "r" << "\033[0m";
        break;
      }
   } else {
      switch (disease) {
       case Susceptible:
        os << "\033[00;34m" << "S" << "\033[0m";
        break;
       case Infected:
        os << "\033[00;31m" << "I" << "\033[0m";
        break;
       case Recovered:
        os << "\033[00;32m" << "R" << "\033[0m";
        break;
      }
   }      
}

/******************************************************************/
// overloading of << operator for a VertexState so that it
// can be used in std output pipes
/******************************************************************/
std::ostream& operator<<(std::ostream& os, const VertexState& v)
{
   v.print(os);
   return os;
}


/******************************************************************/
// read_dbl_val from ode_io_utils.cc, needed for InitFromFile
/******************************************************************/
double read_dbl_val(std::ifstream& ifile)
{
   std::string line, str_tmp;
   
   getline(ifile,line); 
   int i=line.find("="); // looking for the =  
   if(i>=0) // check if = was found 
   {
      str_tmp=line.substr(i+2,line.length()-i); // extracting substring
   }
   else
   {
      std::cout << "... problem reading params, no = was found\n";
      exit(1);
   }

   return atof(str_tmp.c_str());

}

