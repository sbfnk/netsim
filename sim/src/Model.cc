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
   gamma[Uninformed] = 1.;
   gamma[Informed] = 1.;
   delta[Uninformed] = 1.;
   delta[Informed] = 1.;
   beta[Uninformed][Uninformed] = 1.;
   beta[Uninformed][Informed] = 1.;
   beta[Informed][Uninformed] = 1.;
   beta[Informed][Informed] = 1.;
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
   beta[Uninformed][Uninformed]=read_dbl_val(ifile);
   
   // gamma_d 
   gamma[Uninformed]=read_dbl_val(ifile);
  
   // delta_d 
   delta[Uninformed]=read_dbl_val(ifile);
      
   // beta_--
   beta[Informed][Informed]=read_dbl_val(ifile);
   
   // gamma_i 
   gamma[Informed]=read_dbl_val(ifile);
   
   // delta_i 
   delta[Informed]=read_dbl_val(ifile);
   
   // beta_m1 
   beta[Uninformed][Informed]=read_dbl_val(ifile);
   
   // beta_m2 
   beta[Informed][Uninformed]=read_dbl_val(ifile);
   
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
                            VertexState state, EdgeType edge,
                            VertexState nbState) const
{
   double rateSum(.0);
   if (edge.getType() == Disease) {
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
   } else if (edge.getType() == Information) {
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

std::vector<VertexState> Model::getPossibleStates()
{
   std::vector<VertexState> v;
   v.push_back(VertexState(Susceptible, Uninformed));
   v.push_back(VertexState(Infected, Uninformed));
   v.push_back(VertexState(Recovered, Uninformed));
   v.push_back(VertexState(Susceptible, Informed));
   v.push_back(VertexState(Infected, Informed));
   v.push_back(VertexState(Recovered, Informed));
   return v;
}

std::vector<EdgeType> Model::getPossibleEdgeTypes()
{
   std::vector<EdgeType> e;
   e.push_back(Disease);
   e.push_back(Information);
   return e;
}

/******************************************************************/
// EdgeType::print
// for easier printout of the edge type
/******************************************************************/
void EdgeType::print(std::ostream& os) const
{
   switch (edgeType) {
       case Disease:
        os << "d";
        break;
       case Information:
        os << "i";
        break;
   }
}

/******************************************************************/
// VertexState::set
// set VertexState to state defined by string (S,s,I,i,R,r)
/******************************************************************/
void VertexState::set(std::string s)
{
   if (s == "S-") {
      setDisease(Susceptible);
      setInfo(Uninformed);
   } else if (s == "S+") {
      setDisease(Susceptible);
      setInfo(Informed);
   } else if (s == "I-") {
      setDisease(Infected);
      setInfo(Uninformed);
   } else if (s == "I+") {
      setDisease(Infected);
      setInfo(Informed);
   } else if (s == "R-") {
      setDisease(Recovered);
      setInfo(Uninformed);
   } else if (s == "R+") {
      setDisease(Recovered);
      setInfo(Informed);
   }
}

/******************************************************************/
// VertexState::getString
// convert state to string
/******************************************************************/
std::string VertexState::getString() const
{
   if (info == Informed) {
      switch (disease) {
       case Susceptible:
        return "S+";
       case Infected:
        return "I+";
       case Recovered:
        return "R+";
      }
   } else {
      switch (disease) {
       case Susceptible:
        return "S-";
       case Infected:
        return "I-";
       case Recovered:
        return "R-";
      }
   }
   return "";
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
        os << "\033[01;34m" << "S+" << "\033[0m";
        break;
       case Infected:
        os << "\033[01;31m" << "I+" << "\033[0m";
        break;
       case Recovered:
        os << "\033[01;32m" << "R+" << "\033[0m";
        break;
      }
   } else {
      switch (disease) {
       case Susceptible:
        os << "\033[00;34m" << "S-" << "\033[0m";
        break;
       case Infected:
        os << "\033[00;31m" << "I-" << "\033[0m";
        break;
       case Recovered:
        os << "\033[00;32m" << "R-" << "\033[0m";
        break;
      }
   }      
}

/******************************************************************/
// overloading of << operator for EdgeType so that it
// can be used in std output pipes
/******************************************************************/
std::ostream& operator<<(std::ostream& os, const EdgeType& e)
{
   e.print(os);
   return os;
}

/******************************************************************/
// overloading of << operator for VertexState so that it
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

