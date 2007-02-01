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

/******************************************************************/
// Model constructor
/******************************************************************/
Model::Model()
{
  possibleStates.push_back(VertexState(Susceptible, Uninformed));
  possibleStates.push_back(VertexState(Infected, Uninformed));
  possibleStates.push_back(VertexState(Recovered, Uninformed));
  possibleStates.push_back(VertexState(Susceptible, Informed));
  possibleStates.push_back(VertexState(Infected, Informed));
  possibleStates.push_back(VertexState(Recovered, Informed));

  possibleEdgeTypes.push_back(Disease);
  possibleEdgeTypes.push_back(Information);
}

/******************************************************************/
// Vertex destructor
/******************************************************************/
Model::~Model()
{
}

/******************************************************************/
// Model::InitFromFile
// initialize parameters from file
/******************************************************************/
void Model::Init(po::variables_map& vm)
{
  if (vm.count("beta--")) {
    beta[0][0]=vm["beta--"].as<double>();
  } else {
    std::cerr << "WARNING: no beta-- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    beta[0][0]=0;
  }
  if (vm.count("beta+-")) {
    beta[1][0]=vm["beta+-"].as<double>();
  } else {
    std::cerr << "WARNING: no beta+- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    beta[1][0]=0;
  }
  if (vm.count("beta-+")) {
    beta[0][1]=vm["beta-+"].as<double>();
  } else {
    std::cerr << "WARNING: no beta-+ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    beta[0][1]=0;
  }
  if (vm.count("beta++")) {
    beta[1][1]=vm["beta++"].as<double>();
  } else {
    std::cerr << "WARNING: no beta++ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    beta[1][1]=0;
  }
  if (vm.count("gamma-")) {
    gamma[0]=vm["gamma-"].as<double>();
  } else {
    std::cerr << "WARNING: no gamma- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    gamma[0]=0;
  }
  if (vm.count("gamma+")) {
    gamma[1]=vm["gamma+"].as<double>();
  } else {
    std::cerr << "WARNING: no gamma+ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    gamma[1]=0;
  }
  if (vm.count("delta-")) {
    delta[0]=vm["delta-"].as<double>();
  } else {
    std::cerr << "WARNING: no delta- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    delta[0]=0;
  }
  if (vm.count("delta+")) {
    delta[1]=vm["delta+"].as<double>();
  } else {
    std::cerr << "WARNING: no delta+ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    delta[1]=0;
  }
  if (vm.count("alpha")) {
    alpha=vm["alpha"].as<double>();
  } else {
    std::cerr << "WARNING: no alpha given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    alpha=0;
  }
  if (vm.count("nu")) {
    nu=vm["nu"].as<double>();
  } else {
    std::cerr << "WARNING: no nu given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    nu=0;
  }
  if (vm.count("omega")) {
    omega=vm["omega"].as<double>();
  } else {
    std::cerr << "WARNING: no omega given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    omega=0;
  }
  if (vm.count("lambda")) {
    lambda=vm["lambda"].as<double>();
  } else {
    std::cerr << "WARNING: no lambda given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    lambda=0;
  }
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
   if (state.getDisease() == Infected) {
      // recovery
      event recovery;
      recovery.rate = gamma[state.getInfo()];
      recovery.newState.setDisease(Recovered);
      recovery.newState.setInfo(state.getInfo());
      events.push_back(recovery);
      rateSum += recovery.rate;
      if (state.getInfo() == Uninformed) {
        // local information generation
        event localInfo;
        localInfo.rate = omega;
        localInfo.newState.setDisease(state.getDisease());
        localInfo.newState.setInfo(Informed);
        events.push_back(localInfo);
        rateSum += localInfo.rate;
      }
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

const std::vector<VertexState>& Model::getPossibleStates() const
{
  return possibleStates;
}

const std::vector<EdgeType>& Model::getPossibleEdgeTypes() const
{
  return possibleEdgeTypes;
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
