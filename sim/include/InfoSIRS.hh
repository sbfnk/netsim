/******************************************************************/
// Model.hh
// 
// this contains all classes related to the model:
//   -- VertexState
//        for the state a vertex can assume
//   -- event
//        for the events that can happen in the simulation
//   -- Model
//        for the model parameters and the events associated with it
/******************************************************************/
#ifndef INFOSIRS_HH
#define INFOSIRS_HH

#include "Model.hh"

class InfoSIRS :
  public Model
{

  enum diseaseStatesEnum {Susceptible,Infected,Recovered};
  enum infoStatesEnum {Uninformed, Informed};
  enum edgeTypesEnum {Disease, Information};
  
public:

  InfoSIRS(unsigned int v = 0);
  ~InfoSIRS();

  void Init(po::variables_map& vm);
  
  double getNodeEvents(eventList& events, unsigned int state) const;
  double getEdgeEvents(eventList& events, unsigned int state,
                       unsigned int edge, unsigned int nbState) const;

  bool isInfection(unsigned int before_state, unsigned int after_state) const
  { return ((getDisease(before_state) == Susceptible) &&
            (getDisease(after_state) == Infected)); }

  unsigned int getDisease(unsigned int state) const
  { return (state % 3); }
  unsigned int getInfo(unsigned int state) const
  { return (state / 3); }
  unsigned int getState(unsigned int dState, unsigned int iState) const
  { return dState+iState*3; }
    
  
private:

  double gamma[2], delta[2]; // model parameters
  double beta[2][2], alpha, nu, lambda, omega; // model parameters
  double sigma;
  
};

#endif
