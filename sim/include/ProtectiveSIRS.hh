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
#ifndef PROTECTIVESIRS_HH
#define PROTECTIVESIRS_HH

#include "Model.hh"

class ProtectiveSIRS :
  public Model
{
  
  enum diseaseStatesEnum {Susceptible,Infected,Recovered};
  enum infoStatesEnum {Uninformed, Informed};
  enum edgeTypesEnum {Disease, Information};

public:

  ProtectiveSIRS();
  ~ProtectiveSIRS();
  
  double getNodeEvents(eventList& events, unsigned int state,
                       unsigned int nb) const;
  double getEdgeEvents(eventList& events, unsigned int state,
                       unsigned int edge, unsigned int nbState,
                       unsigned int nb) const;

  bool isInfection(unsigned int before_state, unsigned int after_state) const
  { return ((getDisease(before_state) == Susceptible) &&
            (getDisease(after_state) == Infected)); }
  
  unsigned int getDisease(unsigned int state) const
  { return (state < 2 ? 0 : state - 1); }
  unsigned int getInfo(unsigned int state) const
  { return (state == 1); }
  unsigned int getState(unsigned int dState, unsigned int iState) const
  { return (dState == 0 ? iState : dState + 1); }
    
  
private:

  double gamma, delta; // model parameters
  double beta, nu, lambda, sigma; // model parameters
  
};

#endif
