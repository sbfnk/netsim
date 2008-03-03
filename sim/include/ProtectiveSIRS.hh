/*! \file ProtectiveSIRS.hh
  \brief The Models::ProtectiveSIRS class.
*/
#ifndef PROTECTIVESIRS_HH
#define PROTECTIVESIRS_HH

#include "Model.hh"

namespace Models {

  /*! \brief Class implementing an SIRS model with individuals protecting
    themselves. 

  This defines 4 vertex states (S-,S+,I,R) and 2 edge types (d and i),
  with information-dependent susceptibility, as well as recovery, loss of
  immunity, information generation over i links and forgetting.
  
  */
  class ProtectiveSIRS :
    public Model
  {
  
    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible information states.
    enum infoStatesEnum {Uninformed, Informed};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};

  public:

    ProtectiveSIRS(unsigned int v = 0);
    ~ProtectiveSIRS() {;}
  
    unsigned int getNodeEvents(eventList& events, State state,
                               unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State state,
                               unsigned int edge, State nbState,
                               unsigned int nb) const;

    bool isInfection(State before_state, State after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }
  
    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State state) const
    { return (state.base < 2 ? 0 : state.base - 1); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State state) const
    { return (state.base == 1); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(State dState, State iState) const
    { return (dState.base == 0 ? iState.base : dState.base + 1); }
    
  
  private:

    unsigned int gamma; //!< Recovery rate.
    unsigned int delta; //!< Loss of immunity rate.
    unsigned int beta; //!< Infection rate.
    unsigned int nu; //!< Information generateion rate over i-edges.
    unsigned int lambda; //!< Rate of forgetting.
    double sigma; //!< Ratio between informed and uninformed susceptibility.
  
  };
  
}

#endif
