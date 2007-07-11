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
  
    double getNodeEvents(eventList& events, unsigned int state,
                         unsigned int nb) const;
    double getEdgeEvents(eventList& events, unsigned int state,
                         unsigned int edge, unsigned int nbState,
                         unsigned int nb) const;

    bool isInfection(unsigned int before_state, unsigned int after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }
  
    //! Get the disease part of a full vertex state.
    unsigned int getDisease(unsigned int state) const
    { return (state < 2 ? 0 : state - 1); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(unsigned int state) const
    { return (state == 1); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(unsigned int dState, unsigned int iState) const
    { return (dState == 0 ? iState : dState + 1); }
    
  
  private:

    double gamma; //!< Recovery rate.
    double delta; //!< Loss of immunity rate.
    double beta; //!< Infection rate.
    double nu; //!< Information generateion rate over i-edges.
    double lambda; //!< Rate of forgetting.
    double sigma; //!< Ratio between informed and uninformed susceptibility.
  
  };
  
}

#endif
