/*! \file InfoSIRS.hh
  \brief The Models::InfoSIRS class.
*/
#ifndef INFOSIRS_HH
#define INFOSIRS_HH

#include "Model.hh"

namespace Models {

  /*! \brief Class implementing an SIRS model with information.

  This defines 6 vertex states (S-,I-,R-,S+,I+,R+) and 2 edge types (d and i),
  with information-dependent infection, recovery and loss of immunity, as well as
  information transmission and forgetting, local information generation and
  information generation over i links.
  
  */
  class InfoSIRS :
    public Model
  {

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible information states.
    enum infoStatesEnum {Uninformed, Informed};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};
    
  public:

    InfoSIRS(unsigned int v = 0);
    ~InfoSIRS() {;}

    void Init(po::variables_map& vm);
  
    unsigned int getNodeEvents(eventList& events, State state,
                         unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State state,
                         unsigned int edge, State nbState, 
                         unsigned int nb) const;

    bool isInfection(State before_state, State after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }
    bool isRecovery(State before_state, State after_state) const
    { return ((getDisease(before_state) == Infected) &&
              (getDisease(after_state) != Infected)); }
    bool isInformation(State before_state, State after_state) const
    { return ((getInfo(before_state) == Uninformed) &&
              (getInfo(after_state) == Informed)); }
    bool isForgetting(State  before_state, State after_state) const
    { return ((getInfo(before_state) == Informed) &&
              (getInfo(after_state) == Uninformed)); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State state) const
    { return (state.base % 3); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State state) const
    { return (state.base / 3); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(State dState, State iState) const
    { return dState.base+iState.base*3; }
    
  
  private:

    unsigned int gamma[2]; //!< Recovery rates.
    unsigned int delta[2]; //!< Loss of immunity rates.
    unsigned int beta[2][2]; //!< Infection rates.
    unsigned int alpha; //!< Information transmission rate.
    unsigned int lambda; //!< Rate of forgetting.
    unsigned int omega; //!< Local infromation generation rate.
    unsigned int nu; //!< Information generation rate over i-edges.
    double sigma; //!< Ratio between informed and uninformed susceptibility. 
  
  };

}

#endif
