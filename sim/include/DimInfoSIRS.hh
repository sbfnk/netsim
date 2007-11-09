/*! \file DimInfoSIRS.hh
  \brief The Models::DimInfoSIRS class.
*/
#ifndef DIMINFOSIRS_HH
#define DIMINFOSIRS_HH

#include "Model.hh"

namespace Models {

  /*! \brief Class implementing an SIRS model with information.

  This defines 6 vertex states (S-,I-,R-,S+,I+,R+) and 2 edge types (d and i),
  with information-dependent infection, recovery and loss of immunity, as well as
  information transmission and forgetting, local information generation and
  information generation over i links. Information diminishes at each
  transmission by a factor rho.
  
  */
  class DimInfoSIRS :
    public Model
  {

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible information states.
    enum infoStatesEnum {Informed};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};
    
  public:

    DimInfoSIRS(unsigned int v = 0);
    ~DimInfoSIRS() {;}

    void Init(po::variables_map& vm);
  
    double getNodeEvents(eventList& events, State state,
                         unsigned int nb) const;
    double getEdgeEvents(eventList& events, State state,
                         unsigned int edge, State nbState,
                         unsigned int nb) const;

    bool isInfection(State before_state, State after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }
    bool isRecovery(State before_state, State after_state) const
    { return ((getDisease(before_state) == Infected) &&
              (getDisease(after_state) != Infected)); }
    bool isInformation(State before_state, State after_state) const
    { return (before_state.detail < after_state.detail); }
    bool isForgetting(State before_state, State after_state) const
    { return (before_state.detail > after_state.detail); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State state) const
    { return (state.base % 3); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State state) const
    { return 0; }
    
  private:

    double gamma; //!< Recovery rates.
    double delta; //!< Loss of immunity rates.
    double beta; //!< Infection rates.
    double alpha; //!< Information transmission rate.
    double lambda; //!< Rate of forgetting.
    double omega; //!< Local infromation generation rate.
    double nu; //!< Information generation rate over i-edges.
    double rho; //!< Ratio between informed and uninformed susceptibility.

    double threshold; //!< Threshold below which information will not be passed

  };

}

#endif
