/*! \file VaccinationSIRS.hh
  \brief The Models::VaccinationSIRS class.
*/
#ifndef VACCINATIONSIRS_HH
#define VACCINATIONSIRS_HH

#include "Model.hh"

namespace Models {

  /*! \brief Class implementing the SIRS model with information and vaccination.

  This defines 6 vertex states (S-,I-,R-,S+,I+,R+) and 2 edge types (d and i),
  with information-dependent infection, recovery and loss of immunity, as well as
  information transmission and forgetting, local information generation,
  information generation over i links and vaccination.
  
  */
  class VaccinationSIRS :
    public Model
  {

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible information states.
    enum infoStatesEnum {Uninformed, Informed};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};
  
  public:

    VaccinationSIRS(unsigned int v = 0);
    ~VaccinationSIRS() {;}

    void Init(po::variables_map& vm);
  
    double getNodeEvents(eventList& events, unsigned int state,
                         unsigned int nb) const;
    double getEdgeEvents(eventList& events, unsigned int state,
                         unsigned int edge, unsigned int nbState,
                         unsigned int nb) const;

    bool isInfection(unsigned int before_state, unsigned int after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }

    bool isRecovery(unsigned int before_state, unsigned int after_state) const
    { return ((getDisease(before_state) == Infected) &&
              (getDisease(after_state) != Infected)); }

    bool isInformation(unsigned int before_state, unsigned int after_state) const
    { return ((getInfo(before_state) == Uninformed) &&
              (getInfo(after_state) == Informed)); }

    bool isForgetting(unsigned int before_state, unsigned int after_state) const
    { return ((getInfo(before_state) == Informed) &&
              (getInfo(after_state) == Uninformed)); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(unsigned int state) const
    { return (state % 3); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(unsigned int state) const
    { return (state / 3); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(unsigned int dState, unsigned int iState) const
    { return dState+iState*3; }
    
  
  private:

    double gamma[2]; //!< Recovery rates.
    double delta[2]; //!< Loss of immunity rates.
    double beta[2][2]; //!< Infection rates.
    double alpha; //!< Information transmission rate.
    double lambda; //!< Rate of forgetting.
    double omega; //!< Local infromation generation rate.
    double nu; //!< Information generation rate over i-edges.
    double sigma; //!< Ratio between informed and uninformed susceptibility.
    double theta; //!< Vaccination rate of informed susceptibles.
  
  };
  
}

#endif
