/*! \file SingleSIRS.hh
  \brief The Models::SingleSIRS class.
*/
#ifndef SINGLESIRS_HH
#define SINGLESIRS_HH

#include "Model.hh"

namespace Models {

  /*! \brief Class implementing a simple SIRS model spreading over multiple edges

  This defines 3 vertex states (S,I,R) and two edge types,
  with infection happening over any of the edge types
  
  */
  class SingleSIRS :
    public Model
  {

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    
  public:

    SingleSIRS(unsigned int v = 0);
    ~SingleSIRS() {;}

    double getNodeEvents(eventList& events, State state,
                         unsigned int nb) const;
    double getEdgeEvents(eventList& events, State state,
                         unsigned int edge, State nbState, 
                         unsigned int nb) const;

    bool isInfection(State before_state, State after_state) const
    { return ((before_state.base == Susceptible) &&
              (after_state.base == Infected)); }
    bool isRecovery(State before_state, State after_state) const
    { return ((before_state.base == Infected) &&
              (after_state.base != Infected)); }

  private:

    double gamma; //!< Recovery rates.
    double delta; //!< Loss of immunity rates.
    double beta; //!< Infection rates.
  
  };

}

#endif
