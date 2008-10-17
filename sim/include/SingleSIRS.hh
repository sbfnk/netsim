/*! \file SingleSIRS.hh
  \brief The Models::SingleSIRS class.
*/
#ifndef SINGLESIRS_HH
#define SINGLESIRS_HH

#include "EpiModel.hh"

namespace Models {

  /*! \brief Class implementing a simple SIRS model spreading over multiple edges

  This defines 3 vertex states (S,I,R) and two edge types,
  with infection happening over any of the edge types
  
  */
  template <class Graph>
  class SingleSIRS :
    public EpiModel<State, Graph>
  {

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    
  public:

    SingleSIRS(unsigned int v = 0);
    ~SingleSIRS() {;}

    virtual Model<Graph>* clone() { return new SingleSIRS<Graph>(*this); }

    unsigned int getNodeEvents(eventList& events, State* currentState,
                               unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State* currentState,
                               unsigned int edge, State* currentNbState, 
                               unsigned int nb) const;

    bool isInfection(State before_state, State after_state) const
    { return ((before_state.getState() == Susceptible) &&
              (after_state.getState() == Infected)); }
    bool isRecovery(State before_state, State after_state) const
    { return ((before_state.getState() == Infected) &&
              (after_state.getState() != Infected)); }

    std::vector<StatRecorder<Graph>*> 
    getStatRecorders(const po::variables_map& vm) const
    { }

  private:

    unsigned int gamma; //!< Recovery rates.
    unsigned int delta; //!< Loss of immunity rates.
    unsigned int beta; //!< Infection rates.
  
  };

}

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
template <class Graph>
Models::SingleSIRS<Graph>::SingleSIRS(unsigned int v)
  : EpiModel<State, Graph>("SingleSIRS", v)
{
  /*************************************/
  // define vertex classes
  /************************************/
//   // susceptible
//   vertexStates.push_back(Label("S","01;32", 0, "", Label::rgbColour(0,0,255)));
//   // infected 
//   vertexStates.push_back(Label("I","01;31", 1, "", Label::rgbColour(255,0,0)));
//   // recovered
//   vertexStates.push_back(Label("R","01;34", 2, "", Label::rgbColour(0,255,0)));

  /*************************************/
  // define edge types
  /************************************/
  this->edgeTypes.push_back(Label("0", "", 0, "style=\"solid\""));
  this->edgeTypes.push_back(Label("1", "", 1, "style=\"dashed\""));

  /*************************************/
  // define model parameters
  /************************************/
  this->model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate")
    ("gamma", po::value<double>(),
     "recovery rate")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ;

  /*************************************/
  // assign model parameters to variables
  /************************************/
  this->rates.insert(std::make_pair("beta", &beta));
  this->rates.insert(std::make_pair("gamma", &gamma));
  this->rates.insert(std::make_pair("delta", &delta));
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::SingleSIRS<Graph>::getNodeEvents(eventList& events,
                                               State* currentState,
                                               unsigned int nb) const
{
   State* state = currentState;

   unsigned int rateSum(0);

   if (state->getState() == Recovered) {
     // loss of immunity
      Event immunityLoss;
      immunityLoss.rate = delta;
      immunityLoss.newState = 
	this->newState(State(Susceptible));
      immunityLoss.nb = nb;
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding loss of immunity event with rate " 
	            << immunityLoss.rate/1e+4 << std::endl;
        }
      }
   } else if (state->getState() == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma;
      recovery.newState = this->newState(State(Recovered));
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate/1e+4
                    << std::endl;
        }
      }
   }
   return rateSum;
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::SingleSIRS<Graph>::getEdgeEvents(eventList& events,
                                               State* currentState,
                                               unsigned int edge,
                                               State* currentNbState,
                                               unsigned int nb) const
{
   State* state = currentState;
   State* nbState = currentNbState;

   unsigned int rateSum(0);
   // infection
   if (state->getState() == Susceptible &&
       nbState->getState() == Infected) {
     Event infection;
     infection.rate = beta;
     infection.newState = this->newState(State(Infected));
     infection.nb = nb;
     infection.et = edge;
     if (infection.rate > 0) {
       events.push_back(infection);
       rateSum += infection.rate;
       if (this->getVerbose() >= 2) {
         std::cout << "Adding infection event with rate "
                   << infection.rate/1e+4 << std::endl;
       }
     }
   }

   return rateSum;
}

#endif
