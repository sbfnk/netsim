/*! \file ProtectiveSIRS.hh
  \brief The Models::ProtectiveSIRS class.
*/
#ifndef PROTECTIVESIRS_HH
#define PROTECTIVESIRS_HH

#include "EpiModel.hh"

namespace Models {

  /*! \brief Class implementing an SIRS model with individuals protecting
    themselves.

  This defines 4 vertex states (S-,S+,I,R) and 2 edge types (d and i),
  with information-dependent susceptibility, as well as recovery, loss of
  immunity, information generation over i links and forgetting.

  */
  template <class Graph>
  class ProtectiveSIRS :
    public EpiModel<State, Graph>
  {

    using EpiModel<State, Graph>::newState;

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible information states.
    enum infoStatesEnum {Uninformed, Informed};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};

  public:

    ProtectiveSIRS(unsigned int v = 0);
    ~ProtectiveSIRS() {;}

    virtual Model<Graph>* clone() const { return new ProtectiveSIRS<Graph>(*this); }

    virtual State* newState(unsigned int disease, unsigned int info) const
    { State* s = new State(info*3+disease); return s; }

    unsigned int getNodeEvents(eventList& events, State* currentState,
                               unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State* currentState,
                               unsigned int edge, State* currentNbState,
                               unsigned int nb) const;

    bool isInfection(State* before_state, State* after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State* state) const
    { return (state->getState() < 2 ? 0 : state->getState() - 1); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State* state) const
    { return (state->getState() == 1); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(State* dState, State* iState) const
    { return (dState->getState() == 0 ? iState->getState() : dState->getState() + 1); }

    std::vector<StatRecorder<Graph>*>
    getStatRecorders(const po::variables_map& vm) const
    {;}

  private:

    unsigned int gamma; //!< Recovery rate.
    unsigned int delta; //!< Loss of immunity rate.
    unsigned int beta; //!< Infection rate.
    unsigned int nu; //!< Information generateion rate over i-edges.
    unsigned int lambda; //!< Rate of forgetting.
    double sigma; //!< Ratio between informed and uninformed susceptibility.

  };

}

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
template <class Graph>
Models::ProtectiveSIRS<Graph>::ProtectiveSIRS(unsigned int v)
  : EpiModel<State, Graph>("ProtectiveSIRS", v)
{
//   // susceptible uninformed
//   vertexStates.push_back(Label("S-","00;32", 0, "fillcolor=\"royalblue4\""));
//   // susceptible informed
//   vertexStates.push_back(Label("S+","01;32", 1, "fillcolor=\"royalblue\""));
//   // infected
//   vertexStates.push_back(Label("I","00;31", 2, "fillcolor=\"red\""));
//   // recovered
//   vertexStates.push_back(Label("R","00;34", 3, "fillcolor=\"green\""));

  this->edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  this->edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

  this->model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate")
    ("sigma", po::value<double>(),
     "reduction of transmission rate due to information")
    ("gamma", po::value<double>(),
     "recovery rate")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ;

  this->rates.insert(std::make_pair("beta", &beta));
  this->rates.insert(std::make_pair("gamma", &gamma));
  this->rates.insert(std::make_pair("delta", &delta));
  this->rates.insert(std::make_pair("nu", &nu));
  this->rates.insert(std::make_pair("lambda", &lambda));
  this->params.insert(std::make_pair("sigma", &sigma));
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::ProtectiveSIRS<Graph>::getNodeEvents(eventList& events,
                                                   State* currentState,
                                                   unsigned int nb) const
{
   State* state = currentState;

   unsigned int rateSum(0);

   // loss of immunity
   if (getDisease(state) == Recovered) {
      Event immunityLoss;
      immunityLoss.rate = delta;
      immunityLoss.newState =
	this->newState(Susceptible, getInfo(state));
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
      }
   } else
   if (getDisease(state) == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma;
      recovery.newState =
	this->newState(Recovered, getInfo(state));
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
      }
   }
   // information loss
   if (getInfo(state) == Informed) {
      Event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState =
	this->newState(getDisease(state), Uninformed);
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
      }
   }

   return rateSum;
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::ProtectiveSIRS<Graph>::getEdgeEvents(eventList& events,
                                                   State* currentState,
                                                   unsigned int edge,
                                                   State* currentNbState,
                                                   unsigned int nb) const
{
   State* state = currentState;
   State* nbState = currentNbState;

   unsigned int rateSum(0);
   if (edge == Disease) {
      // infection
      if (getDisease(state) == Susceptible &&
          getDisease(nbState) == Infected) {
         Event infection;
         infection.rate = beta;
         if (getInfo(state) == 1) {
           infection.rate =
             static_cast<unsigned int>(infection.rate * sigma + .5);
         }
         infection.newState =
	   this->newState(Infected, getInfo(state));
         if (infection.rate > 0) {
           events.push_back(infection);
           rateSum += infection.rate;
         }
      }
   } else if (edge == Information) {
     if (getDisease(state) == Susceptible &&
         getInfo(state) == Uninformed &&
	 getDisease(nbState) == Infected) {
      // information generation
         Event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState =
	   this->newState(getDisease(state), Informed);
         if (infoGeneration.rate > 0) {
           events.push_back(infoGeneration);
           rateSum += infoGeneration.rate;
         }
      }
   }

   return rateSum;
}

#endif
