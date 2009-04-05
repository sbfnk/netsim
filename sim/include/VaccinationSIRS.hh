/*! \file VaccinationSIRS.hh
  \brief The Models::VaccinationSIRS class.
*/
#ifndef VACCINATIONSIRS_HH
#define VACCINATIONSIRS_HH

#include "EpiModel.hh"

namespace Models {

  /*! \brief Class implementing the SIRS model with information and vaccination.

  This defines 6 vertex states (S-,I-,R-,S+,I+,R+) and 2 edge types (d and i),
  with information-dependent infection, recovery and loss of immunity, as well as
  information transmission and forgetting, local information generation,
  information generation over i links and vaccination.
  
  */
  template <class Graph>
  class VaccinationSIRS :
    public EpiModel<State, Graph>
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

    virtual Model<Graph>* clone() const { return new VaccinationSIRS<Graph>(*this); }

    void Init(const po::variables_map& vm,
              std::vector<StatRecorder<Graph>*>& rec);
  
    unsigned int getNodeEvents(eventList& events, State* currentState,
                               unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State* currentState,
                               unsigned int edge, State* currentNbState,
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

    bool isForgetting(State before_state, State after_state) const
    { return ((getInfo(before_state) == Informed) &&
              (getInfo(after_state) == Uninformed)); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State state) const
    { return (state.getState() % 3); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State state) const
    { return (state.getState() / 3); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(State dState, State iState) const
    { return dState.getState()+iState.getState()*3; }
    
    std::vector<StatRecorder<Graph>*> 
    getStatRecorders(const po::variables_map& vm) const
    {;}
  
  private:

    unsigned int gamma[2]; //!< Recovery rates.
    unsigned int delta[2]; //!< Loss of immunity rates.
    unsigned int beta[2][2]; //!< Infection rates.
    unsigned int alpha; //!< Information transmission rate.
    unsigned int lambda; //!< Rate of forgetting.
    unsigned int omega; //!< Local infromation generation rate.
    unsigned int nu; //!< Information generation rate over i-edges.
    unsigned int theta; //!< Vaccination rate of informed susceptibles.
    double sigma; //!< Ratio between informed and uninformed susceptibility.
  
  };
  
}

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
template <class Graph>
Models::VaccinationSIRS<Graph>::VaccinationSIRS(unsigned int v)
  : EpiModel<State, Graph>("VaccinationSIRS", v)
{
  /*************************************/
  // define vertex classes
  /************************************/
//   // susceptible uninformed
//   vertexStates.push_back(Label("S-","00;32", 0, "fillcolor=\"royalblue4\""));
//   // infected uninformed
//   vertexStates.push_back(Label("I-","00;31", 1, "fillcolor=\"red4\""));
//   // recovered uninformed
//   vertexStates.push_back(Label("R-","00;34", 2, "fillcolor=\"green4\""));
//   // susceptible informed
//   vertexStates.push_back(Label("S+","01;32", 3, "fillcolor=\"royalblue\""));
//   // infected informed
//   vertexStates.push_back(Label("I+","01;31", 4, "fillcolor=\"red\""));
//   // recovered informed
//   vertexStates.push_back(Label("R+","01;34", 5, "fillcolor=\"green\""));

  /*************************************/
  // define edge types
  /************************************/
  this->edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  this->edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

  /*************************************/
  // define model parameters
  /************************************/
  this->model_options.add_options()
    ("beta--", po::value<double>(),
     "disease transmission rate uninformed->uninformed")
    ("beta+-", po::value<double>(),
     "disease transmission rate informed->uninformed")
    ("beta-+", po::value<double>(),
     "disease transmission rate uninformed->informed")
    ("beta++", po::value<double>(),
     "disease transmission rate informed->informed")
    ("gamma-", po::value<double>(),
     "recovery rate of uninformed")
    ("gamma+", po::value<double>(),
     "recovery rate of informed")
    ("delta-", po::value<double>(),
     "loss of immunity rate of uninformed")
    ("delta+", po::value<double>(),
     "loss of immunity rate of informed")
    ("alpha", po::value<double>(),
     "information transmission rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("omega", po::value<double>(),
     "local information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ("theta", po::value<double>(),
     "vaccination rate")
    ("sigma", po::value<double>(),
     "ratio between uninformed/uninformed susceptibilities")
    ;

  this->rates.insert(std::make_pair("beta--", &beta[0][0]));
  this->rates.insert(std::make_pair("beta+-", &beta[1][0]));
  this->rates.insert(std::make_pair("beta-+", &beta[0][1]));
  this->rates.insert(std::make_pair("beta++", &beta[1][1]));
  this->rates.insert(std::make_pair("gamma-", &gamma[0]));
  this->rates.insert(std::make_pair("gamma+", &gamma[1]));
  this->rates.insert(std::make_pair("delta-", &delta[0]));
  this->rates.insert(std::make_pair("delta+", &delta[1]));
  this->rates.insert(std::make_pair("alpha", &alpha));
  this->rates.insert(std::make_pair("nu", &nu));
  this->rates.insert(std::make_pair("omega", &omega));
  this->rates.insert(std::make_pair("lambda", &lambda));
  this->rates.insert(std::make_pair("theta", &theta));
  this->params.insert(std::make_pair("sigma", &sigma));

}

//----------------------------------------------------------
template <class Graph>
void Models::VaccinationSIRS<Graph>::Init(const po::variables_map& vm,
                                          std::vector<StatRecorder<Graph>*>& rec)
{
  Model<Graph>::Init(vm, rec);
  if (vm.count("sigma")) {
    // if sigma is defined, beta+- and beta++ are overwritten
    beta[1][0]=static_cast<unsigned int>(sigma * beta[0][0] + .5);
    beta[1][1]=static_cast<unsigned int>(sigma * beta[0][1] + .5);
    if (this->getVerbose() >= 1) {
      std::cout << "sigma given, setting beta+-=" << beta[1][0]
                << " and beta++=" << beta[1][1] << std::endl;
    }
  }
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::VaccinationSIRS<Graph>::getNodeEvents(eventList& events,
                                                    State* currentState,
                                                    unsigned int nb) const
{
   State* state = currentState;

   unsigned int rateSum(0);

   // loss of immunity
   if (getDisease(state->getState()) == Recovered) {
      Event immunityLoss;
      immunityLoss.rate = delta[getInfo(state->getState())];
      immunityLoss.newState = 
	this->newState(State(getState(Susceptible, getInfo(state->getState()))));
      immunityLoss.nb = nb;
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding loss of immunity event with rate " 
	            << immunityLoss.rate/1e+4 << std::endl;
        }
      }
   } else
   if (getDisease(state->getState()) == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma[getInfo(state->getState())];
      recovery.newState = 
	this->newState(State(getState(Recovered, getInfo(state->getState()))));
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate/1e+4
                    << std::endl;
        }
      }
      if (getInfo(state->getState()) == Uninformed) {
        // local information generation
        Event localInfo;
        localInfo.rate = omega;
        localInfo.newState = 
	  this->newState(State(getState(getDisease(state->getState()), Informed)));
        localInfo.nb = nb;
        if (localInfo.rate > 0) {
          events.push_back(localInfo);
          rateSum += localInfo.rate;
          if (this->getVerbose() >= 2) {
            std::cout << "Adding local information event with rate " 
	              << localInfo.rate/1e+4 << std::endl;
          }
        }
      }
   }
   // information loss
   if (getInfo(state->getState()) == Informed) {
      Event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState = 
	this->newState(State(getState(getDisease(state->getState()), Uninformed)));
      infoLoss.nb = nb;
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding information loss event with rate " 
	            << infoLoss.rate/1e+4 << std::endl;
        }
      }
      // vaccination
      if (getDisease(state->getState()) == Susceptible) {
        Event vaccination;
        vaccination.rate = theta;
        vaccination.newState = 
	  this->newState(State(getState(Recovered, Informed)));
        vaccination.nb = nb;
        if (vaccination.rate > 0) {
          events.push_back(vaccination);
          rateSum += vaccination.rate;
          if (this->getVerbose() >= 2) {
            std::cout << "Adding vaccination event with rate "
                      << vaccination.rate/1e+4 << std::endl;
          }
        }
      }
   }

   return rateSum;
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::VaccinationSIRS<Graph>::getEdgeEvents(eventList& events,
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
      if (getDisease(state->getState()) == Susceptible &&
          getInfo(nbState->getState()) == Infected) {
         Event infection;
         infection.rate = beta[getInfo(state->getState())][getInfo(nbState->getState())];
         infection.newState = 
	   this->newState(State(getState(Infected, getInfo(state->getState()))));
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
   } else if (edge == Information) {
      // information transmission
      if (getInfo(state->getState()) == Uninformed && 
	  getInfo(nbState->getState()) == Informed) {
         Event infoTransmission;
         infoTransmission.rate = alpha;
         infoTransmission.newState =
           this->newState(State(getState(getDisease(state->getState()), Informed)));
         infoTransmission.nb = nb;
         infoTransmission.et = edge;
         if (infoTransmission.rate > 0) {
           events.push_back(infoTransmission);
           rateSum += infoTransmission.rate;
           if (this->getVerbose() >= 2) {
             std::cout << "Adding information transmission event with rate " 
                       << infoTransmission.rate/1e+4 << std::endl;
           }
         }
      }
      // information generation
      if (getInfo(state->getState()) == Uninformed && 
	  getInfo(nbState->getState()) == Infected) {
         Event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState = 
	   this->newState(State(getState(getDisease(state->getState()), Informed)));
         infoGeneration.nb = nb;
         infoGeneration.et = edge;
         if (infoGeneration.rate > 0) {
           events.push_back(infoGeneration);
           rateSum += infoGeneration.rate;
           if (this->getVerbose() >= 2) {
             std::cout << "Adding information generation event with rate " 
                       << infoGeneration.rate/1e+4 << std::endl;
           }
         }
      }
   }

   return rateSum;
}

#endif
