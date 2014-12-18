/*! \file SIRInfo.hh
  \brief The Models::SIRInfo class.
*/
#ifndef SIRINFO_HH
#define SIRINFO_HH

#include <algorithm>

#include "EpiModel.hh"

namespace Models {

  /*! \brief Class implementing an SIRS model with information.

  This defines 6 vertex states (S-,I-,R-,S+,I+,R+) and 2 edge types (d and i),
  with information-dependent infection, recovery and loss of immunity, as well as
  information transmission and forgetting, local information generation and
  information generation over i links.

  */
  template <class Graph>
  class SIRInfo :
    public EpiModel<State, Graph>
  {

    using EpiModel<State, Graph>::newState;

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible information states.
    enum infoStatesEnum {Uninformed, Informed, Stifling};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};

  public:

    SIRInfo(unsigned int v = 0);
    virtual ~SIRInfo() {;}

    virtual Model<Graph>* clone() const { return new SIRInfo<Graph>(*this); }

    virtual State* newState(unsigned int disease, unsigned int info) const
    { State* s = new State(info*3+disease); return s; }

    void Init(const po::variables_map& vm,
              std::vector<StatRecorder<Graph>*>& rec);

    unsigned int getNodeEvents(eventList& events, State* currentState,
                         unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State* currentState,
                         unsigned int edge, State* currentNbState,
                         unsigned int nb) const;

    bool isInfection(State* before_state, State* after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }
    bool isRecovery(State* before_state, State* after_state) const
    { return ((getDisease(before_state) == Infected) &&
              (getDisease(after_state) != Infected)); }
    bool isInformation(State* before_state, State* after_state) const
    { return ((getInfo(before_state) == Uninformed) &&
              (getInfo(after_state) == Informed)); }
    bool isForgetting(State*  before_state, State* after_state) const
    { return ((getInfo(before_state) == Informed) &&
              (getInfo(after_state) == Uninformed)); }
    bool isInfected(State* state) const
    { return (getDisease(state) == Infected); }

    bool isInformed(State* state) const
    { return (getInfo(state) == Informed); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State* state) const
    { return (state->getState() % 3); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State* state) const
    { return (state->getState() / 3); }
    //! Get the full vertex state from disease and information parts.
    unsigned int getState(State* dState, State* iState) const
    { return dState->getState()+iState->getState()*3; }

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
    double sigma; //!< Ratio between informed and uninformed susceptibility.

  };

}

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
template <class Graph>
Models::SIRInfo<Graph>::SIRInfo(unsigned int v)
  : EpiModel<State, Graph>("SIRInfo", v)
{
  /*************************************/
  // define vertex classes
  /************************************/
  // susceptible uninformed
  this->vertexStates.push_back(Label("S-","00;32", 0, "fillcolor=\"royalblue4\"", Label::rgbColour(0,0,51)));
  // infected uninformed
  this->vertexStates.push_back(Label("I-","00;31", 1, "fillcolor=\"red4\"", Label::rgbColour(51,0,0)));
  // recovered uninformed
  this->vertexStates.push_back(Label("R-","00;34", 2, "fillcolor=\"green4\"", Label::rgbColour(0,51,0)));
  // susceptible informed
  this->vertexStates.push_back(Label("S+","01;32", 3, "fillcolor=\"royalblue\"", Label::rgbColour(0,0,255)));
  // infected informed
  this->vertexStates.push_back(Label("I+","01;31", 4, "fillcolor=\"red\"", Label::rgbColour(255,0,0)));
  // recovered informed
  this->vertexStates.push_back(Label("R+","01;34", 5, "fillcolor=\"green\"", Label::rgbColour(0,255,0)));
  // susceptible stifling
  this->vertexStates.push_back(Label("S%","04;32", 6, "fillcolor=\"royalblue2\"", Label::rgbColour(0,0,153)));
  // infected stifling
  this->vertexStates.push_back(Label("I%","04;31", 7, "fillcolor=\"red2\"", Label::rgbColour(153,0,0)));
  // recovered stifling
  this->vertexStates.push_back(Label("R%","04;34", 8, "fillcolor=\"green2\"", Label::rgbColour(0,153,0)));

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
    ("sigma", po::value<double>(),
     "ratio between uninformed/uninformed susceptibilities")
    ;

  /*************************************/
  // assign model parameters and rates to variables
  /************************************/
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
  this->params.insert(std::make_pair("sigma", &sigma));
}

//----------------------------------------------------------
template <class Graph>
void Models::SIRInfo<Graph>::Init(const po::variables_map& vm,
                                   std::vector<StatRecorder<Graph>*>& rec)
{
  Model<Graph>::Init(vm, rec);
  if (vm.count("sigma")) {
    // if sigma is defined, beta+- and beta++ and beta-+ are overwritten
    beta[1][0]=static_cast<unsigned int>(sigma * beta[0][0] + .5);
    beta[1][1]=static_cast<unsigned int>(sigma * beta[0][0] + .5);
    beta[0][1]=static_cast<unsigned int>(sigma * beta[0][0] + .5);
    if (this->getVerbose() >= 1) {
      std::cout << "sigma given, setting beta+-=" << beta[1][0]/1e+4
                << " and beta++=" << beta[1][1]/1e+4 << " and beta-+="
                << beta[0][1]/1e+4 << std::endl;
    }
  }
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::SIRInfo<Graph>::getNodeEvents(eventList& events,
                                             State* currentState,
                                             unsigned int nb) const
{
   State* state = currentState;
   unsigned int rateSum(0);

   if (getDisease(state) == Recovered) {
     // loss of immunity
      Event immunityLoss;
      unsigned int infoState = std::min(getInfo(state), static_cast<unsigned int>(1));
      immunityLoss.rate = delta[infoState];
      immunityLoss.newState = this->newState(Susceptible, getInfo(state));
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
   if (getDisease(state) == Infected) {
      // recovery
      Event recovery;
      unsigned int infoState = std::min(getInfo(state), static_cast<unsigned int>(1));
      recovery.rate = gamma[infoState];
      recovery.newState = this->newState(Recovered, getInfo(state));
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate/1e+4
                    << std::endl;
        }
      }
      if (getInfo(state) == Uninformed) {
        // local information generation
        Event localInfo;
        localInfo.rate = omega;
        localInfo.newState = this->newState(getDisease(state), Informed);
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
   if (getInfo(state) == Informed) {
      Event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState =
        this->newState(getDisease(state), Stifling);
      infoLoss.nb = nb;
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding information loss event with rate "
	            << infoLoss.rate/1e+4 << std::endl;
        }
      }
   }

   return rateSum;
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::SIRInfo<Graph>::getEdgeEvents(eventList& events,
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
         unsigned int infoState = std::min(getInfo(state), static_cast<unsigned int>(1));
         unsigned int nbInfoState = std::min(getInfo(nbState), static_cast<unsigned int>(1));
         infection.rate = beta[infoState][nbInfoState];
         infection.newState = this->newState(Infected, getInfo(state));
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
      if (getInfo(state) == Uninformed &&
	  getInfo(nbState) == Informed) {
         Event infoTransmission;
         infoTransmission.rate = alpha;
         infoTransmission.newState =
           this->newState(getDisease(state), Informed);
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
      if (getInfo(state) == Uninformed &&
	  getDisease(nbState) == Infected) {
         Event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState =
	   this->newState(getDisease(state), Informed);
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
