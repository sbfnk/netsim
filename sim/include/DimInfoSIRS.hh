/*! \file DimInfoSIRS.hh
  \brief The Models::DimInfoSIRS class.
*/
#ifndef DIMINFOSIRS_HH
#define DIMINFOSIRS_HH

#include "EpiModel.hh"

class DimInfoState :
  virtual public State
{
public:
  
  DimInfoState(unsigned int b = 0, double i = 0.)
    : State(b), info_quality(i) {;}
  ~DimInfoState() {;}

  virtual State* clone() { return new DimInfoState(*this);}

  double getInfo() const { return info_quality; }
  void setInfo(double i) { info_quality = i; }

private:
  double info_quality;
};
  
namespace Models {

  /*! \brief Class implementing an SIRS model with information.

  This defines 3 vertex states (S,I,R), a real value associated with quality
  of information, and 2 edge types (d and i), with information-dependent
  infection, as well as  information transmission and forgetting and local
  information generation. Information diminishes at each transmission by a factor rho.
  
  */
  template <class Graph>
  class DimInfoSIRS :
    public EpiModel<DimInfoState, Graph>
  {

    //! Possible disease states.
    enum diseaseStatesEnum {Susceptible,Infected,Recovered};
    //! Possible edge types.
    enum edgeTypesEnum {Disease, Information};

  public:

    DimInfoSIRS(unsigned int v = 0);
    ~DimInfoSIRS() {;}

    virtual Model<Graph>* clone() { return new DimInfoSIRS<Graph>(*this); }

    void Init(const po::variables_map& vm, std::vector<StatRecorder<Graph>*>& rec);

    virtual DimInfoState* newState() const
    { DimInfoState* d = new DimInfoState(); d->setInfo(defaultInfo); return d; }
  
    virtual DimInfoState* newState(unsigned int red, unsigned int green,
                                   unsigned int blue) const
    {
      DimInfoState* state =
        new DimInfoState
        (this->getStateFromColour((red > 0)*1 + (green > 0)*2 + (blue > 0) * 4));
      double detail =
        1 - 5/4. * (1 - (red+green+blue)/ static_cast<double>
                    (this->getVertexStates()[state->getState()].getRGB(0)*256 +
                     this->getVertexStates()[state->getState()].getRGB(1)*256 +
                     this->getVertexStates()[state->getState()].getRGB(2)*256 +
                     255));
      // brush up rounding errors
      detail = detail > 0 ? detail : 0;
      state->setInfo(detail);
      return state;
    }
    
    unsigned int getNodeEvents(eventList& events, State* state,
                               unsigned int nb) const;
    unsigned int getEdgeEvents(eventList& events, State* state,
                               unsigned int edge, State* nbState,
                               unsigned int nb) const;

    bool isInfection(State* before_state, State* after_state) const
    { return ((getDisease(before_state) == Susceptible) &&
              (getDisease(after_state) == Infected)); }
    bool isRecovery(State* before_state, State* after_state) const
    { return ((getDisease(before_state) == Infected) &&
              (getDisease(after_state) != Infected)); }
    bool isInformation(State* before_state, State* after_state) const
    { return (dynamic_cast<DimInfoState*>(before_state)->getInfo() == 0 &&
              dynamic_cast<DimInfoState*>(after_state)->getInfo() > 0); }
    bool isForgetting(State* before_state, State* after_state) const
    { return (dynamic_cast<DimInfoState*>(before_state)->getInfo() > 0 &&
              dynamic_cast<DimInfoState*>(after_state)->getInfo() == 0); }

    bool isInfected(State* state) const
    { return (getDisease(state) == Infected); }

    bool isInformed(State* state) const
    { return (dynamic_cast<DimInfoState*>(state)->getInfo() > 0); }

    //! Get the disease part of a full vertex state.
    unsigned int getDisease(State* state) const
    { return state->getState(); }
    //! Get the information part of a full vertex state.
    unsigned int getInfo(State* state) const
    { return dynamic_cast<DimInfoState*>(state)->getInfo(); }

    virtual std::string printState(State* s) const
    {
      DimInfoState* myState = dynamic_cast<DimInfoState*>(s);
      std::stringstream ss;
      std::streamsize prec = ss.precision();
      ss << this->getVertexState(myState->getState()) << ","
         << std::setprecision(2) << myState->getInfo();
      ss << std::setprecision(prec);
      ss.unsetf(std::ios::fixed);
      return ss.str();
    }

    std::vector<unsigned int> getRGB(State* s) const
    {
      DimInfoState* myState = dynamic_cast<DimInfoState*>(s);
      std::vector<unsigned int> rgb;
      double darkening = 1 - (1 - myState->getInfo())*4./5.;
      
      rgb.push_back
        (static_cast<unsigned int>
         (this->getVertexStates()[myState->getState()].getRGB(0) * darkening));
      rgb.push_back
        (static_cast<unsigned int>
         (this->getVertexStates()[myState->getState()].getRGB(1) * darkening));
      rgb.push_back
        (static_cast<unsigned int>
         (this->getVertexStates()[myState->getState()].getRGB(2) * darkening));
      
      return rgb;
    }

  private:

    unsigned int gamma; //!< Recovery rates.
    unsigned int delta; //!< Loss of immunity rates.
    unsigned int beta; //!< Infection rates.
    unsigned int alpha; //!< Information transmission rate.
    unsigned int lambda; //!< Rate of forgetting.
    unsigned int omega; //!< Local infromation generation rate.
    unsigned int nu; //!< Information generation rate over i-edges.

    double rho; //!< Ratio between informed and uninformed susceptibility.
    double threshold; //!< Threshold below which information will not be passed
    double dim_alpha; //!< Whether to diminish alpha or not

    double defaultInfo;

  };

}

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
template <class Graph>
Models::DimInfoSIRS<Graph>::DimInfoSIRS(unsigned int v)
  : EpiModel<DimInfoState, Graph>("DimInfoSIRS", v)
{
  
  /*************************************/
  // define vertex classes
  /************************************/
  // susceptible
  this->vertexStates.push_back(Label("S","01;32", 0, "", Label::rgbColour(0,0,255)));
  // infected 
  this->vertexStates.push_back(Label("I","01;31", 1, "", Label::rgbColour(255,0,0)));
  // recovered
  this->vertexStates.push_back(Label("R","01;34", 2, "", Label::rgbColour(0,255,0)));

  /*************************************/
  // define edge types
  /************************************/
  this->edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  this->edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

  /*************************************/
  // define model parameters
  /************************************/
  this->model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate uninformed->uninformed")
    ("rho", po::value<double>(),
     "loss of quality in information")
    ("gamma", po::value<double>(),
     "recovery rate")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ("alpha", po::value<double>(),
     "information transmission rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("omega", po::value<double>(),
     "local information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ("threshold", po::value<double>(),
     "threshold below which not to consider information transmission")
    ("dim_alpha", po::value<double>()->default_value(0.),
     "whether alpha diminishes or not")
    ("effective,e", po::value<double>(), 
     "write effective data at arg timesteps")
    ("info-dist,i", po::value<double>(),
     "create information distribution in the dist directory at arg timesteps")
    ("info-dis-corr,c", po::value<double>(),
     "create correlation between information and infection in the corr directory at arg timesteps")
    ("risk-info,r",po::value<double>(),
     "create informedness of population at risk in risk directory at arg timesteps")
    ("default-info", po::value<double>()->default_value(0.),
     "default information quality in the population")
    ;

  /*************************************/
  // assign model parameters to variables
  /************************************/
  this->rates.insert(std::make_pair("beta", &beta));
  this->rates.insert(std::make_pair("gamma", &gamma));
  this->rates.insert(std::make_pair("delta", &delta));
  this->rates.insert(std::make_pair("alpha", &alpha));
  this->rates.insert(std::make_pair("nu", &nu));
  this->rates.insert(std::make_pair("omega", &omega));
  this->rates.insert(std::make_pair("lambda", &lambda));

  this->params.insert(std::make_pair("rho", &rho));
  this->params.insert(std::make_pair("threshold", &threshold));
  this->params.insert(std::make_pair("dim_alpha", &dim_alpha));

}

//----------------------------------------------------------
template <class Graph>
void Models::DimInfoSIRS<Graph>::Init(const po::variables_map& vm,
                                      std::vector<StatRecorder<Graph>*>& rec)
{
  EpiModel<DimInfoState, Graph>::Init(vm, rec);
  if (vm.count("effective")) {
    rec.push_back
      (new StatRecorder<Graph>
       (new write_effective_data<Graph, DimInfoSIRS>(*this), 
        vm["effective"].as<double>()));
  }
  if (vm.count("info-dist")) {
    rec.push_back(new StatRecorder<Graph>
                  (new write_detail_dist<Graph, DimInfoSIRS>,
                   vm["info-dist"].as<double>()));
  }
  if (vm.count("info-dis-corr")) {
    rec.push_back(new StatRecorder<Graph>
                  (new write_info_dis_corr<Graph, DimInfoSIRS>(*this, 0),
                   vm["info-dis-corr"].as<double>()));
  }
  if (vm.count("risk-info")) {
    rec.push_back(new StatRecorder<Graph>
                  (new write_risk_info<Graph, DimInfoSIRS>,
                   vm["risk-info"].as<double>()));
  }
  defaultInfo = vm["default-info"].as<double>();
}

//----------------------------------------------------------
template <class Graph>
unsigned int Models::DimInfoSIRS<Graph>::getNodeEvents(eventList& events,
                                                State* currentState,
                                                unsigned int nb) const
{

   DimInfoState* state = dynamic_cast<DimInfoState*>(currentState);

   unsigned int rateSum(0);

   if (getDisease(state) == Recovered) {
     // loss of immunity
      Event immunityLoss;
      immunityLoss.rate = delta;
      DimInfoState* targetState = newState();
      targetState->setState(Susceptible);
      targetState->setInfo(state->getInfo());
      immunityLoss.newState = targetState;
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
      recovery.rate = gamma;
      DimInfoState* targetState = newState();
      targetState->setState(Recovered);
      targetState->setInfo(state->getInfo());
      recovery.newState = targetState;
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate/1e+4
                    << std::endl;
        }
      }
      if (state->getInfo() < 1) {
        // local information generation
        Event localInfo;
        localInfo.rate = omega;
        DimInfoState* targetState = newState();
        targetState->setState(state->getState());
        targetState->setInfo(1.);
        localInfo.newState = targetState;
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
   if (state->getInfo() > 0) {
      Event infoLoss;
      infoLoss.rate = lambda;
      double detailUpdate = state->getInfo() * rho;
      if (detailUpdate <= threshold) detailUpdate = 0;
      DimInfoState* targetState = newState();
      targetState->setState(state->getState());
      targetState->setInfo(detailUpdate);
      infoLoss.newState = targetState;
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
unsigned int Models::DimInfoSIRS<Graph>::getEdgeEvents(eventList& events,
                                                State* currentState,
                                                unsigned int edge,
                                                State* currentNbState,
                                                unsigned int nb) const
{
   DimInfoState* state = dynamic_cast<DimInfoState*>(currentState);   
   DimInfoState* nbState = dynamic_cast<DimInfoState*>(currentNbState);
   
   unsigned int rateSum(0);
   if (edge == Disease) {
     // infection
     if (getDisease(state) == Susceptible &&
         getDisease(nbState) == Infected) {
         Event infection;
         infection.rate =
           static_cast<unsigned int>((1 - state->getInfo()) * beta + .5);
         DimInfoState* targetState = newState();
         targetState->setState(Infected);
         targetState->setInfo(state->getInfo());
         infection.newState = targetState;
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
      if (state->getInfo() < nbState->getInfo() * rho && 
	  nbState->getInfo() > threshold) {
         Event infoTransmission;
         if (dim_alpha) {
           infoTransmission.rate =
             static_cast<unsigned int>(nbState->getInfo() * alpha + .5);
         } else {
           infoTransmission.rate = alpha;
         }
         DimInfoState* targetState = newState();
         targetState->setState(state->getState());
         targetState->setInfo(nbState->getInfo() * rho);
         infoTransmission.newState = targetState;
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
      if (state->getInfo() < 1 && 
	  getDisease(state) == Susceptible && 
          getDisease(nbState) == Infected) {
         Event infoGeneration;
         infoGeneration.rate = nu;
         DimInfoState* targetState = newState();
         targetState->setState(state->getState());
         targetState->setInfo(1);
         infoGeneration.newState = targetState;
         infoGeneration.et = edge;
         infoGeneration.nb = nb;
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
