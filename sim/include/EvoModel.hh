/*! \file EvoModel.hh
  \brief The Models::EvoModel class.
*/
#ifndef EVOMODEL_HH
#define EVOMODEL_HH

#include "Model.hh"

class EvoModelState :
  virtual public State
{
public:
  
  EvoModelState(double trait = .0)
    : State(), trait(trait), pendingConnection(false), newNb(0) {;}
  ~EvoModelState() {;}
  
  virtual State* clone() const { return new EvoModelState(*this); }
  
  double getTrait() const { return trait; }
  void setTrait(double t) { trait = t; }
  void setPending() { pendingConnection = true; }
  void clearPending() { pendingConnection = false; }
  bool getPending() const { return pendingConnection; }
  unsigned int getNewNb() const { return newNb; }
  void setNewNb(unsigned int n) { newNb = n; }
  
private:
  double trait;
  bool pendingConnection;
  unsigned int newNb;
};

namespace Models {
  
  /*! \brief Class implementing the network evolution model.

  */
  template <class Graph>
  class EvoModel :
    public Model<Graph>
  {

  public:

    //! Possible edge types (base layer / interaction layer)
    enum edgeTypesEnum {Base, Interaction};
    
    EvoModel(unsigned int v = 0);
    ~EvoModel() {;}
    
    virtual Model<Graph>* clone() const { return new EvoModel<Graph>(*this); }

    void Init(const po::variables_map& vm,
              std::vector<StatRecorder<Graph>*>& rec);

    virtual EvoModelState* newState() const
    { EvoModelState* d = new EvoModelState(); return d; }

    unsigned int getEdgeEvents(eventList& events, State* state,
                               unsigned int edge, State* nbState,
                               unsigned int nb) const;
    
    virtual std::string printState(State* s) const;
    std::vector<unsigned int> getRGB(State* s) const;
    
    
  private:
    double assortativity; //!< strength of assortativity
    bool cluster;
  };

  template <class Graph>
  EvoModel<Graph>::EvoModel(unsigned int)
  {
    /*************************************/
    // define edge type
    /************************************/
    // just one edgetype
    this->edgeTypes.push_back(Label("b", "", 0, "style=\"invis\""));
    this->edgeTypes.push_back(Label("i", "", 1, "style=\"solid\""));

    /*************************************/
    // define model parameters
    /************************************/
    this->model_options.add_options()
      ("assortativity,a",po::value<double>()->default_value(0.),
       "assortativity in rewiring")
      ("clustering,c",
       "base assortativity on clustering")
      ;
  }

  template <class Graph>
  void EvoModel<Graph>::Init
  (const po::variables_map& vm, std::vector<StatRecorder<Graph>*>& rec)
  {
    Model<Graph>::Init(vm, rec);

    assortativity = vm["assortativity"].as<double>();
    cluster = vm.count("clustering");

    // make vertices red
    this->vertexStates.push_back
      (Label("0", "01;31", 0, "", Label::rgbColour(255,0,0)));
  }

  //----------------------------------------------------------
  template <class Graph>
  unsigned int Models::EvoModel<Graph>::getEdgeEvents(eventList& events,
                                                      State* currentState,
                                                      unsigned int edge,
                                                      State* currentNbState,
                                                      unsigned int nb) const
  {
    EvoModelState* state = dynamic_cast<EvoModelState*>(currentState);   
    EvoModelState* nbState = dynamic_cast<EvoModelState*>(currentNbState);
    
    unsigned int rateSum(0);
    if (edge == Base) { // nothing happens on the interaction layer
      Event interaction;
      interaction.rate =
        static_cast<unsigned int>
        (1e+4*(pow(1-fabs(state->getTrait() - nbState->getTrait()),
                     assortativity)));
      EvoModelState* targetState = newState();
      targetState->setTrait(state->getTrait());
      targetState->setPending();
      targetState->setNewNb(nb);
      interaction.newState = targetState;
      interaction.nb = nb;
      interaction.et = edge;
      if (interaction.rate > 0) {
        events.push_back(interaction);
        rateSum += interaction.rate;
        if (this->getVerbose() >= 2) {
          std::cout << "Adding interaction event with rate "
                    << interaction.rate/1e+4 << std::endl;
        }
      }
    }

    return rateSum;
  }
  

  template <class Graph>
  std::vector<unsigned int> EvoModel<Graph>::
  getRGB(State* s) const
  {
    EvoModelState* myState = dynamic_cast<EvoModelState*>(s);
    std::vector<unsigned int> rgb;
    double darkening = 1 - (1 - myState->getTrait())*4./5.;
    
    rgb.push_back
      (static_cast<unsigned int>
       (this->getVertexState(myState->getState()).getRGB(0) * darkening));
    rgb.push_back
      (static_cast<unsigned int>
       (this->getVertexState(myState->getState()).getRGB(1) * darkening));
    rgb.push_back
      (static_cast<unsigned int>
       (this->getVertexState(myState->getState()).getRGB(2) * darkening));

    return rgb;
  }


  template <class Graph>
  std::string EvoModel<Graph>::
  printState(State* s) const
  {
    EvoModelState* myState = dynamic_cast<EvoModelState*>(s);
    std::stringstream ss;
    std::streamsize prec = ss.precision();
    ss << std::setprecision(2) << myState->getTrait();
    ss << std::setprecision(prec);
    ss << "<" << myState->getNewNb() << ">";
    ss.unsetf(std::ios::fixed);
    return ss.str();
  }
  
}

#endif
