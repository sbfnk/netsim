/*! \file DimInfoSIRS.hh
  \brief The Models::GroupFormModel class.
*/
#ifndef GROUPFORMMODEL_HH
#define GROUPFORMMODEL_HH

#include "Model.hh"

#include <math.h>
#include <algorithm>

class GroupFormState :
  public State
{
public:
  
  GroupFormState(unsigned int b = 0, std::vector<double> t = std::vector<double>(0),
                 double v = 1., double a = 1.)
    : State(b), trait_vector(t), volatility(v), acceptance(a) {;}
  ~GroupFormState() {;}
  
  virtual State* clone() const { return new GroupFormState(*this); }
  
  const std::vector<double>& getTraits() const { return trait_vector; }
  const double getTrait(unsigned int i) const
  {
    if (i+1 > trait_vector.size()) {
      return 0.;
    } else {
      return trait_vector[i];
    }
  }

  const double getVolatility() const
  {
    return volatility;
  }
  
  const double getAcceptance() const
  {
    return acceptance;
  }
  
  void setTrait(std::vector<double> t) { trait_vector = t; }
  void setVolatility(double v) { volatility = v; }
  void setAcceptance(double a) { acceptance = a; }

private:
  std::vector<double> trait_vector;
  double volatility;
  double acceptance;
};

namespace Models {
  
  /*! \brief Class implementing the group formation model.

  */
  template <class Graph>
  class GroupFormModel :
    public Model<Graph>
  {

  public:

    typedef GroupFormState StateType;

    GroupFormModel(unsigned int v = 0);
    ~GroupFormModel() {;}
    
    virtual Model<Graph>* clone() const { return new GroupFormModel<Graph>(*this); }

    virtual StateType* newState() const
    { return new StateType(); }
    
    virtual void Init(const po::variables_map& vm,
                      std::vector<StatRecorder<Graph>*>& rec);
  
    virtual std::string printState(State* s) const;
    std::vector<unsigned int> getRGB(State* s) const;

    void addState()
    {
      unsigned i = this->vertexStates.size();
      std::stringstream s;
      s << i;
      switch (i % 3) {
       case 0: 
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;32", i, "", Label::rgbColour(0,0,255)));
        break;
       case 1: 
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;31", i, "", Label::rgbColour(255,0,0)));
        break;
       case 2:
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;34", i, "", Label::rgbColour(0,255,0)));
        break;
      }
    }
  
    unsigned int getStates() const
    { return states; }

    unsigned int getTraitDim() const
    { return trait_dim; }

    double getAcceptance() const
    { return acceptance; }

    double accept(State* s1, State* s2) const
    {
      StateType* state1 = dynamic_cast<StateType*>(s1);
      StateType* state2 = dynamic_cast<StateType*>(s2);
      return acceptance*(1-distance(state1, state2));
    }

    double distance(State* s1, State* s2) const
    {
      StateType* state1 = dynamic_cast<StateType*>(s1);
      StateType* state2 = dynamic_cast<StateType*>(s2);
      return distance(state1, state2);
    }

    double getVolatility(State* s) const
    {
      StateType* state = dynamic_cast<StateType*>(s);
      return getVolatility(state);
    }

    double getAcceptance(State* s) const
    {
      StateType* state = dynamic_cast<StateType*>(s);
      return getAcceptance(state);
    }

  private:

    double distance(StateType* s1, StateType* s2) const
    {
      double dist = 0.;
      for (unsigned int i = 0; i < s1->getTraits().size(); ++i) {
        double diff = 0.;
        if (close) {
          diff = std::min(fabs(s1->getTrait(i)-s2->getTrait(i)),
                          1 - fabs(s1->getTrait(i)-s2->getTrait(i)));
        } else {
          diff = s1->getTrait(i) - s2->getTrait(i);
        }                              
        dist += pow(diff, 2);
      }
      return sqrt(dist);
    }

    double getVolatility(StateType* s) const
    {
      return s->getVolatility();
    }

    double getAcceptance(StateType* s) const
    {
      return s->getAcceptance();
    }

    unsigned int states; //!< number of states (groups).
    unsigned int trait_dim; //!< Dimension of the trait vector
    double acceptance; //!< Base property of acceptance.
    bool close;
  };

  template <class Graph>
  GroupFormModel<Graph>::GroupFormModel(unsigned int)
  {
    this->model_options.add_options()
      ("states", po::value<unsigned int>()->default_value(0),
       "number of possible states (groups?)")
      ("traitdim", po::value<unsigned int>()->default_value(1),
       "dimensionality of the trait vector")
      ("alpha", po::value<double>()->default_value(1.),
       "alpha")
      ("avg-trait", po::value<double>(),
       "write average trait for each state at arg timesteps")
      ("close",
       "close ring at [0,1] when calculating distances")
      ;

    this->edgeTypes.push_back(Label("x", "", 0, "style=\"solid\""));

    this->params.insert(std::make_pair("alpha", &acceptance));
    this->intParams.insert(std::make_pair("states", &states));
    this->intParams.insert(std::make_pair("traitdim", &trait_dim));
  }

  template <class Graph>
  void GroupFormModel<Graph>::Init
  (const po::variables_map& vm, std::vector<StatRecorder<Graph>*>& rec)
  {
    Model<Graph>::Init(vm, rec);

    close = vm.count("close");

    if (vm.count("avg-trait")) {
      rec.push_back
        (new StatRecorder<Graph>
         (new write_avg_trait<Graph, GroupFormModel>(*this), 
          vm["avg-trait"].as<double>()));
    }
    
    /*************************************/
    // define vertex classes
    /************************************/
    // ungrouped state
    this->vertexStates.push_back
      (Label("0", "01;37", 0, "", Label::rgbColour(0,0,0)));
    for (unsigned int i = 1; i < states+1; ++i) {
      addState();
    }
  }

  template <class Graph>
  std::vector<unsigned int> GroupFormModel<Graph>::
  getRGB(State* s) const
  {
    StateType* myState = dynamic_cast<StateType*>(s);
    std::vector<unsigned int> rgb;
    
    double darkening = 1 - (1 - distance(myState, new StateType()))*4./5.;
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
  std::string GroupFormModel<Graph>::
  printState(State* s) const
  {
    StateType* myState = dynamic_cast<StateType*>(s);
    std::stringstream ss;
    std::streamsize prec = ss.precision();
    ss << this->getVertexState(myState->getState()) << " <"
       << std::setprecision(2);
    for (unsigned int i = 0; i < myState->getTraits().size(); ++i) {
      if (i > 0) ss << ",";
      ss << myState->getTrait(i);
    }
    ss << std::setprecision(prec) << ">";
    ss.unsetf(std::ios::fixed);
    return ss.str();
  }
  
}

#endif
