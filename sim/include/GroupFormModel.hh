/*! \file DimInfoSIRS.hh
  \brief The Models::GroupFormModel class.
*/
#ifndef GROUPFORMMODEL_HH
#define GROUPFORMMODEL_HH

#include "Model.hh"

class GroupFormState :
  public State
{
public:
  
  GroupFormState(unsigned int b = 0, std::vector<double> t = std::vector<double>(0))
    : State(b), trait_vector(t) {;}
  ~GroupFormState() {;}
  
  virtual State* clone() { return new GroupFormState(*this); }
  
  const std::vector<double>& getTraits() const { return trait_vector; }
  const double getTrait(unsigned int i) const
  {
    if (i+1 > trait_vector.size()) {
      return 0.;
    } else {
      return trait_vector[i];
    }
  }
  void setTrait(std::vector<double> t) { trait_vector = t; }

private:
  std::vector<double> trait_vector;
};

namespace Models {
  
  /*! \brief Class implementing the group formation model.

  */
  template <class Graph>
  class GroupFormModel :
    public Model<Graph>
  {

  public:

    GroupFormModel(unsigned int v = 0);
    ~GroupFormModel() {;}
    
    virtual Model<Graph>* clone() { return new GroupFormModel<Graph>(*this); }

    virtual GroupFormState* newState() const
    { return new GroupFormState(); }
    
    virtual void Init(const po::variables_map& vm,
                      std::vector<StatRecorder<Graph>*>& rec);
  
    virtual std::string printState(State* s) const;
    std::vector<unsigned int> getRGB(State* s) const;

    unsigned int getStates() const
    { return states; }

    unsigned int getTraitDim() const
    { return trait_dim; }

    double getAcceptance() const
    { return acceptance; }

    double accept(State* s1, State* s2) const
    {
      GroupFormState* state1 = dynamic_cast<GroupFormState*>(s1);
      GroupFormState* state2 = dynamic_cast<GroupFormState*>(s2);
      return acceptance*(1-distance(state1, state2));
    }

    double distance(State* s1, State* s2) const
    {
      GroupFormState* state1 = dynamic_cast<GroupFormState*>(s1);
      GroupFormState* state2 = dynamic_cast<GroupFormState*>(s2);
      return distance(state1, state2);
    }

  private:

    double distance(GroupFormState* s1, GroupFormState* s2) const
    {
      double dist = 0.;
      for (unsigned int i = 0; i < s1->getTraits().size(); ++i) {
        dist += pow(s1->getTrait(i)-s2->getTrait(i), 2);
      }
      return sqrt(dist);
    }

    unsigned int states; //!< number of states (groups).
    unsigned int trait_dim; //!< Dimension of the trait vector
    double acceptance; //!< Base property of acceptance.
  };

  template <class Graph>
  GroupFormModel<Graph>::GroupFormModel(unsigned int)
  {
    this->model_options.add_options()
      ("states", po::value<unsigned int>()->default_value(1),
       "number of possible states (groups?)")
      ("traitdim", po::value<unsigned int>()->default_value(1),
       "dimensionality of the trait vector")
      ("alpha", po::value<double>()->default_value(1.),
       "alpha")
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
    
    /*************************************/
    // define vertex classes
    /************************************/
    // ungrouped state
    this->vertexStates.push_back
      (Label("0", "01;37", 0, "", Label::rgbColour(0,0,0)));
    for (unsigned int i = 1; i < states+1; ++i) {
      std::stringstream s;
      s << i;
      switch (i) {
       case 1: 
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;32", i, "", Label::rgbColour(0,0,255)));
        break;
       case 2: 
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;31", i, "", Label::rgbColour(255,0,0)));
        break;
       case 3:
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;34", i, "", Label::rgbColour(0,255,0)));
        break;
       default:
        this->vertexStates.push_back
          (Label(s.str().c_str(), "01;37", i, "", Label::rgbColour(0,0,0)));
        break;
      }
    }
  }

  template <class Graph>
  std::vector<unsigned int> GroupFormModel<Graph>::
  getRGB(State* s) const
  {
    GroupFormState* myState = dynamic_cast<GroupFormState*>(s);
    std::vector<unsigned int> rgb;
    
    double darkening = 1 - (1 - distance(myState, new GroupFormState()))*4./5.;
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
    GroupFormState* myState = dynamic_cast<GroupFormState*>(s);
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
