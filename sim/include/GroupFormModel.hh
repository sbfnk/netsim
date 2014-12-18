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
  
  GroupFormState(unsigned int b = 0)
    : State(b) {;}
  ~GroupFormState() {;}
  
  virtual State* clone() const { return new GroupFormState(*this); }
  
private:
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
  
  private:
  };

  template <class Graph>
  GroupFormModel<Graph>::GroupFormModel(unsigned int)
  {
    this->edgeTypes.push_back(Label("x", "", 0, "style=\"solid\""));
    this->model_options.add_options()
      ("nstates", po::value<unsigned int>()->default_value(1),
       "number of initial states")
      ;
  }

  template <class Graph>
  void GroupFormModel<Graph>::Init
  (const po::variables_map& vm, std::vector<StatRecorder<Graph>*>& rec)
  {
    Model<Graph>::Init(vm, rec);

    // ungrouped state
    this->vertexStates.push_back
      (Label("0", "01;37", 0, "", Label::rgbColour(0,0,0)));
    
    for (unsigned int i = 1; i < vm["nstates"].as<unsigned int>(); ++i) {
      addState();
    }
  }

  template <class Graph>
  std::vector<unsigned int> GroupFormModel<Graph>::
  getRGB(State* s) const
  {
    StateType* myState = dynamic_cast<StateType*>(s);
    std::vector<unsigned int> rgb;
    
    rgb.push_back
      (static_cast<unsigned int>
       (this->getVertexState(myState->getState()).getRGB(0)));
    rgb.push_back
      (static_cast<unsigned int>
       (this->getVertexState(myState->getState()).getRGB(1)));
    rgb.push_back
      (static_cast<unsigned int>
       (this->getVertexState(myState->getState()).getRGB(2)));
    
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
    ss << std::setprecision(prec) << ">";
    ss.unsetf(std::ios::fixed);
    return ss.str();
  }
  
}

#endif
