/******************************************************************/
// Model.hh
// 
// this contains the classes related to the model:
//   -- event
//        for the events that can happen in the simulation
//   -- Model
//        class implementing the used model, i.e. parameters and
//        dependance of the rates of various processes on node states
//        and edge properties
//
/******************************************************************/
#ifndef MODEL_HH
#define MODEL_HH

#include <list>
#include <vector>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

struct event
{
      double rate; // rate at which an event occurs (depends on
                   // model parameters) 
      int newState; // state an event will change a vertex to
};

typedef std::list<event> eventList;

class Label
{

public:

  Label(std::string t, std::string c, unsigned int i, std::string d = "") :
    text(t), color(c), id(i), drawOption(d) {}
  
  const std::string& getText() const
  { return text; }

  void print(std::ostream &os) const
  { os << "\033[" << color << "m" << text << "\033[0m"; }

  const std::string& getDrawOption() const
  { return drawOption; }
  
  const unsigned int& getId() const
  { return id; }

private:
  
  std::string text;
  std::string color;
  unsigned int id;
  std::string drawOption;
  
};

// streams a letter associated VertexState/EdgeType, useful for easy printout
std::ostream& operator<<(std::ostream& os, const Label& l);

class Model
{
  
public:

  Model(unsigned int v = 0): verbose(v) {}
  virtual ~Model() {}
  
  virtual void Init(po::variables_map& vm);
  void Print();
  
  virtual double getNodeEvents(eventList& events,
                               unsigned int state) const = 0;
  virtual double getEdgeEvents(eventList& events,
                               unsigned int state, unsigned int edge,
                               unsigned int nbState) const = 0;

  const std::vector<Label>& getVertexStates() const
  { return vertexStates; }

  const Label& getVertexState(unsigned int id) const
  { return vertexStates[id]; }
  
  const std::vector<Label>& getEdgeTypes() const
  { return edgeTypes; }

  const po::options_description& getOptions() const
  { return model_options; }

protected:
  
  std::vector<Label> vertexStates;
  std::vector<Label> edgeTypes;

  po::options_description model_options;
  
  std::map<std::string, double*> params;

  unsigned int verbose;
  
};

#endif
