/*! \file Model.hh
  \brief The Model class.
*/
#ifndef MODEL_HH
#define MODEL_HH

#include <list>
#include <vector>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>

#include "Vertex.hh"
#include "StatRecorder.hh"
#include "sim_statistics.hh"

//! \addtogroup models Models

namespace po = boost::program_options;

//----------------------------------------------------------
/*!
  \brief A label used for displaying a particular vertex state or edge type.
  \ingroup models
*/
class Label
{

public:

  struct rgbColour
  {
    rgbColour(unsigned int r, unsigned int g, unsigned int b):
      red(r), green(g), blue(b) {;}

    unsigned int red;
    unsigned int green;
    unsigned int blue;

  };

  /*! \brief Constructor.
  \param[in] t text initialiser
  \param[in] c colour initialiser
  \param[in] i id initialiser
  \param[in] d drawOption initialiser
  */
  Label(std::string t, std::string c, unsigned int i,
        std::string drawOption = "",
        rgbColour r = rgbColour(255,255,255)):
    text(t), colour(c), id(i), rgb(3)
  {
    rgb[0] = r.red;
    rgb[1] = r.green;
    rgb[2] = r.blue;
  }

  //! Accessor for the text variable.
  const std::string& getText() const
  { return text; }

  //! Print a coloured letter representing the vertex state.
  void print(std::ostream &os) const
  { os << "\033[" << colour << "m" << text << "\033[0m"; }

  //! Accessor for the rgb variable
  unsigned int getRGB(unsigned int i) const
  { return rgb[i]; }

  //! Accessor for the drawOption variable
  const std::string& getDrawOption() const
  { return drawOption; }

  //! Accessor for the colour variable
  const std::string& getColour() const
  { return colour; }

  //! Accessor for the id variable
  const unsigned int& getId() const
  { return id; }

private:

  std::string text; //!< The letter marking the state of the variable
  std::string colour; //!< The colour marking the state of the variable
  unsigned int id; //!< A unique id assigned to the Label

  std::vector<unsigned int> rgb; //!< The rgb intensities used for graphviz
  std::string drawOption; //!< Any extra graphviz drawOption

};

std::ostream& operator<<(std::ostream& os, const Label& l);

//----------------------------------------------------------
/*! \brief Base class for models to be simulated on a network.
  
Classes derived from this class implement a given model, i.e. the model
parameters and the state transitions that can happen according to the model. 

\ingroup models
*/
template <class Graph>
class Model
{
  
public:

  /*! \brief Constructor.
    
  Here, all states, edge types and their corresponding Labels, as well as the
  model parameters are defined by the classes derived from Model.

  \param[in] v verbose initialiser
  */
  Model(std::string name = "", unsigned int v = 0):
    model_options(po::options_description("\nModel options ("+name+")")),
    verbose(v)
  {;}
  //! Destructor.
  virtual ~Model() = 0;
  
  virtual void Init(const po::variables_map& vm,
                    std::vector<StatRecorder<Graph>*>& rec);
  
  void print(std::ostream &os) const;
  void Print() const;

  virtual Model<Graph>* clone() = 0;

  /*! \brief Get node events.

  This is implemented by classes derived from Model. It considers all node
  events that can happen to a node in a given state, according to the Model
  under consideration. These events are then stored with their corresponding
  rates in an eventList. 
  
  \param[out] events The list of events that can happen to the vertex.
  \param[in] state The state of the vertex.
  \param[in] nb The vertex id.
  \return The sum of the rates of all events that have been stored in the
  event list 
  */
  virtual unsigned int getNodeEvents(eventList& events,
                                     State* state,
                                     unsigned int nb) const
  { return 0; }
  /*! \brief Get edge events.

  This is implemented by classes derived from Model. It considers all edge events
  that can happen to a node in a given state over an edge of given type from a
  neighbour of another given state, according to the Model under
  consideration. These events are then stored with their corresponding rates in
  an eventList. 

  \param[out] events The list of edge events that can happen to the vertex.
  \param[in] state The state of the vertex.
  \param[in] detail The detailed (real-valued) state of the vertex.
  \param[in] edge The type of edge.
  \param[in] nbState The state of the neighbour.
  \param[in] nbDetail The detailed (real-valued) state of the neighbour.
  \param[in] nb The vertex id of the neighbour.
  \return The sum of the rates of all events that have been stored in the
  event list
  */
  virtual unsigned int getEdgeEvents(eventList& events, State* state, 
                                     unsigned int edge, State* nbState,
                                     unsigned int nb) const
  { return 0; }

  virtual State* newState() const
  { return new State; }

  //! Accessor for vertexStates
  const std::vector<Label>& getVertexStates() const
  { return vertexStates; }

  //! Accessor for an element of vertexStates
  const Label& getVertexState(unsigned int id) const
  { return vertexStates[id]; }
  
  //! Accessor for edgeTypes
  const std::vector<Label>& getEdgeTypes() const
  { return edgeTypes; }

  //! Accessor for an element of edgeTypes
  const Label& getEdgeType(unsigned int id) const
  { return edgeTypes[id]; }

  //! Accessor for model_options
  const po::options_description& getOptions() const
  { return model_options; }

  //! Get initial value of state_detail for a given state
  virtual double getInitDetail(unsigned int baseState) const
  { return 0.; }

  unsigned int getVerbose() const
  { return verbose; }

  virtual std::string printState(State* s) const
  { std::stringstream ss; ss << getVertexState(s->getState()); return ss.str();}
  
  virtual std::vector<unsigned int> getRGB(State* s) const
  {
    std::vector<unsigned int> rgb;

    rgb.push_back(this->getVertexStates()[s->getState()].getRGB(0));
    rgb.push_back(this->getVertexStates()[s->getState()].getRGB(1));
    rgb.push_back(this->getVertexStates()[s->getState()].getRGB(2));
    
    return rgb;
  }

  std::vector<double> getColour(State* s) const
  {
    std::vector<unsigned int> rgb = this->getRGB(s);
    std::vector<double> colour;
    for (unsigned int i = 0; i < rgb.size(); ++i) {
      colour.push_back(rgb[i] / 255.);
    }
    return colour;
  }

protected:

  /*! \brief The vertex states.
    
  A vector of all states which a vertex can assume in the model, and
  how to visualise them (as defined by Label).
  */
  std::vector<Label> vertexStates;

  /*! \brief The edge types.
    
  A vector of all types which an edge can be of in the model, and
  how to visualise them (as defined by Label).
  */
  std::vector<Label> edgeTypes;

  //! Command line options the model defines and uses.
  po::options_description model_options;
  
  /*! \brief The model parameters
    
  A map of command line options to the model paramters
  */
  std::map<std::string, double*> params;
  /*! \brief The model parameters
    
  A map of command line options to the model paramters
  */
  std::map<std::string, unsigned int*> intParams;
  /*! \brief The model rates
    
  A map of command line options to the model rates
  */
  std::map<std::string, unsigned int*> rates;

private:
  //! Verbosity level
  unsigned int verbose;
  
};

/*! \brief Initialise model paramters.
   
This should be called after the command line parameters of the model haven been
assigned. It initialises the model parameter variables with the values found in
the command line parameters.

\param[in] vm The map of command line parameters
*/

template<class Graph>
void Model<Graph>::Init(const po::variables_map& vm,
                        std::vector<StatRecorder<Graph>*>& rec)
{
  // loop over all model parameters

  for (std::map<std::string, unsigned int*>::iterator it = rates.begin();
       it != rates.end(); it++) {
    if (vm.count(it->first)) {
      // command line parameter has been specified, assign to rate
      if (vm[it->first].as<double>() > 1e+5) {
        std::cerr << "WARNING: rates bigger than 1e+5 not supported."
                  << std::endl;
        std::cerr << "setting " << it->first << " to 1e+5." << std::endl;
        *(it->second) = static_cast<unsigned int>(1e+9);
      } else {
        *(it->second) =
          static_cast<unsigned int>(vm[it->first].as<double>() * 1e+4 + .5);
      }
    } else {
      std::cerr << "WARNING: no " << it->first << " given" << std::endl;
      std::cerr << "setting to 0." << std::endl;
      *(it->second) = 0;
    }
  }
  for (std::map<std::string, double*>::iterator it = params.begin();
       it != params.end(); it++) {
    if (vm.count(it->first)) {
      // command line parameter has been specified, assign to model variable
      *(it->second) = vm[it->first].as<double>();
    } else {
      std::cerr << "WARNING: no " << it->first << " given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      *(it->second) = 0;
    }
  }
  for (std::map<std::string, unsigned int*>::iterator it = intParams.begin();
       it != intParams.end(); it++) {
    if (vm.count(it->first)) {
      // command line parameter has been specified, assign to model variable
      *(it->second) = vm[it->first].as<unsigned int>();
    } else {
      std::cerr << "WARNING: no " << it->first << " given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      *(it->second) = 0;
    }
  }
}

template<class Graph>
inline Model<Graph>::~Model() {;}

//! Print the model parameters to ostream.
template<class Graph>
void Model<Graph>::print(std::ostream &os) const
{
  if ((params.size() + intParams.size()) > 0) {
    os << "Model parameters:" << std::endl;
    os << "=================" << std::endl;
    for (std::map<std::string, double*>::const_iterator it =
           params.begin(); it != params.end(); it++) {
      os << it->first << ": " << *(it->second) << std::endl; 
    }
    for (std::map<std::string, unsigned int*>::const_iterator it =
           intParams.begin(); it != intParams.end(); it++) {
      os << it->first << ": " << *(it->second) << std::endl; 
    }
    os << std::endl;
  }
  if (rates.size() > 0) {
    os << "Model rates:" << std::endl;
    os << "=================" << std::endl;
    for (std::map<std::string, unsigned int*>::const_iterator it =
           rates.begin(); it != rates.end(); it++) {
      os << it->first << ": " << (*(it->second))/1e+4 << std::endl; 
    }
  }
}

/*! \brief Stream operator for Model.

Stream the model paramters.
\param[in, out] os The stream to write the Model to
\param[in] l The Model to stream

\return A reference to the stream written to
*/
template <class Graph>
std::ostream& operator<<(std::ostream& os, const Model<Graph>& m)
{ m.print(os); return os; }

//! Print the model parameters to the screen.
template<class Graph>
void Model<Graph>::Print() const
{
  std::cout << *this;
}

//----------------------------------------------------------
/*! \brief The models of interaction for usage in the simulation.
  \ingroup models
*/
namespace Models {} // just define the namespace

#endif
