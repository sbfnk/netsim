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

//! \addtogroup models Models

namespace po = boost::program_options;

//----------------------------------------------------------
/*! \brief An event which can happen to a vertex.
  
\ingroup models
\ingroup gillespie_simulator
*/
struct Event
{
  /*! \brief Constructor.
  \param[in] r rate initialiser
  \param[in] s newState initialiser
  \param[in] d newDetail initialiser
  \param[in] n nb initialiser
  \param[in] e et initialiser
  */
  Event(double r=0., int s=0, unsigned int n=0, unsigned int e=0,
        double d=0.) :
    rate(r), newState(s), newDetail(d), nb(n), et(e) {}
  
  //! The rate at which an event occurs (depends on model parameters)
  double rate; 
  //! The state an event will change a vertex to
  int newState; 
  //! The refined real value of the state an event will change a vertex to
  double newDetail;
  //! The neighbour "responsible" for the event
  unsigned int nb; 
  //! The edge type over which event is transmitted (if applicable)
  unsigned int et;

};

typedef std::vector<Event> eventList;

//----------------------------------------------------------
/*!
  \brief A label used for displaying a particular vertex state or edge type.
  \ingroup models
*/
class Label
{

public:

  /*! \brief Constructor.
  \param[in] t text initialiser
  \param[in] c colour initialiser
  \param[in] i id initialiser
  \param[in] d drawOption initialiser
  */
  Label(std::string t, std::string c, unsigned int i, std::string d = "") :
    text(t), colour(c), id(i), drawOption(d) {}
  
  //! Accessor for the text variable.
  const std::string& getText() const
  { return text; }

  //! Print a coloured letter representing the vertex state.
  void print(std::ostream &os) const
  { os << "\033[" << colour << "m" << text << "\033[0m"; }

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
  std::string drawOption;   //!< The drawOption used for graphviz output
  
};

std::ostream& operator<<(std::ostream& os, const Label& l);

//----------------------------------------------------------
/*! \brief Base class for models to be simulated on a network.
  
Classes derived from this class implement a given model, i.e. the model
parameters and the state transitions that can happen according to the model. 

\ingroup models
*/
class Model
{
  
public:

  /*! \brief Constructor.
    
  Here, all states, edge types and their corresponding Labels, as well as the
  model parameters are defined by the classes derived from Model.

  \param[in] v verbose initialiser
  */
  Model(unsigned int v = 0):
    model_options(po::options_description("\nModel options:")), verbose(v) {}
  //! Destructor.
  virtual ~Model() {}
  
  virtual void Init(po::variables_map& vm);
  void Print();

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
  virtual double getNodeEvents(eventList& events,
                               unsigned int state,
                               unsigned int nb) const = 0;
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
  virtual double getEdgeEvents(eventList& events,
                               unsigned int state, double detail,
                               unsigned int edge, unsigned int nbState,
                               double nbDetail, unsigned int nb) const = 0;

  /*! \brief Check whether an event is an infection.

  This is implemented by classes derived from Model.

  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
   */
  virtual bool isInfection(unsigned int before_state,
                           unsigned int after_state) const
  { return false; }

  /*! \brief Check whether an event is a recovery.

  This is implemented by classes derived from Model.

  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
   */
  virtual bool isRecovery(unsigned int before_state,
                          unsigned int after_state) const
  { return false; }
  
  /*! \brief Check whether an event is an information event
     
  This is implemented by classes derived from Model.
  
  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
  */
  virtual bool isInformation(unsigned int before_state,
                             unsigned int after_state) const
  { return false; }

  /*! \brief Check whether an event is a forgetting event
     
  This is implemented by classes derived from Model.
  
  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
  */
  virtual bool isForgetting(unsigned int before_state,
                            unsigned int after_state) const
  { return false; }

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
  virtual double getInitDetail(unsigned int state) const
  { return 0.; }

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

  //! Verbosity level
  unsigned int verbose;
  
};

//----------------------------------------------------------
/*! \brief The models of interaction for usage in the simulation.
  \ingroup models
*/
namespace Models {} // just define the namespace

#endif
