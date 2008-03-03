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
/*! \brief The state a vertex can assume
  
\ingroup models
\ingroup gillespie_simulator
*/
  //! Type denoting the state of an individual
  struct State {

    //! Constructor.
    State() : base(0), detail(0.) {;}

    /*! \brief Constructor.
      \param[in] b base state.
    */
    State(unsigned int b) : base(b), detail(0.) {;}

    /*! \brief Constructor.
      \param[in] b base state.
      \param[in] d refined state.
    */
    State(unsigned int b, double d) : base(b), detail(d) {;}
    /*! Integer variable corresponding to the current state of the individual,
      to be defined by the used Model.
    */
    unsigned int base;
    /*! Real-valued variable refining the current state of the individual, to be
      defined by the used Model.
    */
    double detail;
  };

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
  \param[in] n nb initialiser
  \param[in] e et initialiser
  */
  Event(unsigned int r=0, State s = State(), unsigned int n=0, unsigned int e=0) :
    rate(r), newState(s), nb(n), et(e) {}
  
  //! The rate at which an event occurs (depends on model parameters),
  //  multiplied by 10^5 for integer representation
  unsigned int rate; 
  //! The state an event will change a vertex to
  State newState; 
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
  void print(std::ostream &os) const;
  void Print() const;

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
                                     State state,
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
  virtual unsigned int getEdgeEvents(eventList& events, State state, 
                                     unsigned int edge, State nbState,
                                     unsigned int nb) const = 0;

  /*! \brief Check whether an event is an infection.

  This is implemented by classes derived from Model.

  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
   */
  virtual bool isInfection(State before_state, State after_state) const
  { return false; }

  /*! \brief Check whether an event is a recovery.

  This is implemented by classes derived from Model.

  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
   */
  virtual bool isRecovery(State before_state, State after_state) const
  { return false; }
  
  /*! \brief Check whether an event is an information event
     
  This is implemented by classes derived from Model.
  
  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
  */
  virtual bool isInformation(State before_state, State after_state) const
  { return false; }

  /*! \brief Check whether an event is a forgetting event
     
  This is implemented by classes derived from Model.
  
  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
  */
  virtual bool isForgetting(State before_state, State after_state) const
  { return false; }

  /*! \brief Check whether a vertex is infected
     
  This is implemented by classes derived from Model.
  
  \param[in] state The state to check.
  */
  virtual bool isInfected(State state) const
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
  virtual double getInitDetail(unsigned int baseState) const
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
  /*! \brief The model rates
    
  A map of command line options to the model rates
  */
  std::map<std::string, unsigned int*> rates;

  //! Verbosity level
  unsigned int verbose;
  
};

std::ostream& operator<<(std::ostream& os, const Model& m);

//----------------------------------------------------------
/*! \brief The models of interaction for usage in the simulation.
  \ingroup models
*/
namespace Models {} // just define the namespace

#endif
