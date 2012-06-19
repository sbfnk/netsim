/*! \file EpiModel.hh
  \brief The Model class.
*/
#ifndef EPIMODEL_HH
#define EPIMODEL_HH

#include <list>
#include <vector>
#include <string>
#include <iostream>

#include <Model.hh>

//! \addtogroup models Models

namespace po = boost::program_options;

template <class Graph>
class EpiModel_base
  : public Model<Graph>
{
public:
  EpiModel_base(std::string n = "", unsigned int v = 0) : Model<Graph>(n, v)
  {;}
  virtual ~EpiModel_base() = 0;
  
  virtual bool isInfection(State* before_state, State* after_state) const
  { return false; }

  /*! \brief Check whether an event is a recovery.

  This is implemented by classes derived from Model.

  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
   */
  virtual bool isRecovery(State* before_state, State* after_state) const
  { return false; }
  
  /*! \brief Check whether an event is an information event
     
  This is implemented by classes derived from Model.
  
  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
  */
  virtual bool isInformation(State* before_state, State* after_state) const
  { return false; }

  /*! \brief Check whether an event is a forgetting event
     
  This is implemented by classes derived from Model.
  
  \param[in] before_state The state before the event happens.
  \param[in] after_state The state after the event happens.
  */
  virtual bool isForgetting(State* before_state, State* after_state) const
  { return false; }

  /*! \brief Check whether a vertex is infected
     
  This is implemented by classes derived from Model.
  
  \param[in] state The state to check.
  */
  virtual bool isInfected(State* state) const
  { return false; }

  /*! \brief Check whether a vertex is informed
     
  This is implemented by classes derived from Model.
  
  \param[in] state The state to check.
  */
  virtual bool isInformed(State* state) const
  { return false; }

};

//----------------------------------------------------------
/*! \brief Base class for models to be simulated on a network.
  
Classes derived from this class implement a given model, i.e. the model
parameters and the state transitions that can happen according to the model. 

\ingroup models
*/
template <class ModelState, class Graph>
class EpiModel
  : public EpiModel_base<Graph>
{
  
public:

  typedef ModelState StateType;
  /*! \brief Constructor.
    
  Here, all states, edge types and their corresponding Labels, as well as the
  model parameters are defined by the classes derived from Model.

  \param[in] n name initialiser
  \param[in] v verbose initialiser
  */
  EpiModel(std::string n = "", unsigned int v = 0) : EpiModel_base<Graph>(n, v)
  {;}
  //! Destructor.
  virtual ~EpiModel() = 0;
  
  virtual State* newState() const
  { return new ModelState; }

  State* newState(ModelState s) const
  { ModelState* new_state = new ModelState; *new_state = s; return new_state; }

  virtual State* newState(unsigned int red, unsigned int green,
                          unsigned int blue) const
  {
    return 
      new ModelState
      (this->getStateFromColour((red > 0)*1 + (green > 0)*2 + (blue > 0) * 4));
  } 
};

template <class Graph>
inline EpiModel_base<Graph>::~EpiModel_base() {;}

template <class ModelState, class Graph>
inline EpiModel<ModelState, Graph>::~EpiModel() {;}
//----------------------------------------------------------
/*! \brief The models of interaction for usage in the simulation.
  \ingroup models
*/
namespace Models {} // just define the namespace

#endif
