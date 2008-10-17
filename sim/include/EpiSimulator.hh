/*! \file EpiSimulator.hh
  \brief The Simulators::EpiSimulator class.
*/

#ifndef EPISIMULATOR_HH
#define EPISIMULATOR_HH

#include "GillespieSimulator.hh"
#include "EpiModel.hh"

//! \addtogroup gillespie_simulator Gillespie simulator

namespace Simulators {
  
  /*! \brief The Gillespie simulation
    
  This contains the Gillespie simulation, able to run a simulation on a graph
  using a given model, with a sequential update rule according to
  Gillespie (1977). The rates of the events are stored in a Tree for quick
  access, which is generated using the member function initialise. The actual
  updating of the graph is done in updateState.
  
  \ingroup gillespie_simulator
  */
  template <typename RandomGenerator, typename Graph>
  class EpiSimulator :
    virtual public GillespieSimulator<RandomGenerator, Graph>
  {

  public:

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    EpiSimulator(RandomGenerator& r, Graph& g, unsigned int v = 0);
    virtual ~EpiSimulator() {;}

    virtual void initialise();

    //! Accessor for the numInfections variable.
    unsigned int getNumInfections() const { return numInfections; }
    //! Accessor for the numInformations variable.
    unsigned int getNumInformations() const { return numInformations; }

    virtual bool parse_options(const po::variables_map& vm);
    virtual bool stopCondition() const;
    virtual void updateEventStats(State* before, State* after);
    
  private:
  
    unsigned int numInfections; //!< A counter for the number of infections.
    unsigned int numInformations; //!< A counter for the number of informations.

    unsigned int stopInfections;
    unsigned int stopInformations;
    unsigned int infLimit;

    bool stopOutbreak;
    bool stopInfoOutbreak;
  };

  template <typename RandomGenerator, typename Graph>
  EpiSimulator<RandomGenerator, Graph>::
  EpiSimulator(RandomGenerator& r, Graph& g, unsigned int v) :
    GillespieSimulator<RandomGenerator, Graph>(r, g, v)
  {
    this->stop_options.add_options()
      ("imax", po::value<unsigned int>(),
       "number of infections after which to stop (if >0)")
      ("pmax", po::value<unsigned int>(),
       "number of informations after which to stop (if >0)")
      ("stop-disease", 
       "stop when outbreak has ended")
      ("stop-info", 
       "stop when info outbreak has ended")
      ("limit", po::value<unsigned int>(),
       "limit to the number of infectives (stop if reached) in addition to initial infectives")
      ;
    this->recorder_options.add_options()
      ("stats", po::value<double>(),
       "write cumulative stats at arg timesteps")
      ;
  }

  template <typename RandomGenerator, typename Graph>
  void EpiSimulator<RandomGenerator, Graph>::
  initialise()
  {
    GillespieSimulator<RandomGenerator, Graph>::initialise();
    numInfections = numInformations = 0;
  }
  
  template <typename RandomGenerator, typename Graph>
  bool EpiSimulator<RandomGenerator, Graph>::
  parse_options(const po::variables_map& vm)
  {
    bool ret = GillespieSimulator<RandomGenerator, Graph>::parse_options(vm);
    if (vm.count("stats")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_epi_stats<Graph, EpiSimulator<RandomGenerator, Graph> >
          (*this, this->getVerbose()), 
          vm["stats"].as<double>()));
    }
    if (vm.count("imax")) {
      stopInfections = vm["imax"].as<unsigned int>();
    } else {
      stopInfections = 0;
    }
    if (vm.count("pmax")) {
      stopInformations = vm["pmax"].as<unsigned int>();
    } else {
      stopInformations = 0;
    }
    if (vm.count("limit")) {
      infLimit = vm["limit"].as<unsigned int>();
    } else {
      infLimit = 0;
    }
    stopOutbreak = vm.count("stop-disease");
    stopInfoOutbreak = vm.count("stop-info");
    return ret;
  }


  template <typename RandomGenerator, typename Graph>
  void EpiSimulator<RandomGenerator, Graph>::
  updateEventStats(State* before, State* after)
  {
    const EpiModel_base<Graph>* model =
      dynamic_cast<const EpiModel_base<Graph>*>(this->getModel());
  
    if (model) {
      if (model->isInfection(before, after)) ++numInfections;
      if (model->isInformation(before, after)) ++numInformations;
    }
  }

  //! Check if conditions to stop a run have been reached. 
  template <typename RandomGenerator, typename Graph>
  bool EpiSimulator<RandomGenerator, Graph>::
  stopCondition() const
  {
    const EpiModel_base<Graph>* model =
      dynamic_cast<const EpiModel_base<Graph>*>(this->getModel());
    
    bool ret = Simulator<Graph>::stopCondition();
    unsigned int numInfected = 0;
    unsigned int numInformed = 0;
    if (infLimit>0 || stopOutbreak || stopInfoOutbreak) {
      typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
      vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(this->getGraph()); vi != vi_end; ++vi) {
        if (model->isInformed(this->getGraph()[*vi].state)) ++numInformed;
        if (model->isInfected(this->getGraph()[*vi].state)) ++numInfected;
      }
    }
    return (ret ||
            ((stopInfections > 0) && (numInfections+1 > stopInfections)) ||
            ((stopOutbreak) && (numInfected == 0)) ||
            ((stopInformations > 0) && (numInformed+1 > stopInformations)) ||
            ((stopInfoOutbreak) && (numInformed == 0)) || 
            ((infLimit > 0) && (numInformed+1 > infLimit)));
  }
}
#endif
