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

    typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    EpiSimulator(RandomGenerator& r, Graph& g, unsigned int v = 0);
    virtual ~EpiSimulator() {;}

    virtual bool initialise();

    //! Accessor for the numInfections variable.
    unsigned int getNumInfections() const { return numInfections; }
    //! Accessor for the numInformations variable.
    unsigned int getNumInformations() const { return numInformations; }

    virtual bool parse_options(const po::variables_map& vm);
    virtual bool stopCondition() const;
    virtual void updateEventStats(State* before, State* after,
                                  vertex_descriptor v, vertex_descriptor nb);
    
  private:
  
    unsigned int numInfections; //!< A counter for the number of infections.
    unsigned int numInformations; //!< A counter for the number of informations.

    unsigned int numInfected;
    unsigned int numInformed;

    unsigned int stopInfections;
    unsigned int stopInformations;
    unsigned int infLimit;

    bool stopOutbreak;
    bool stopInfoOutbreak;

    std::vector<std::vector<unsigned int> > generations;
    std::vector<unsigned int> genInf;

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
      ("stats", 
       "write cumulative stats at the end of a run")
      ("r0", po::value<unsigned int>(),
       "calculate r0 for arg generations")
      ;
  }

  template <typename RandomGenerator, typename Graph>
  bool EpiSimulator<RandomGenerator, Graph>::
  initialise()
  {
    bool result = true;

    result &= GillespieSimulator<RandomGenerator, Graph>::initialise();

    const EpiModel_base<Graph>* model =
      dynamic_cast<const EpiModel_base<Graph>*>(this->getModel());

    numInfections = numInformations = numInformed = numInfected = 0;

    for (std::vector<std::vector<unsigned int> >::iterator it = generations.begin();
         it != generations.end(); it++) {
      it->clear();
    }
    for (std::vector<unsigned int>::iterator it = genInf.begin();
         it != genInf.end(); it++) {
      (*it) = 0;
    }
    typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(this->getGraph()); vi != vi_end; ++vi) {
      if (model->isInformed(this->getGraph()[*vi].state)) ++numInformed;
      if (model->isInfected(this->getGraph()[*vi].state)) {
        ++numInfected;
        if (generations.size() > 0) generations[0].push_back(*vi);
      }
    }

    return result;
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
          (*this, this->getVerbose()), 0.));
    }
    if (vm.count("r0")) {
      generations.resize(vm["r0"].as<unsigned int>());
      genInf.resize(vm["r0"].as<unsigned int>(), 0);
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_r0<Graph, EpiSimulator<RandomGenerator, Graph> >
          (generations, genInf, this->getVerbose()), 0.));
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
  updateEventStats(State* before, State* after,
                   vertex_descriptor v, vertex_descriptor nb)
  {
    const EpiModel_base<Graph>* model =
      dynamic_cast<const EpiModel_base<Graph>*>(this->getModel());
  
    if (model) {
      if (model->isInfection(before, after)) {
        ++numInfections;
        if (generations.size() > 0) {
          // see if originator is already in a generation
          unsigned int i = 0;
          std::vector<unsigned int>::iterator res =
            std::find(generations[i].begin(), generations[i].end(), nb);
          while (i < generations.size() && res == generations[i].end()) {
            ++i;
            if (i < generations.size()) {
              res = std::find(generations[i].begin(), generations[i].end(), nb);
            }
          }
          if (i < generations.size()) {
            if (i < (generations.size() - 1)) {
              generations[i+1].push_back(v);
              if (this->getVerbose() >= 2) {
                std::cout << "R0: vertex " << v << " recorded as an infected"
                          << " of generation " << i+1 << std::endl;
              }
              ++numInfected;
            } // else this->getGraph()[nb].state->setState(2);
            ++genInf[i];
          } // else {
//	    this->getGraph()[nb].state->setState(2);
//	  }
        } else {
          ++numInfected;
        }
      }
      if (model->isInformation(before, after)) {
        ++numInformations;
        ++numInformed;
      }
      if (model->isRecovery(before, after)) {
        if (generations.size() > 0) {
          unsigned int i = 0;
          std::vector<unsigned int>::iterator res =
            std::find(generations[i].begin(), generations[i].end(), nb);
          while (i < generations.size() && res == generations[i].end()) {
            ++i;
            if (i < generations.size()) {
              res = std::find(generations[i].begin(), generations[i].end(), nb);
            }
          }
          if (i < generations.size()) --numInfected;
        } else {
          --numInfected;
        }
      }
      if (model->isForgetting(before, after)) --numInformed;
    }
  }

  //! Check if conditions to stop a run have been reached. 
  template <typename RandomGenerator, typename Graph>
  bool EpiSimulator<RandomGenerator, Graph>::
  stopCondition() const
  {

    bool ret = GillespieSimulator<RandomGenerator, Graph>::stopCondition();

    return (ret ||
            ((stopInfections > 0) && (numInfections+1 > stopInfections)) ||
            ((stopOutbreak) && (numInfected == 0)) ||
            ((stopInformations > 0) && (numInformed+1 > stopInformations)) ||
            ((stopInfoOutbreak) && (numInformed == 0)) || 
            ((infLimit > 0) && (numInformed+1 > infLimit)));
  }
}
#endif
