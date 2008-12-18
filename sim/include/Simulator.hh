/*! \file Simulator.hh
  \brief The Simulator class.
*/
#ifndef SIMULATOR_HH
#define SIMULATOR_HH

#include <vector>

#include <boost/program_options.hpp>

#include "Model.hh"

namespace po = boost::program_options;

//! \addtogroup epi_simulator Epidemic simulator
//! \addtogroup simulator Simulator

/*! \brief Base class for epidemic simulation classes
  
Classes derived from this class implement a given simulation algorithm, e.g. the
Gillespie algorithm or, possibly, other algorithms with different waiting times
etc. Most importantly, this class provides a declaration for the updateState and
getTime methods used by the main simulation code.

\ingroup simulator
\ingroup epi_simulator
*/

template <typename Graph>
class Simulator
{

public:

  /*! \brief Constructor.
  \param[in] v verbose intialiser
  */
  Simulator(Graph& g, unsigned int v = 0) :
    graph(g), model(0), verbose(v), time(0.)
  {
    po::options_description new_options1("\nSimulator options");
    new_options1.add_options()
      ("model,m", po::value<std::string>(),
       "model to use")
        ;
    po::options_description new_options2("\nStats options");
    new_options2.add_options()
      ("data,d", po::value<double>(), 
       "write output data at arg timesteps")
      ("pairs",
       "count pairs")
      ("triples",
       "count triples")
      ("graphviz,g", po::value<double>(),
       "create graphviz output in the images directory at arg timesteps")
      ("lattice,l", po::value<double>(),
       "create pixelised lattice output at arg timesteps")
      ;
    po::options_description new_options3("\nStop condition options");
    new_options3.add_options()
    ("tmax", po::value<double>()->default_value(0.),
     "time after which to stop\n(use tmax=0 to run until extinction)")
      ;
    simulator_options.add(new_options1);
    recorder_options.add(new_options2);
    stop_options.add(new_options3);
  }

  //! Destructor.
  virtual ~Simulator() { if (model) delete model;}

  //! Initialise the simulation. Implemented by derived classes.
  virtual void initialise()
  {
    time = 0.;
    for (unsigned int i = 0; i < statRecorders.size(); ++i) {
      statRecorders[i]->reset(dir);
      statRecorders[i]->update(graph, time, true);
    }
  }
  
  //! Perform an update. Implemented by derived classes.
  virtual bool updateState() = 0;

  //! Print the state of the simulation. Implemented by derived classes.
  virtual void print() {;}

  //! Accessor for the time variable.
  double getTime() const { return time; };
  //! Accessor for the verbose variable.
  unsigned int getVerbose() const { return verbose; };

  /*! \brief Update the current time in the simulation
  \param[in] t The timestep to advance the time by.
  */
  void updateTime(double t) { time += t; };

  //! Check if conditions to stop a run have been reached. 
  virtual bool stopCondition() const
  { return ((stopTime > 0 && getTime() >= stopTime)); }

  virtual bool parse_options(const po::variables_map& vm)
  {
    if (vm.count("model")) {
      std::string modelString = vm["model"].as<std::string>();
      while (knownModels.size() > 0) {
        if (modelString == knownModels[0].first) {
          model = knownModels[0].second;
        } else {
          delete knownModels[0].second;
        }
        knownModels.erase(knownModels.begin());
      }
      if (!model) {
	std::cerr << "ERROR: unknown model: " << modelString << std::endl;
	return false;
      }
    } else {
      if (knownModels.size() > 0) {
        std::cerr << "WARNING: no model specified, using first available: "
                  << knownModels[0].first << std::endl;
      }
      model = knownModels[0].second;
    }

    stopTime = vm["tmax"].as<double>();
    bool pairs = (vm.count("pairs"));
    bool triples = (vm.count("triples"));
    if (model) {
      if (vm.count("data")) {
        statRecorders.push_back(new StatRecorder<Graph>
                                (new write_sim_data<Graph, Model<Graph> >
                                 (*model, pairs, triples),
                                 vm["data"].as<double>()));
      }
      if (vm.count("graphviz")) {
        statRecorders.push_back(new StatRecorder<Graph>
                                (new write_sim_graph<Graph, Model<Graph> >(*model),
                                 vm["graphviz"].as<double>()));
      }
      if (vm.count("lattice")) {
        statRecorders.push_back(new StatRecorder<Graph>
                                (new write_sim_lattice<Graph, Model<Graph> >(*model),
                                 vm["lattice"].as<double>()));
      }
      if (verbose>=1) {
        double freq;
        if (verbose >=2) {
          freq = -1;
        } else {
          freq = (vm.count("data") ? vm["data"].as<double>() : 0.);
        }
        statRecorders.push_back(new StatRecorder<Graph>
                                (new print_sim_status<Graph, Model<Graph> >
                                 (*model, pairs, triples),
                                 freq));
      }
    }
    
    return true;
  }

  void InitModel(const po::variables_map& vm)
  {
    if (model) {
      model->Init(vm, statRecorders);
      if (verbose >= 1) model->Print();
    }
  }

  void updateStats()
  {
    bool force = stopCondition();
    for (unsigned int i = 0; i < statRecorders.size(); ++i) {
      statRecorders[i]->update(graph, time, force);
    }
  }

  Graph& getGraph() 
  { return graph; }

  const Graph& getGraph() const
  { return graph; }

  po::options_description getOptions() const
  {
    po::options_description all_options;
    all_options.add(simulator_options).add(recorder_options).add(stop_options);
    return all_options;
  }

  const Model<Graph>* getModel() const
  { return model; }

  bool doIO() const
  { return (statRecorders.size() > 0); }

  void setDir(std::string s)
  { dir = s; }

protected:

  void stopRun()
  { stopTime = time; updateStats(); }
  
  po::options_description simulator_options;
  po::options_description recorder_options;
  po::options_description stop_options;
  std::vector<StatRecorder<Graph>*> statRecorders;
  Graph& graph; //!< The graph determining how vertices can effect another

  std::vector<std::pair<std::string, Model<Graph>*> > knownModels;
  
private:

  Model<Graph>* model;

  unsigned int verbose; //!< The verbosity level.
  double time; //!< The current time of the simulation.
  double stopTime; //!< The time to stop the simulation

  std::string dir;

};

//----------------------------------------------------------
/*! \brief Simulation algorithms.
*/
namespace Simulators {} // just define the namespace

#endif
