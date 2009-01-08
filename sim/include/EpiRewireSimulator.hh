/*! \file EpiRewireSimulator.hh
  \brief The Simulators::EpiRewireSimulator class.
*/

#ifndef EPIREWIRESIMULATOR_HH
#define EPIREWIRESIMULATOR_HH

#include <boost/graph/random.hpp>

#include "EpiSimulator.hh"
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
  class EpiRewireSimulator :
    virtual public EpiSimulator<RandomGenerator, Graph>
  {

  public:

    typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
    edge_descriptor;
    typedef typename boost::edge_property_type<Graph>::type
    g_edge_property_type;
    typedef typename boost::variate_generator
    <RandomGenerator&, boost::uniform_real<> > uniform_gen;
    

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    EpiRewireSimulator(RandomGenerator& r, Graph& g, unsigned int v = 0);
    virtual ~EpiRewireSimulator() {;}

    virtual bool parse_options(const po::variables_map& vm);
    virtual void updateEventStats(State* before, State* after,
                                  vertex_descriptor v, vertex_descriptor nb);
    
  private:

    uniform_gen rewireGen; //!< The random generator for rewiring

    double randRewire;
    double localRewire;
    bool infoRewire;

  };

  template <typename RandomGenerator, typename Graph>
  EpiRewireSimulator<RandomGenerator, Graph>::
  EpiRewireSimulator(RandomGenerator& r, Graph& g, unsigned int v) :
    GillespieSimulator<RandomGenerator, Graph>(r, g, v),
    EpiSimulator<RandomGenerator, Graph>(r, g, v),
    rewireGen(r, boost::uniform_real<> (0,1)),
    infoRewire(false)
  {
    this->simulator_options.add_options()
      ("rewire-random", po::value<double>(),
       "base probability of random rewiring after infection")
      ("rewire-local", po::value<double>(),
       "base probability of local rewiring after infection")
      ("rewire-info",
       "try to rewire to info network")
      ;
  }

  
  template <typename RandomGenerator, typename Graph>
  bool EpiRewireSimulator<RandomGenerator, Graph>::
  parse_options(const po::variables_map& vm)
  {
    bool ret = EpiSimulator<RandomGenerator, Graph>::parse_options(vm);
    
    if (vm.count("rewire-random")) {
      randRewire = vm["rewire-random"].as<double>();
      if (randRewire > 1.) randRewire = 1.;
      else if (randRewire < 0.) randRewire = 0.;
    }
    if (vm.count("rewire-local")) {
      if (randRewire > 0.) {
        std::cerr << "WARNING: cannot rewire randomly and locally, "
                  << " will rewire only randomly" << std::endl;
      } else {
        localRewire = vm["rewire-local"].as<double>();
        if (localRewire > 1.) localRewire = 1.;
        else if (localRewire < 0.) localRewire = 0.;
      }
    }
    if (vm.count("info-based")) {
      if (randRewire + localRewire == 0.) {
        std::cerr << "WARNING: info-based makes sense only if rewiring, "
                  << "ignoring option" << std::endl;
      }
      infoRewire = true;
    }
    return ret;
  }


  template <typename RandomGenerator, typename Graph>
  void EpiRewireSimulator<RandomGenerator, Graph>::
  updateEventStats(State* before, State* after,
                   vertex_descriptor v, vertex_descriptor nb)
  {
    const EpiModel_base<Graph>* model =
      dynamic_cast<const EpiModel_base<Graph>*>(this->getModel());
    const DimInfoState* state =
      dynamic_cast<DimInfoState*>(after);
    Graph& g = this->getGraph();

  
    if (model && state) {
      if (model->isInformation(before, after) && v != nb) {
        if (randRewire > 0 &&
            rewireGen() < randRewire * state->getInfo()) {
          // random rewiring
          edge_descriptor del_edge;
          out_edge_iterator oi, oi_end;
          bool found = false;
          for (tie(oi, oi_end) = out_edges(v, g);
               oi != oi_end && !found; ++oi) {
            if (target(*oi, g) == nb &&
                g[*oi].type == 0) {
              del_edge = *oi;
              found = true;
            }
          }
          if (found) {
            // find somewhere to rewire to
            bool foundNew = false;
            vertex_descriptor randVertex;
            do {
              randVertex =
                static_cast<unsigned int>(rewireGen() * num_vertices(g));
              if (randVertex != v) {
                std::pair<edge_descriptor, bool> e = 
                  edge(v, randVertex, g);
                if (e.second) {
                  foundNew = (g[e.first].type != 0);
                } else {
                  foundNew = true;
                }
              }
            } while (!foundNew);
            if (this->getVerbose() >=2) {
              std::cout << "Rewiring " << del_edge << " to ("
                        << v << "," << randVertex << ")" << std::endl;
            }
            boost::remove_edge(del_edge, g);
            boost::add_edge(v, randVertex, g_edge_property_type(0), g);
          }
        } else if (localRewire > 0 &&
                   rewireGen() < localRewire * state->getInfo()) {
          // local rewiring
        }
      }
    }
    EpiSimulator<RandomGenerator, Graph>::updateEventStats
      (before, after, v, nb);
  }
}
#endif
