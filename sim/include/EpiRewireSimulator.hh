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

    vertex_descriptor* random_walk
    (vertex_descriptor original_node,  vertex_descriptor source_node,
     std::vector<vertex_descriptor>* previous_nodes = 0,
     unsigned int rewireType = 0, unsigned int baseType = 0);
    
  private:

    uniform_gen rewireGen; //!< The random generator for rewiring

    double rewireProb;
    bool infoRewire;
    bool localRewire;

  };

  template <typename RandomGenerator, typename Graph>
  EpiRewireSimulator<RandomGenerator, Graph>::
  EpiRewireSimulator(RandomGenerator& r, Graph& g, unsigned int v) :
    GillespieSimulator<RandomGenerator, Graph>(r, g, v),
    EpiSimulator<RandomGenerator, Graph>(r, g, v),
    rewireGen(r, boost::uniform_real<> (0,1)),
    infoRewire(false), localRewire(false)
  {
    this->simulator_options.add_options()
      ("dynamic-rewiring", po::value<double>(),
       "base probability of random rewiring after infection")
      ("rewire-local",
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
    
    if (vm.count("dynamic-rewiring")) {
      rewireProb = vm["dynamic-rewiring"].as<double>();
      if (rewireProb > 1.) rewireProb = 1.;
      else if (rewireProb < 0.) rewireProb = 0.;
    }
    if (vm.count("rewire-local")) {
      if (rewireProb == 0.) {
        std::cerr << "WARNING: local rewiring makes sense only if rewiring "
                  << "with probability >0, ignoring option" << std::endl;
      }
      localRewire = true;
    }
    if (vm.count("rewire-info")) {
      if (rewireProb == 0.) {
        std::cerr << "WARNING: info-based rewiring makes sense only if rewiring"
                  << " with probability >0, ignoring option" << std::endl;
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
        if (rewireProb > 0 &&
            rewireGen() < rewireProb * state->getInfo()) {
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
            vertex_descriptor targetVertex;
            if (localRewire) {
              // local rewiring (through random walk)
              unsigned int baseType = infoRewire ? 1 : 0;
              vertex_descriptor* findVertex =
                random_walk(v, v, 0, 0, baseType);
              if (findVertex) {
                foundNew = true;
                targetVertex = *findVertex;
              } 
            } else {
              // random rewiring
              do {
                targetVertex =
                  static_cast<unsigned int>(rewireGen() * num_vertices(g));
                if (targetVertex != v) {
                  std::pair<edge_descriptor, bool> e = 
                    edge(v, targetVertex, g);
                  if (e.second) {
                    foundNew = (g[e.first].type != 0 &&
                                !g[e.first].parallel);
                  } else {
                    foundNew = true;
                  }
                }
              } while (!foundNew);
            }
            if (foundNew) {
              if (this->getVerbose() >=2) {
                std::cout << "Rewiring " << del_edge << " to ("
                          << v << "," << targetVertex << ")" << std::endl;
              }
              boost::remove_edge(del_edge, g);
              boost::add_edge(v, targetVertex, g_edge_property_type(0), g);
            }
          }
        } 
      }
    }
    EpiSimulator<RandomGenerator, Graph>::updateEventStats
      (before, after, v, nb);
  }
  
  template <typename RandomGenerator, typename Graph>
  typename EpiRewireSimulator<RandomGenerator, Graph>::vertex_descriptor*
  EpiRewireSimulator<RandomGenerator, Graph>::random_walk
  (vertex_descriptor original_node, vertex_descriptor source_node,
   std::vector<vertex_descriptor>* previous_nodes,
   unsigned int rewireType, unsigned int baseType)
  {
    Graph& graph = this->getGraph();

    bool delPrevious = false;
    if (!previous_nodes) {
      previous_nodes = new std::vector<vertex_descriptor>;
      delPrevious = true;
    }

    if (this->getVerbose() >=2) {
      std::cout << "random_walk: original_node "
                << original_node << ", previous_nodes";
      if (previous_nodes->size() > 0) {
        for (unsigned int i = 0;
             i < previous_nodes->size(); ++i) {
          std::cout << " " << (*previous_nodes)[i];
        }
      } else {
        std::cout << " none";
      }
      std::cout << ", source_node "
                << source_node << std::endl;
    }
    
    vertex_descriptor* target_node = 0;
    std::vector<vertex_descriptor> possible_nodes;
    std::vector<vertex_descriptor> all_neighbours;
    out_edge_iterator oi, oi_end;;

    for (tie(oi, oi_end) =
           boost::out_edges(source_node, graph);
         oi != oi_end; ++oi) {
      if (graph[*oi].type == baseType) {
        bool previous = false;
        for (unsigned int i = 0;
             i < previous_nodes->size() && !previous; ++i) {
          if (target(*oi, graph) == (*previous_nodes)[i]) {
            previous = true;
          }
        }
        if (!previous) {
          std::pair<edge_descriptor, bool> e = 
            edge(original_node, target(*oi, graph), graph);
          if (!e.second ||
              (graph[e.first].type != rewireType &&
               !graph[e.first].parallel)) {
            possible_nodes.push_back(target(*oi, graph));
            if (this->getVerbose() >=2) {
              std::cout << "considering neighbour " << target(*oi, graph)
                        << std::endl;
            }
          } else {
            previous_nodes->push_back(target(*oi, graph));
          }
          all_neighbours.push_back(target(*oi, graph));
        }
      }
    }

    if (possible_nodes.size() > 0) {
      target_node = new vertex_descriptor;
      *target_node =
        possible_nodes[static_cast<unsigned int>
                       (rewireGen() * possible_nodes.size())];
    } else if (all_neighbours.size() > 0) {
      vertex_descriptor step =
        static_cast<unsigned int>(rewireGen() * all_neighbours.size());
      previous_nodes->push_back(source_node);
      target_node = random_walk(original_node, step, previous_nodes,
                                rewireType, baseType);
    } else {
      if (this->getVerbose() >=2) {
        std::cout << "no nodes available" << std::endl;
      }
    }
    
    if (this->getVerbose() >=2) {
      std::cout << "returning ";
      if (target_node) {
        std::cout << *target_node;
      } else {
        std::cout << "nil";
      }
      std::cout << std::endl;
    }
    if (delPrevious) { delete previous_nodes; }
    return target_node;
  }    
}

#endif
