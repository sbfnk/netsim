/*! \file RewireSimulator.hh
  \brief The Simulators::RewireSimulator class.
*/

#ifndef REWIRESIMULATOR_HH
#define REWIRESIMULATOR_HH

#include <math.h>
#include <algorithm>

#include "Simulator.hh"
#include "Vertex.hh"
#include "GroupFormModel.hh"

#include "network/include/community_structure.hh"

#include <boost/graph/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//! \addtogroup simulator Simulator

namespace Simulators {
  
  /*! \brief The Group formation simulation
    
  \ingroup simulator
  */
  template <typename RandomGenerator, typename Graph>
  class RewireSimulator :
    virtual public Simulator<Graph>
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
    vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
    edge_descriptor;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
    typedef typename boost::graph_traits<Graph>::adjacency_iterator
    adjacency_iterator;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type
    vertex_index_type;
    typedef typename boost::vertex_property_type<Graph>::type::value_type
    vertex_property_type;
    typedef typename boost::variate_generator
    <RandomGenerator&, boost::uniform_real<> > uniform_gen;

  public:

    /*! \brief Constructor
    
    \param[in] r randGen initialiser
    \param[in] g graph initialiser
    \param[in] m model initialiser (in Simulator::Simulator)
    \param[in] v verbose initialiser (in Simulator::Simulator)
    */
    RewireSimulator(RandomGenerator& r, Graph& g,
                       unsigned int v = 0);
    ~RewireSimulator() {;}
      
    void initialise();
    bool updateState();
  
    vertex_descriptor* random_state_walk
    (vertex_descriptor original_node,  vertex_descriptor source_node,
     std::vector<vertex_descriptor>* previous_nodes = 0,
     double distanceSum = 0.);
   
   void print();

    bool parse_options(const po::variables_map& vm);

    //! Check if conditions to stop a run have been reached. 
    bool stopCondition() const
    {
      bool ret = Simulator<Graph>::stopCondition();
      unsigned int largest_component = stopComponent;
      if (stopComponent > 0) {
        // find component dist writer if available
        bool found = false;
        for (unsigned int i = 0;
             ((!false) && i < (this->statRecorders.size())); ++i) {
          const write_comp_dist<Graph>* w =
            dynamic_cast<const write_comp_dist<Graph>*>
            (this->statRecorders[i]->getStatFunc());
          if (w) {
            largest_component = w->getLargest();
            found = true;
          }
        }
        if (!found) largest_component = write_component_dist(this->getGraph());
      }
      return (ret || (largest_component < stopComponent));
    }

  private:

    uniform_gen randGen; //!< The random generator to be used for choosing events

    std::vector<vertex_descriptor> active; //!< Active nodes

    unsigned int verbose;
    unsigned int stopComponent;

    double rewireProb;
    double updateProb;
    double randomiseProb;
  
  };

  template <typename RandomGenerator, typename Graph>
  RewireSimulator<RandomGenerator, Graph>::
  RewireSimulator(RandomGenerator& r, Graph& g,
		     unsigned int v) :
    Simulator<Graph>(g, v), randGen(r, boost::uniform_real<> (0,1)),
    active(), verbose(v)
  {
    this->simulator_options.add_options()
      ("rewire-prob,p",po::value<double>()->default_value(0.),
       "probability of rewiring")
      ("state-update-prob,q",po::value<double>()->default_value(0.),
       "probability of state update")
      ("state-randomise-prob,r",po::value<double>()->default_value(0.),
       "probability of random state assignment")
      ;
    this->recorder_options.add_options()
      ("component-dist,t",po::value<double>(),
       "write component distribution in comp directory at arg timesteps")
      ("modularity",po::value<double>(),
       "write modularity at arg timesteps")
      ("community",po::value<double>(),
       "write community distribution in comm directory at arg timesteps")
      ;
    this->stop_options.add_options()
      ("cmin", po::value<unsigned int>()->default_value(0),
       "limit to the size of the largest component (stop if drops to that value")
      ;
    this->knownModels.push_back
      (std::make_pair("GroupForm",
                      new Models::GroupFormModel<Graph>
                      (this->getVerbose())));
  }
  

  template <typename RandomGenerator, typename Graph>
  bool RewireSimulator<RandomGenerator, Graph>::
  parse_options(const po::variables_map& vm)
  { 
    bool ret = Simulator<Graph>::parse_options(vm);
    stopComponent = vm["cmin"].as<unsigned int>();
    rewireProb = vm["rewire-prob"].as<double>();
    if (rewireProb == 0.) {
      std::cerr << "WARNING: rewiring probability 0" << std::endl;
    }
    updateProb = vm["state-update-prob"].as<double>();
    if (rewireProb == 0.) {
      std::cerr << "WARNING: state update probability 0" << std::endl;
    }
    randomiseProb = vm["state-randomise-prob"].as<double>();
    if (rewireProb == 0.) {
      std::cerr << "WARNING: state randomise probability 0" << std::endl;
    }
    if (vm.count("component-dist")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_comp_dist<Graph>(this->getGraph()),
          vm["component-dist"].as<double>()));
    }
    if (vm.count("modularity") || vm.count("community")) {
      double rate = -1;
      if (vm.count("modularity")) {
        if (vm.count("community") &&
            vm["modularity"].as<double>() != vm["community"].as<double>()) {
          rate = std::min(vm["modularity"].as<double>(),
                          vm["community"].as<double>());
          std::cerr << "WARNING: different rates given for community and "
                    << "modularity, using smaller rate " << rate << std::endl;
        } else {
          rate = vm["modularity"].as<double>();
        }
      } else {
        rate = vm["community"].as<double>();
      }
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_community_structure<Graph, Model<Graph> >
          (this->getGraph(), *(this->getModel()),
           vm.count("community"), vm.count("modularity")),
          rate));
    }

    return ret;
  }


  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
    
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  void RewireSimulator<RandomGenerator, Graph>::initialise()
  {
    Simulator<Graph>::initialise();
    Graph& graph = this->getGraph();
    const Models::GroupFormModel<Graph>* model =
      dynamic_cast<const Models::GroupFormModel<Graph>*>(this->getModel());
    
    // set trait for each vertex 
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      std::vector<double> randTraits;
      for (unsigned int i = 0; i < model->getTraitDim(); ++i) {
	randTraits.push_back((randGen)());
      }
      GroupFormState* myState = dynamic_cast<GroupFormState*>(graph[*vi].state);
      myState->setTrait(randTraits);
      if (verbose >= 2) {
        std::cout << "Associating trait (";
	for (unsigned int i = 0; i < model->getTraitDim(); ++i) {
	  if (i>0) std::cout << ",";
	  std::cout << myState->getTrait(i);
	}
	std::cout  << ") with vertex " << *vi << std::endl;
      }
    }
  }
  
  template <typename RandomGenerator, typename Graph>
  typename RewireSimulator<RandomGenerator, Graph>::vertex_descriptor*
  RewireSimulator<RandomGenerator, Graph>::random_state_walk
  (vertex_descriptor original_node, vertex_descriptor source_node,
   std::vector<vertex_descriptor>* previous_nodes,
   double distanceSum)
  {
    Graph& graph = this->getGraph();
    const Models::GroupFormModel<Graph>* model =
      dynamic_cast<const Models::GroupFormModel<Graph>*>(this->getModel());
    bool delPrevious = false;
    if (!previous_nodes) {
      previous_nodes = new std::vector<vertex_descriptor>;
      delPrevious = true;
    }

    if (verbose >=2) {
      std::cout << "random_state_walk: original_node "
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

    out_edge_iterator oi, oi_end;;

    if (distanceSum == 0.) {
      for (tie(oi, oi_end) =
             boost::out_edges(source_node, graph);
           oi != oi_end; ++oi) {
        if (graph[source_node].state->getState() ==
            graph[target(*oi, graph)].state->getState()) {
          bool previous = false;
          for (unsigned int i = 0;
               i < previous_nodes->size() && !previous; ++i) {
            if (target(*oi, graph) == (*previous_nodes)[i]) {
              previous = true;
            }
          }
          if (!previous) {
            if (verbose >=2) {
              std::cout << "checking neighbour " << target(*oi, graph)
                        << std::endl;
            }
            distanceSum +=
              1 - (model->distance(graph[original_node].state,
                                   graph[target(*oi, graph)].state));
          }
        }
      }
    }

    if (distanceSum > 0.) {
      double randSameStateNeighbour = (randGen)() * distanceSum;
      tie(oi, oi_end) = boost::out_edges(source_node, graph);      
      double compareDistance = 0.;
      while (compareDistance == 0. || compareDistance < randSameStateNeighbour) {
        if (graph[source_node].state->getState() ==
            graph[target(*oi, graph)].state->getState()) {
          bool previous = false;
          for (unsigned int i = 0;
               i < previous_nodes->size() && !previous; ++i) {
            if (target(*oi, graph) == (*previous_nodes)[i]) {
              previous = true;
            }
          }
          if (!previous) {
            if (verbose >=2) {
              std::cout << "checking neighbour " << target(*oi, graph)
                        << std::endl;
            }
            compareDistance +=
              1 - model->distance(graph[original_node].state,
                                  graph[target(*oi, graph)].state);
          }
        }
        ++oi;
      }
      --oi;
      if (edge(original_node, target(*oi, graph), graph).second) {
        previous_nodes->push_back(source_node);
        target_node = random_state_walk(original_node, target(*oi, graph),
                                        previous_nodes);
      } else {
        target_node = new vertex_descriptor;
        *target_node = target(*oi, graph);
      }
    }

    if (verbose >=2) {
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

  //----------------------------------------------------------
  /*! \brief Perform an update and process one event.
  
  \return true if an event is processed, false if no events can happen or
  something goes wrong 
  */
  template <typename RandomGenerator, typename Graph>
  bool RewireSimulator<RandomGenerator, Graph>::updateState()
  {
    Graph& graph = this->getGraph();
    const Models::GroupFormModel<Graph>* model =
      dynamic_cast<const Models::GroupFormModel<Graph>*>(this->getModel());

    double randRewire = randGen();
    if (randRewire < rewireProb) {
      // rewiring stage
      if (verbose >= 2) {
        std::cout << "Rewiring" << std::endl;
      }
      
      // find nodes which have a state
      std::vector<vertex_descriptor> state_nodes;
      
      vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        if (graph[*vi].state->getState() > 0) state_nodes.push_back(*vi);
      }

      if (state_nodes.size() > 0) {

        // pick one at random
        unsigned int randStateNode =
          static_cast<unsigned int>((randGen)() * state_nodes.size());
        
        if (verbose >=2) {
          std::cout
            << "Randomly picked node " << state_nodes[randStateNode] << " ("
            << model->printState(graph[state_nodes[randStateNode]].state)
            << ")" << std::endl;
        }
        
        // pick a neighbour of a different state and large distance
        double differentDistanceSum = 0;
        double sameDistanceSum = 0;

        // loop over all member of group and update ties
        out_edge_iterator oi, oi_end;;
        
        // find least similar neighbour and same state neighbours
          
        for (tie(oi, oi_end) =
               boost::out_edges(state_nodes[randStateNode], graph);
             oi != oi_end; ++oi) {
          if (graph[target(*oi, graph)].state->getState() !=
              graph[state_nodes[randStateNode]].state->getState()) {
            differentDistanceSum +=
              (model->distance(graph[state_nodes[randStateNode]].state,
                               graph[target(*oi, graph)].state));
          } else {
            sameDistanceSum += 1 -
              (model->distance(graph[state_nodes[randStateNode]].state,
                               graph[target(*oi, graph)].state));
          }
        }

        // pick a little similar neighbour at random
        if (differentDistanceSum > 0. && sameDistanceSum > 0.) {
          // pick a little similar neighbour at random
          double randNeighbour = ((randGen)() * differentDistanceSum);
          tie(oi, oi_end) = boost::out_edges(state_nodes[randStateNode], graph);
          double compareDistance =
            model->distance(graph[state_nodes[randStateNode]].state,
                            graph[target(*oi, graph)].state);
          while (compareDistance < randNeighbour) {
            compareDistance +=
              model->distance(graph[state_nodes[randStateNode]].state,
                              graph[target(*(++oi), graph)].state);
          }
          vertex_descriptor cut_node = target(*oi, graph);

          // find a new node of the same state by random walk
          vertex_descriptor* target_node =
            random_state_walk
            (state_nodes[randStateNode], state_nodes[randStateNode], 0,
             sameDistanceSum);

          if (target_node) {
            if (verbose >=2) {
              std::cout << "removing link to " << cut_node << " ("
                        << model->printState(graph[cut_node].state) << ")"
                        << std::endl;
            }
            // remove weakest link
            boost::remove_edge
              (edge(state_nodes[randStateNode], cut_node, graph).first, graph);
            // link to group node
            if (verbose >=2) {
              std::cout << "adding link to " << *target_node << " ("
                        << model->printState(graph[*target_node].state) << ")"
                        << std::endl;
            }
            boost::add_edge(state_nodes[randStateNode], *target_node, graph);
            delete target_node;
          } else {
            if (verbose >=2) {
              std::cout << "Random walk for same state neigbhours of "
                        << state_nodes[randStateNode] << " failed."
                        << std::endl;
            }
          }
        } else {
          if (verbose >= 2) {
            if (differentDistanceSum == 0.) { 
              std::cout << "No neighbours of different state" << std::endl;
            }
            if (sameDistanceSum == 0.) {            
              std::cout << "No neighbours of same state" << std::endl;
            }
          }
        }
      } else {
        if (verbose >= 2) {
          std::cout << "No node with state" << std::endl;
        }
      }
    }

    double randUpdate = randGen();
    if (randUpdate < updateProb) {
      // update stage
      if (verbose >= 2) {
        std::cout << "Updating" << std::endl;
      }

      // find nodes which have a state
      std::vector<vertex_descriptor> state_nodes;
      
      vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        if (graph[*vi].state->getState() > 0) state_nodes.push_back(*vi);
      }

      if (state_nodes.size() > 0) {

        // pick node with a state at random
        unsigned int randStateNode =
          static_cast<unsigned int>((randGen)() * state_nodes.size());
        

        if (verbose >=2) {
          std::cout
            << "Randomly picked node " << state_nodes[randStateNode] << " ("
            << model->printState(graph[state_nodes[randStateNode]].state)
            << ")" << std::endl;
        }
        
        // pick a neighbour of a different state and low distance
        double closenessSum = 0;

        // loop over all member of group and update ties
        out_edge_iterator oi, oi_end;;
        
        // find most similar neighbour
          
        for (tie(oi, oi_end) =
               boost::out_edges(state_nodes[randStateNode], graph);
             oi != oi_end; ++oi) {
          if (graph[target(*oi, graph)].state->getState() !=
              graph[state_nodes[randStateNode]].state->getState()) {
            closenessSum += 1 -
              (model->distance(graph[state_nodes[randStateNode]].state,
                               graph[target(*oi, graph)].state));
          }
        }

        if (closenessSum > 0) {
          // pick most similar neighbour at random
          
          double randNeighbour = ((randGen)() * closenessSum);
          tie(oi, oi_end) = boost::out_edges(state_nodes[randStateNode], graph);
          double compareDistance =
            1 - model->distance(graph[state_nodes[randStateNode]].state,
                                graph[target(*oi, graph)].state);
          while (compareDistance < randNeighbour) {
            compareDistance += 1 -
            model->distance(graph[state_nodes[randStateNode]].state,
                            graph[target(*(++oi), graph)].state);
          }

          vertex_descriptor update_node = target(*oi, graph);
          
          double randAccept = randGen();
          if (verbose >=2) {
            std::cout << "Update request to " << update_node << " ("
                      << model->printState(graph[update_node].state)
                      << ")";
          }
          if (randAccept < model->getAcceptance()) {
            // accept invitation
            GroupFormState* myState =
              dynamic_cast<GroupFormState*>(graph[update_node].state);
            myState->setState
              (graph[state_nodes[randStateNode]].state->getState());
            // add to active nodes
            if (verbose >=2) {
              std::cout << " accepted." << std::endl;
            }
          } else {
            if (verbose >=2) {
              std::cout << " refused." << std::endl;
            }
          }
        } else {
          if (verbose >=2) {
            std::cout << "no close neighbours" << std::endl;
          }
        }
      } else {
        if (verbose >= 2) {
          std::cout << "No node with state" << std::endl;
        }
      }
    }

    double randRandomise = randGen();
    if (randRandomise < randomiseProb) {
      // randomisation stage
      if (verbose >= 2) {
        std::cout << "Randomising" << std::endl;
      }

      vertex_descriptor v = random_vertex(graph, randGen);
      if (verbose >=2) {
        std::cout << "Randomly picked vertex " << v << std::endl;
      }

      unsigned int randState =
        static_cast<unsigned int>((randGen)() * (model->getStates())) + 1;
      if (verbose >=2) {
        std::cout << "Randomly picked state " << randState << std::endl;
      }

      GroupFormState* myState = dynamic_cast<GroupFormState*>(graph[v].state);
      myState->setState(randState);
    }

    this->updateTime(1.);

    return true;
  }

  //----------------------------------------------------------
  /*! \brief Print the state of the simulation. 
  
  Prints the state of all vertices, as well as the events that can happen to
  them and at which rate.
  */
  template <typename RandomGenerator, typename Graph>
  void RewireSimulator<RandomGenerator, Graph>::print()
  {
    std::cout << std::endl;
  }
  
}

#endif
