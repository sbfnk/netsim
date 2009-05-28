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
    double randomRewiring;

    unsigned int rewireCounter;
    unsigned int updateCounter;
    unsigned int randomiseCounter;
    unsigned int randomRewireCounter;
    unsigned int counter;

    bool rewireEdges;
    bool updateEdges;
    bool volatility;
    bool traits;
    bool updatingVolatility;
    bool acceptance;
    bool pullUpdating;
    bool randomiseNew;
    bool randomWalk;
    unsigned int recordEffectiveRates;

    unsigned int highestState;
  
  };

  template <typename RandomGenerator, typename Graph>
  RewireSimulator<RandomGenerator, Graph>::
  RewireSimulator(RandomGenerator& r, Graph& g,
		     unsigned int v) :
    Simulator<Graph>(g, v), randGen(r, boost::uniform_real<> (0,1)),
    active(), verbose(v), rewireEdges(false), updateEdges(false),
    volatility(false), traits(false), updatingVolatility(false),
    acceptance(false), pullUpdating(false), randomiseNew(false),
    randomWalk(true), recordEffectiveRates(false)
  {
    this->simulator_options.add_options()
      ("rewire-prob,p",po::value<double>()->default_value(0.),
       "probability of rewiring")
      ("state-update-prob,q",po::value<double>()->default_value(0.),
       "probability of state update")
      ("state-randomise-prob,r",po::value<double>()->default_value(0.),
       "probability of random state assignment")
      ("random-rewiring,w",po::value<double>()->default_value(0.),
       "fraction of rewiring which is random")
      ("edge-based-rewiring",
       "pick random edges rather than nodes in rewiring")
      ("edge-based-updating",
       "pick random edges rather than nodes in updating")
      ("volatility,o",
       "have nodes differ in tendency to rewire")
      ("traits",
       "have nodes differ in traits")
      ("updating-volatility",
       "volatility used also in updating")
      ("acceptance",
       "have nodes differ in tendency to accept updating")
      ("pull-updating",
       "pull rather than push on updating")
      ("randomise-new",
       "randomise to new states (invalidate --states)")
      ("no-random-walk",
       "don't do random walk on rewiring")
      ;
    this->recorder_options.add_options()
      ("component-dist,t",po::value<double>(),
       "write component distribution in CompDist directory at arg timesteps")
      ("components",po::value<double>(),
       "write actual components comp directory at arg timesteps")
      ("modularity",po::value<double>(),
       "write modularity at arg timesteps")
      ("same-state",po::value<double>(),
       "write same-state fraction in components at arg timesteps")
      ("community",po::value<double>(),
       "write community distribution in comm directory at arg timesteps")
      ("effective-rates",po::value<unsigned int>(),
       "write effective rates to rates.sim.dat at arg timesteps")
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

    const Models::GroupFormModel<Graph>* model =
      dynamic_cast<const Models::GroupFormModel<Graph>*>(this->getModel());
    highestState = model->getStates();

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
    randomRewiring = vm["random-rewiring"].as<double>();
    if (randomRewiring == 0.) {
      std::cerr << "WARNING: random rewiring fraction 0" << std::endl;
    }
    if (vm.count("component-dist")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_comp_dist<Graph>(this->getGraph()),
          vm["component-dist"].as<double>()));
    }
    if (vm.count("components")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_components<Graph>(this->getGraph()),
          vm["components"].as<double>()));
    }
    if (vm.count("modularity") || vm.count("community")) {
      double rate = -1;
      if (vm.count("modularity") || vm.count("same-state")) {
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
           vm.count("community"), vm.count("modularity"), verbose),
          rate));
    }
    if (vm.count("same-state")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_same_state_components<Graph>(this->getGraph(), verbose),
          vm["same-state"].as<double>()));
    }
    if (vm.count("effective-rates")) {
      recordEffectiveRates = vm["effective-rates"].as<unsigned int>();
      rewireCounter = 0;
      updateCounter = 0;
      randomiseCounter = 0;
      randomRewireCounter = 0;
      counter = 0;
    } else {
      recordEffectiveRates = 0;
    }
    if (vm.count("volatility")) volatility = true;
    if (vm.count("traits")) traits = true;
    if (vm.count("updating-volatility")) updatingVolatility = true;
    if (vm.count("acceptance")) {
      if (updatingVolatility) {
        std::cerr << "WARNING: acceptance makes no sense with "
                  << "updating-volatility" << std::endl;
      } else {
        acceptance = true;
      }
    }
    if (vm.count("edge-based-rewiring")) rewireEdges = true;
    if (vm.count("edge-based-updating")) updateEdges = true;
    if (vm.count("pull-updating")) {
      if (rewireEdges) {
        std::cerr << "WARNING: pull-updating makes no sense with "
                  << "edge-based-updating" << std::endl;
      } else {
        pullUpdating = true;
      }
    }
    if (vm.count("randomise-new")) {
      if (highestState > 0) {
        std::cerr << "WARNING: randomise-new makes no sense with "
                  << "states > 0" << std::endl;
      } else {
        randomiseNew = true;
      }
    }
    if (vm.count("no-random-walk")) {
      randomWalk = false;
    } else {
      randomWalk = true;
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
    Graph& graph = this->getGraph();
    Models::GroupFormModel<Graph>* model =
      dynamic_cast<Models::GroupFormModel<Graph>*>(this->getModel());
    vertex_iterator vi, vi_end;

    // check if we need to add states to model
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      while (graph[*vi].state->getState() > model->getVertexStates().size()-1) {
        model->addState();
      }
    }
    
    // set trait for each vertex 
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      std::vector<double> randTraits;
      for (unsigned int i = 0; i < model->getTraitDim(); ++i) {
        if (traits) {
          randTraits.push_back((randGen)());
        } else {
          randTraits.push_back(.0);
        }
      }
      GroupFormState* myState = dynamic_cast<GroupFormState*>(graph[*vi].state);
      myState->setTrait(randTraits);
      if (verbose >= 2) {
        std::cout << "Associating trait (";
	for (unsigned int i = 0; i < model->getTraitDim(); ++i) {
	  if (i>0) std::cout << ",";
	  std::cout << myState->getTrait(i);
	}
	std::cout  << ")";
      }
      if (volatility) {
        double randVol = randGen();
        myState->setVolatility(randVol);
        if (verbose >= 2) {
          if (acceptance) {
            std::cout << ",";
          } else {
            std::cout << " and";
          }
          std::cout << " volatility " << randVol;
        }
      }
      if (acceptance) {
        double randAcc = randGen();
        myState->setAcceptance(randAcc);
        if (verbose >= 2) {
          std::cout << " and acceptance " << randAcc;
        }
      }
      if (verbose >=2) {
        std::cout << " with vertex " << *vi << std::endl;
      }
    }

    //initialise simulator
    Simulator<Graph>::initialise();
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
    
    vertex_descriptor* target_node(0);

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
            if (traits) {
              distanceSum +=
                1 - (model->distance(graph[original_node].state,
                                     graph[target(*oi, graph)].state));
            } else {
              ++distanceSum;
            }
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
            if (traits) {
              compareDistance +=
                1 - model->distance(graph[original_node].state,
                                    graph[target(*oi, graph)].state);
            } else {
              ++compareDistance;
            }
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
    Models::GroupFormModel<Graph>* model =
      dynamic_cast<Models::GroupFormModel<Graph>*>(this->getModel());

    double randRewire = randGen();

    if (randRewire < rewireProb) {
      // rewiring stage
      if (verbose >= 2) {
        std::cout << "Rewiring" << std::endl;
      }

      bool rewireRandomly = false;

      vertex_descriptor* source_node(0);
      vertex_descriptor cut_node;
      vertex_descriptor* target_node(0);

      // pick edge to rewire
      if (rewireEdges) {
        edge_descriptor e = random_edge(graph, randGen);

        if (volatility) {
          double volSum =
            model->getVolatility(graph[source(e, graph)].state) +
            model->getVolatility(graph[target(e, graph)].state);
          if (!(randGen() * volSum <
                model->getVolatility(graph[source(e, graph)].state))) {
            // swap edge around
            e = edge(target(e, graph), source(e, graph), graph).first;
          }
        }
        
        source_node = new vertex_descriptor;
	*source_node = source(e, graph);
        cut_node = target(e, graph);

      } else {

        // find nodes which have a state
        std::vector<vertex_descriptor> state_nodes;
        double volatilitySum(0.);

        
        vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
          if (graph[*vi].state->getState() > 0) {
            state_nodes.push_back(*vi);
            if (volatility) {
              volatilitySum +=
                model->getVolatility(graph[*vi].state);
            }
          }
        }
        
        if (state_nodes.size() > 0) {
          
          // pick one at random
          if (volatility) {
            double compareVolatility =
              model->getVolatility(graph[state_nodes[0]].state);
            double randVolatility = ((randGen)() * volatilitySum);
            unsigned int i = 0;
            while (compareVolatility < randVolatility) {
              compareVolatility +=
                model->getVolatility(graph[state_nodes[++i]].state);
            }
            source_node = new vertex_descriptor;
            *source_node = state_nodes[i];
          } else {
            unsigned int randStateNode = static_cast<unsigned int>
              ((randGen)() * state_nodes.size());
            source_node = new vertex_descriptor;
            *source_node = state_nodes[randStateNode];
          }
          if (verbose >=2) {
            std::cout
              << "Randomly picked node " << *source_node << " ("
              << model->printState(graph[*source_node].state)
              << ")";
            if (volatility) {
              std::cout << " [volatility "
                        << model->getVolatility(graph[*source_node].state)
                        << "]" << std::endl;
            }
            std::cout << std::endl;
          }
          
          // pick a neighbour of a different state and large distance
          double differentDistanceSum = 0;
          // loop over all member of group and update ties
          out_edge_iterator oi, oi_end;;
          
          // find least similar neighbour
          
          for (tie(oi, oi_end) =
                 boost::out_edges(*source_node, graph);
               oi != oi_end; ++oi) {
            if (graph[target(*oi, graph)].state->getState() !=
                graph[*source_node].state->getState()) {
              if (traits) {
                differentDistanceSum +=
                  (model->distance(graph[*source_node].state,
                                   graph[target(*oi, graph)].state));
              } else {
                ++differentDistanceSum;
              }
            }
          }
          
          // pick a little similar neighbour at random
          if (differentDistanceSum > 0.) {
            // pick a little similar neighbour at random
            double randNeighbour = ((randGen)() * differentDistanceSum);
            tie(oi, oi_end) = boost::out_edges(*source_node, graph);
            double compareDistance;
            if (traits) {
              compareDistance =
                model->distance(graph[*source_node].state,
                                graph[target(*oi, graph)].state);
            } else {
              compareDistance = 1;
            }
            while (compareDistance < randNeighbour) {
              if (traits) {
                compareDistance +=
                  model->distance(graph[*source_node].state,
                                  graph[target(*(++oi), graph)].state);
              } else {
                ++compareDistance;
              }
            }
            cut_node = target(*oi, graph);
          } else {
            if (verbose >= 2) {
              std::cout << "No neighbours of different state" << std::endl;
            }
            delete source_node;
	    source_node = 0;
          }
        } else {
          if (verbose >= 2) {
            std::cout << "No node with state" << std::endl;
          }
        }
      }

      if (source_node) {
        if (randomRewiring > 0 && randGen() < randomRewiring) {
          // random rewiring
          target_node = new vertex_descriptor;
          do {
            *target_node = random_vertex(graph, randGen);
          } while (*target_node == *source_node ||
                   edge(*source_node, *target_node, graph).second);
	  rewireRandomly = true;
        } else {
          // local rewiring
          // only do something if the two nodes don't have the same state
	  
          if (graph[*source_node].state->getState() > 0 &&
              (graph[*source_node].state->getState() !=
               graph[cut_node].state->getState())) {
  
	    if (randomWalk) {
	      double sameDistanceSum(0.);
	      out_edge_iterator oi, oi_end;;
	      for (tie(oi, oi_end) =
		     boost::out_edges(*source_node, graph);
		   oi != oi_end; ++oi) {
		if (graph[target(*oi, graph)].state->getState() ==
		    graph[*source_node].state->getState()) {
		  if (traits) {
		    sameDistanceSum += 1 -
		      (model->distance(graph[*source_node].state,
				       graph[target(*oi, graph)].state));
		  } else {
		    ++sameDistanceSum;
		  }
		}
	      }
	      
	      if (sameDistanceSum > 0.) {
		// find a new node of the same state by random walk
		target_node =
		  random_state_walk
		  (*source_node, *source_node, 0,
		   sameDistanceSum);
		if (!target_node && verbose >=2) {
		  std::cout << "Random walk for same state neighbours of "
			    << *source_node << " failed."
			    << std::endl;
		}            
	      } else {
		if (verbose >= 2) {
                  std::cout << "No neighbours of same state" << std::endl;
		}
	      }
	    } else {
	      // no random walk
	      // find vertices of the same state
	      std::vector<vertex_descriptor> v;
	      vertex_iterator vi, vi_end;
	      
	      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
		if ((*vi != *source_node) &&
		    (graph[*vi].state->getState() == 
		     graph[*source_node].state->getState()) &&
		    (!edge(*source_node, *vi, graph).second)) {
		  v.push_back(*vi);
		}
	      }
	      if (v.size() > 0) {
		//pick a random number
		unsigned int randInt = 
		  static_cast<unsigned int>(randGen() * v.size());
		target_node = new vertex_descriptor;
		*target_node = v[randInt];
	      }
	    }
	  }
        }
        if (target_node) {
          if (verbose >=2) {
            std::cout << "removing link to " << cut_node << " ("
                      << model->printState(graph[cut_node].state) << ")"
                      << std::endl;
          }
          // remove weakest link
          boost::remove_edge
            (edge(*source_node, cut_node, graph).first, graph);
          // link to group node
          if (verbose >=2) {
            std::cout << "adding link to " << *target_node << " ("
                      << model->printState(graph[*target_node].state) << ")"
                      << std::endl;
          }
          boost::add_edge(*source_node, *target_node, graph);
          delete target_node;
	  if (rewireRandomly) {
	    ++randomRewireCounter;
	  } else {
	    ++rewireCounter;
	  }
	}
	delete source_node;
      }   
    } 

    double randUpdate = randGen();
    if (randUpdate < updateProb) {
      // update stage
      if (verbose >= 2) {
        std::cout << "Updating" << std::endl;
      }

      vertex_descriptor source_node;
      vertex_descriptor* target_node(0);

      if (updateEdges) {
        // edge-based updating
        edge_descriptor e = random_edge(graph, randGen);

        if (updatingVolatility) {
          double volSum =
            model->getVolatility(graph[source(e, graph)].state) +
            model->getVolatility(graph[target(e, graph)].state);
          if (randGen() * volSum <
              model->getVolatility(graph[source(e, graph)].state)) {
            // swap edge around
            e = edge(target(e, graph), source(e, graph), graph).first;
          }
        } else if (acceptance) {
          double accSum =
            model->getAcceptance(graph[source(e, graph)].state) +
            model->getAcceptance(graph[target(e, graph)].state);
          if (randGen() * accSum <
                model->getAcceptance(graph[source(e, graph)].state)) {
            // swap edge around
            e = edge(target(e, graph), source(e, graph), graph).first;
          }
        }

        if (graph[source(e, graph)].state->getState() > 0 &&
            (graph[source(e, graph)].state->getState() !=
             graph[target(e, graph)].state->getState())) {
              
          source_node = source(e, graph);
          target_node = new vertex_descriptor;
          *target_node = target(e, graph);
        }
      } else {
        // state-based updating
      
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

          source_node = state_nodes[randStateNode];
          
          if (verbose >=2) {
            std::cout
              << "Randomly picked node " << source_node << " ("
              << model->printState(graph[source_node].state)
              << ")" << std::endl;
          }
        
          // pick a neighbour of a different state and low distance
          double closenessSum = 0;
          
          // loop over all member of group and update ties
          out_edge_iterator oi, oi_end;;
          
          // find most similar neighbour
          
          for (tie(oi, oi_end) =
                 boost::out_edges(source_node, graph);
               oi != oi_end; ++oi) {
            if (graph[target(*oi, graph)].state->getState() !=
                graph[source_node].state->getState()) {
              if (traits) {
                closenessSum += 1 -
                  (model->distance(graph[source_node].state,
                                   graph[target(*oi, graph)].state));
              } else {
                ++closenessSum;
              }
            }
          }

          if (closenessSum > 0) {
            // pick most similar neighbour at random
            
            double randNeighbour = ((randGen)() * closenessSum); 
            tie(oi, oi_end) = boost::out_edges(source_node, graph);
            double compareDistance;
            if (traits) {
              compareDistance =
                1 - model->distance(graph[source_node].state,
                                    graph[target(*oi, graph)].state);
            } else {
              compareDistance = 1;
            }
            while (compareDistance < randNeighbour) {
              if (graph[target(*oi, graph)].state->getState() !=
                  graph[source_node].state->getState()) {
                if (traits) {
                  compareDistance += 1 -
                    model->distance(graph[state_nodes[randStateNode]].state,
                                    graph[target(*(++oi), graph)].state);
                } else {
                  ++compareDistance;
                }
              }
              ++oi;
            }
            
            target_node = new vertex_descriptor;
            if (pullUpdating) {
              *target_node = source_node;
              source_node = target(*oi, graph);
            } else {
              *target_node  = target(*oi, graph);
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
      if (target_node) {
        double randAccept = randGen();
        if (verbose >=2) {
          std::cout << "Update request to " << *target_node << " ("
                    << model->printState(graph[*target_node].state)
                    << ")";
        }
        if (randAccept < model->getAcceptance()) {
          // accept invitation
          GroupFormState* myState =
            dynamic_cast<GroupFormState*>(graph[*target_node].state);
          myState->setState
            (graph[source_node].state->getState());
	  ++updateCounter;
          if (verbose >=2) {
            std::cout << " accepted." << std::endl;
          }
        } else {
          if (verbose >=2) {
            std::cout << " refused." << std::endl;
          }
        }
        delete target_node;
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

      unsigned int newState;
      if (randomiseNew) {
        newState = ++highestState;
        model->addState();
      } else {
        newState =
          static_cast<unsigned int>((randGen)() * (model->getStates())) + 1;
      }

      ++randomiseCounter;
      if (verbose >=2) {
        std::cout << "Assigning state " << model->getVertexState(newState) 
<< std::endl;
      }
      GroupFormState* myState = dynamic_cast<GroupFormState*>(graph[v].state);
      myState->setState(newState);
    }

    ++counter;
    if (recordEffectiveRates > 0 && counter % recordEffectiveRates == 0) {
      std::ofstream outputFile;
      std::string fileName = (this->getDir() + "/rates.sim.dat");
      
      try {
        outputFile.open(fileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::ate);
      }
      catch (std::exception &e) {
	std::cerr << "Unable to open output file: " << e.what() << std::endl;
	std::cerr << "Will not write effective rates to file." << std::endl;
      }
      
      outputFile << this->getTime() << '\t'
		 << (rewireCounter / static_cast<double>(counter)) << '\t'
		 << (updateCounter / static_cast<double>(counter)) << '\t'
		 << (randomiseCounter / static_cast<double>(counter)) << '\t'
		 << (randomRewireCounter / static_cast<double>(counter));

      outputFile << std::endl;
      
      outputFile.close();

      counter = 0;
      rewireCounter = 0;
      updateCounter = 0;
      randomiseCounter = 0;
      randomRewireCounter = 0;
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
