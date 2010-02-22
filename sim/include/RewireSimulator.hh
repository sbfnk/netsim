/*! \file RewireSimulator.hh
  \brief The Simulators::RewireSimulator class.
*/

#ifndef REWIRESIMULATOR_HH
#define REWIRESIMULATOR_HH

#include <math.h>
#include <errno.h>
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
    typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
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

    enum {REWIRING, STATESPREAD, INNOVATION, RANDOMREWIRING};
  
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
      
    bool initialise();
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
      if (ret == true) {
        ret &= ((burnWait == 0) || (tempStates.find(0) == tempStates.end()) ||
                (static_cast<int>(tempStates.at(0).size()) < burnWait) ||
                ((burnWait < 0) && (this->getTime() < burnTime)));
        if (ret == true) {
          if (stopComponent > 0) {
            unsigned int largest_component = stopComponent;
            // find component dist writer if available
            bool found = false;
            for (unsigned int i = 0;
                 ((!found) && i < (this->statRecorders.size())); ++i) {
              if (this->statRecorders[i]->getName() == "component-dist") {
                const write_comp_dist<Graph>* w =
                  dynamic_cast<const write_comp_dist<Graph>*>
                  (this->statRecorders[i]->getStatFunc());
                if (w) {
                  largest_component = w->getLargest();
                  found = true;
                }
              }
            }
            if (!found) {
              largest_component = write_component_dist(this->getGraph());
            }
            ret &= (largest_component < stopComponent); 
          }
        }
      }
      return ret;
    }

  private:

    uniform_gen randGen; //!< The random generator to be used for choosing events

    std::vector<vertex_descriptor> active; //!< Active nodes

    unsigned int verbose;
    unsigned int stopComponent;
    int burnWait;
    double burnTime;

    double rewireProb;
    double updateProb;
    double randomiseProb;
    double randomRewiring;

    double switchOff;

    unsigned int rewireCounter;
    unsigned int rewirefcCounter;
    unsigned int rewirefnCounter;
    unsigned int rewiressCounter;
    unsigned int rewiredsCounter;
    unsigned int updateCounter;
    unsigned int randomiseCounter;
    unsigned int randomRewireCounter;
    unsigned int counter;

    unsigned int numEdges;
    unsigned int numVertices;

    bool recordInitiator;
    bool saveInitiator;
    bool groupLifeTimes;

    double recordEffectiveRates;
    double recordEffectiveTimer;

    unsigned int highestState;

    std::vector<double> rates;

    std::vector<edge_descriptor> tempEdges; // for faster access
    std::vector<unsigned int> tempVertices; // for faster access
    std::map<unsigned int, std::set<unsigned int> > tempStates; // faster

    // the times at which groups are initiated 
    std::vector<double> groupInitiationTimes;
    // the groups which nodes belong to when initiating new groups
    std::vector<unsigned int> groupInitiationPrevious; 

    double rateSum;
  
  };

  template <typename RandomGenerator, typename Graph>
  RewireSimulator<RandomGenerator, Graph>::
  RewireSimulator(RandomGenerator& r, Graph& g,
		     unsigned int v) :
    Simulator<Graph>(g, v), randGen(r, boost::uniform_real<> (0,1)),
    active(), verbose(v), burnWait(0), burnTime(0.),
    recordInitiator(false), saveInitiator(false),
    groupLifeTimes(false), recordEffectiveRates(0.), recordEffectiveTimer(0.), 
    rateSum(0.)
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
      ;
    this->recorder_options.add_options()
      ("component-dist,t",po::value<double>(),
       "write component distribution in CompDist directory at arg timesteps")
      ("components",po::value<double>(),
       "write actual components comp directory at arg timesteps")
      ("modularity",po::value<double>(),
       "write modularity at arg timesteps")
      ("group-modularity",po::value<double>(),
       "write group modularity at arg timesteps")
      ("same-state",po::value<double>(),
       "write same-state fraction in components at arg timesteps")
      ("community",po::value<double>(),
       "write community distribution in comm directory at arg timesteps")
      ("effective-rates",po::value<double>(),
       "write effective rates to rates.sim.dat at arg timesteps")
      ("record-initiators",
       "record initiators in Initiator directory")
      ("group-lifetimes",
       "save group lifetimes to lifetimes.sim.dat")
      ("save-initiations",
       "save data at group initiations")
      ("switchoff", po::value<double>()->default_value(0.),
       "time at which to switch of network processes (0 for never)")
      ;
    this->stop_options.add_options()
      ("cmin", po::value<unsigned int>()->default_value(0),
       "limit to the size of the largest component (stop if drops to that value")
      ("burn-wait", po::value<int>()->default_value(0),
       "wait for burn-in to finish (i.e. wait until noone is in null state")
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

    highestState = 0;

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

    switchOff = vm["switchoff"].as<double>();
    
    if (vm.count("component-dist")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_comp_dist<Graph>(this->getGraph()),
          vm["component-dist"].as<double>(), "component-dist"));
    }
    if (vm.count("components")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_comp<Graph>(this->getGraph()),
          vm["components"].as<double>(), "components"));
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
          rate, "community"));
    }
    if (vm.count("group-modularity")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_group_modularity<Graph, Model<Graph> >
	  (this->getGraph(), *(this->getModel()), verbose),
          vm["group-modularity"].as<double>(), "group-modularity"));
    }
    if (vm.count("same-state")) {
      this->statRecorders.push_back
        (new StatRecorder<Graph>
         (new write_same_state_components<Graph>(this->getGraph(), verbose),
          vm["same-state"].as<double>(), "same-state"));
    }
    if (vm.count("effective-rates")) {
      recordEffectiveRates = vm["effective-rates"].as<double>();
      recordEffectiveTimer = 0;
      rewireCounter = 0;
      rewirefcCounter = 0;
      rewirefnCounter = 0;
      rewiressCounter = 0;
      rewiredsCounter = 0;
      updateCounter = 0;
      randomiseCounter = 0;
      randomRewireCounter = 0;
    } else {
      recordEffectiveRates = 0;
    }
    if (vm.count("record-initiators")) {
      recordInitiator = true;
    }
    if (vm.count("save-initiations")) {
      saveInitiator = true;
    }
    if (vm.count("group-lifetimes")) {
      groupLifeTimes = true;
    }
    if (vm.count("burn-wait")) {
      burnWait = vm["burn-wait"].as<int>();
    }
    return ret;
  }

  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
    
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  bool RewireSimulator<RandomGenerator, Graph>::initialise()
  {
    Graph& graph = this->getGraph();
    Models::GroupFormModel<Graph>* model =
      dynamic_cast<Models::GroupFormModel<Graph>*>(this->getModel());
    vertex_iterator vi, vi_end;
    edge_iterator ei, ei_end;

    // This will not change
    numEdges = num_edges(graph);
    numVertices = num_vertices(graph);
    //tempStates.insert(std::make_pair(0, std::set<unsigned int>()));
    groupInitiationTimes.push_back(0.);
    groupInitiationPrevious.push_back(0);

    for (tie(ei, ei_end) = edges(graph); ei != ei_end; ei++) {
      tempEdges.push_back(*ei);
    }

//    unsigned int max_state = num_vertices(graph)/10;
//    for (unsigned int i = 0; i < max_state; ++i) {
//      groupInitiationTimes.push_back(0.);
//    }
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; vi++) {
//      unsigned int random_state = static_cast<unsigned int>(randGen() * (max_state - 1)) + 1;
//      graph[*vi].state->setState(random_state);
      tempVertices.push_back(graph[*vi].state->getState());
      tempStates[graph[*vi].state->getState()].insert(*vi);
    }

    // check if we need to add states to model
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
      while (tempVertices[*vi] > model->getVertexStates().size()-1) {
        model->addState();
      }
    }
    
    rateSum = 
      numEdges * rewireProb +
      numEdges * updateProb +
      numVertices * randomiseProb +
      numEdges * randomRewiring;

    rates.push_back(numEdges * rewireProb / rateSum);
    rates.push_back(numEdges * updateProb / rateSum);
    rates.push_back(numVertices * randomiseProb / rateSum);
    rates.push_back(numEdges * randomRewiring / rateSum);

    //initialise simulator
    return Simulator<Graph>::initialise();
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

    if (rateSum < 1e-8) {
      if (verbose >=2) {
        std::cout << "Nothing can happen. Stopping run. " << std::endl;
      }
      this->updateStats(true);
      return false;
    }

    // draw a random number from (0,1) for the timestep advance
    double timeStep = 0;
    do { 
      errno = 0;
      timeStep = -log(randGen())/rateSum;
    } while (errno > 0);
    this->updateTime(timeStep);

    if (switchOff > 0 && switchOff < this->getTime()) {
      rewireProb = 0;
      randomRewiring = 0;
    }
    
    recordEffectiveTimer += timeStep;
         
    double randEvent = randGen();
    int event = 0;

    for (; (randEvent - rates[event]) > 0; ++event) {
      randEvent -= rates[event];
    }
    
    switch(event) {

    case REWIRING: 

      {

	// rewiring stage
	if (verbose >= 2) {
	  std::cout << "Rewiring" << std::endl;
	}

	// choose a random edge
        unsigned int randEdge = static_cast<unsigned int>
          (randGen() * numEdges * 2);
        
	vertex_descriptor source_node;
	vertex_descriptor target_node;

        if (randEdge < numEdges) {
	  source_node = source(tempEdges[randEdge], graph);
	  target_node = target(tempEdges[randEdge], graph);
        } else {
          randEdge -= numEdges;
	  source_node = target(tempEdges[randEdge], graph);
	  target_node = source(tempEdges[randEdge], graph);
        }

	if (verbose >= 2) {
	  std::cout << "Chose random edge: " << source_node << " (" 
                    << model->printState(graph[source_node].state) 
                    << ") -- " << target_node << " (" 
                    << model->printState(graph[target_node].state) 
                    << ") " << std::endl;
	}

	if (tempVertices[source_node] > 0) {
	  std::set<unsigned int> nonNbSameStates = 
	    tempStates[tempVertices[source_node]];
	  nonNbSameStates.erase(source_node);
	  out_edge_iterator oi, oi_end;
	  for (tie(oi, oi_end) = out_edges(source_node, graph);
	       oi != oi_end; oi++) {
	    if (tempVertices[target(*oi, graph)] == tempVertices[source_node]) {
              if (verbose >= 2) {
                std::cout << "Connected same state node: " << target(*oi, graph) << std::endl;
              }
	      nonNbSameStates.erase(target(*oi, graph));
	    }
	  }
          if (verbose >= 2) {
            std::cout << source_node << " has " << nonNbSameStates.size()
                      << " non nb same state nodes" << std::endl;
          }
	  if (nonNbSameStates.size() > 0) {
	    // reuse randEvent number to select target node

	    unsigned int randomSameState =
	      static_cast<unsigned int>
	      (randEvent / rates[REWIRING] * 
	       nonNbSameStates.size());
	    std::set<unsigned int>::iterator it = 
	      nonNbSameStates.begin();
	    unsigned int i = 0;
	    while (i++ < randomSameState) {
	      it++;
	    }
	    vertex_descriptor newtarget = *it;
	       
	    boost::remove_edge(tempEdges[randEdge], graph);
	    tempEdges[randEdge] =
              boost::add_edge(source_node, newtarget, graph).first;

	    if (verbose >=2) {
	      std::cout << "EVENT: Rewiring edge " << source_node << " (" 
			<< model->printState(graph[source_node].state) 
			<< ") -- " << target_node << " (" 
			<< model->printState(graph[target_node].state) 
			<< ") "
			<< " to new edge " << source_node << " ("
			<< model->printState(graph[source_node].state) 
			<< ") -- " << newtarget << " (" 
			<< model->printState(graph[newtarget].state) << ")"
			<< std::endl;
	    }

            if (graph[source_node].state->getState() ==
                graph[target_node].state->getState()) {
              ++rewiressCounter;
            } else {
              ++rewiredsCounter;
            }
	    ++rewireCounter;
	  } else {
            // no node of the same state unconnected
            ++rewirefcCounter;
          }
	} else {
          // no other node of the same state
          ++rewirefnCounter;
        }
      }

      break;

    case STATESPREAD: 

      {

	// state spread stage
	if (verbose >= 2) {
	  std::cout << "State spread" << std::endl;
	}

	//choose a random edge

        unsigned int randEdge = static_cast<unsigned int>
          (randGen() * numEdges);
        
	vertex_descriptor source_node;
	vertex_descriptor target_node;

        if (randEdge < numEdges) {
	  source_node = source(tempEdges[randEdge], graph);
	  target_node = target(tempEdges[randEdge], graph);
        } else {
          randEdge -= numEdges;
	  source_node = target(tempEdges[randEdge], graph);
	  target_node = source(tempEdges[randEdge], graph);
        }

	if (tempVertices[source_node] > 0) {
          tempStates[tempVertices[target_node]].erase(target_node);
          tempStates[tempVertices[source_node]].insert(target_node);
	  if (groupLifeTimes && 
	      tempStates[tempVertices[target_node]].size() == 0) {
	    if (verbose >= 2) {
	      std::cout << "time: " << this->getTime() << ", Group " 
			<< tempVertices[target_node] << " dies out, lifetime "
			<< (this->getTime() - 
			    groupInitiationTimes[tempVertices[target_node]])
			<< std::endl;
	    }
            std::stringstream s;
            s << groupInitiationPrevious[tempVertices[target_node]]
              << "\t"
              << (this->getTime() -
                  groupInitiationTimes[tempVertices[target_node]]);
	    write_data((this->getDir() + "/lifetimes.sim.dat"),
		       tempVertices[target_node], s.str());
	    tempStates.erase(tempVertices[target_node]);
	  }
	  
	  if (verbose >=2) {
	    std::cout << "EVENT: Spreading state " 
		      << model->printState(graph[source_node].state)
		      << " from node " << source_node << " (" 
                      << model->printState(graph[source_node].state) 
                      << ") to node " << target_node << " (" 
                      << model->printState(graph[target_node].state) 
                      << ")" << std::endl;
	  }
	  GroupFormState* myState =
	    dynamic_cast<GroupFormState*>(graph[target_node].state);
	  myState->setState(tempVertices[source_node]);
          tempVertices[target_node] = tempVertices[source_node];
	  ++updateCounter;
	}
      }

      break;

    case INNOVATION: 

      {

	// innovation stage
	if (verbose >= 2) {
	  std::cout << "Innovation" << std::endl;
	}

	vertex_descriptor v = random_vertex(graph, randGen);
	if (verbose >=2) {
	  std::cout << "Randomly picked vertex " << v << std::endl;
	}

	unsigned int newState;
	newState = ++highestState;
	model->addState();
	groupInitiationTimes.push_back(this->getTime());
	groupInitiationPrevious.push_back(tempVertices[v]);

	tempStates[tempVertices[v]].erase(v);
	if (groupLifeTimes && tempStates[tempVertices[v]].size() == 0) {
	  if (verbose >= 2) {
	    std::cout << "time: " << this->getTime() << ", Group " 
		      << tempVertices[v] << " dies out, lifetime "
		      << (this->getTime() - 
			  groupInitiationTimes[tempVertices[v]])
		      << std::endl;
	  }
          std::stringstream s;
          s << groupInitiationPrevious[tempVertices[v]]
            << "\t"
            << (this->getTime() -
                groupInitiationTimes[tempVertices[v]]);
	  write_data((this->getDir() + "/lifetimes.sim.dat"),
		     tempVertices[v], s.str());
	  tempStates.erase(tempVertices[v]);
	}
	  
	if (verbose >=2) {
	  std::cout << "EVENT: Assigning randomly picked vertex " << v << " ("
                    << model->printState(graph[v].state) << ")"
		    << " new state " 
		    << model->printState(new GroupFormState(newState)) << std::endl;
	  if (groupLifeTimes) {
	    std::cout << " at time " << this->getTime();
	  }
	  std::cout << std::endl;
	}

	GroupFormState* myState =
          dynamic_cast<GroupFormState*>(graph[v].state);
	myState->setState(newState);

        tempVertices[v] = newState;
	tempStates.insert(std::make_pair(newState,std::set<unsigned int>()));
	tempStates[newState].insert(v);
	++randomiseCounter;

        if (recordInitiator) {
          if (!fs::exists(this->getDir()+"/Initiators")) {
            try {
              mkdir((this->getDir()+"/Initiators").c_str(), 0755);
            } 
            catch (std::exception &e) {
              std::cerr << "... unable to create directory "
                        << this->getDir() << "/Initiators" << std::endl;
              std::cerr << "unsetting record-initiators" << std::endl;
              recordInitiator = false;
            }
          }
          std::string fileName = 
            generateFileName(this->getDir() + "/Initiators/group", 
                             highestState, "graph");
          write_graph(this->getGraph(), fileName, *(this->getModel()), 
                      this->getTime());
        }

        if (saveInitiator) {
	  // find data writer if available
	  bool found = false;
	  for (unsigned int i = 0;
	       ((!found) && i < (this->statRecorders.size())); ++i) {
// 	    const write_comp_dist<Graph>* w =
// 	      dynamic_cast<const write_comp_dist<Graph>*>
// 	      (this->statRecorders[i]->getStatFunc());
// 	    if (w) {
            if (this->statRecorders[i]->getName() == "data") {
              this->statRecorders[i]->update(this->getGraph(),
                                             this->getTime(), true);
	    }
	  }
        }

      }
      
      break;

    case RANDOMREWIRING: 
      
      {

	// random rewiring stage
	if (verbose >= 2) {
	  std::cout << "Random rewiring" << std::endl;
	}
	
	// choose a random edge

        unsigned int randEdge = static_cast<unsigned int>
          (randGen() * numEdges);
        
	vertex_descriptor source_node;
	vertex_descriptor target_node;

        if (randEdge < numEdges) {
	  source_node = source(tempEdges[randEdge], graph);
	  target_node = target(tempEdges[randEdge], graph);
        } else {
          randEdge -= numEdges;
	  source_node = target(tempEdges[randEdge], graph);
	  target_node = source(tempEdges[randEdge], graph);
        }

	boost::remove_edge(tempEdges[randEdge], graph);
	
	// find unconnected vertices
	std::vector<vertex_descriptor> unconnected;
	vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
	  if ((*vi != source_node) && 
	      (edge(source_node, *vi, graph).second == false)) {
	    unconnected.push_back(*vi);
	  }
	}

	// choose a random unconnected edge, reusing previously used
	// random nmber

	vertex_descriptor newtarget = 
	  unconnected[static_cast<unsigned int>
		      (randEvent / rates[RANDOMREWIRING] *
                       unconnected.size())];
      
	tempEdges[randEdge] =
          boost::add_edge(source_node, newtarget, graph).first;

	if (verbose >=2) {
	  std::cout << "EVENT: Randomly rewiring edge " << source_node << " (" 
		    << model->printState(graph[source_node].state) 
		    << ") -- " << target_node << " (" 
		    << model->printState(graph[target_node].state) 
		    << ")" << " to new edge " << source_node << " ("
		    << model->printState(graph[source_node].state) 
		    << ") -- " << newtarget << " (" 
		    << model->printState(graph[newtarget].state) << ")"
		    << std::endl;
	}

	++randomRewireCounter;

      }

      break;

    }

    ++counter;

    if (recordEffectiveRates > 0  && 
	recordEffectiveTimer >= recordEffectiveRates) {

      std::vector<double> effectiveRates;
      
      effectiveRates.push_back(rewireCounter / (recordEffectiveTimer *
						numEdges));
      effectiveRates.push_back(updateCounter / (recordEffectiveTimer *
						numEdges));
      effectiveRates.push_back(randomiseCounter / (recordEffectiveTimer *
						   numVertices));
      effectiveRates.push_back(randomRewireCounter / (recordEffectiveTimer *
						      numEdges));
      effectiveRates.push_back(rewirefnCounter / (recordEffectiveTimer *
                                                  numEdges));
      effectiveRates.push_back(rewirefcCounter / (recordEffectiveTimer *
                                                  numEdges));
      effectiveRates.push_back(rewiressCounter / (recordEffectiveTimer *
                                                  numEdges));
      effectiveRates.push_back(rewiredsCounter / (recordEffectiveTimer *
                                                  numEdges));

      write_data(this->getDir() + "/rates.sim.dat", this->getTime(), 
		 effectiveRates);
      
      recordEffectiveTimer = 0.;
      rewireCounter = 0;
      updateCounter = 0;
      randomiseCounter = 0;
      randomRewireCounter = 0;
      rewirefcCounter = 0;
      rewirefnCounter = 0;
      rewiressCounter = 0;
      rewiredsCounter = 0;
    }
  
    if (burnWait == -1) {
      if (tempStates.find(0) == tempStates.end()) {
        burnWait = -2;
        burnTime = 2 * this->getTime();
      } else {
        burnTime = this->getTime();
      }
    }
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
