/*! \file GroupFormSimulator.hh
  \brief The Simulators::GroupFormSimulator class.
*/

#ifndef GROUPFORMSIMULATOR_HH
#define GROUPFORMSIMULATOR_HH

#include <math.h>

#include "Simulator.hh"
#include "Vertex.hh"
#include "GroupFormModel.hh"

#include <boost/graph/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

//! \addtogroup simulator Simulator

namespace Simulators {
  
  /*! \brief The Group formation simulation
    
  \ingroup simulator
  */
  template <typename RandomGenerator, typename Graph>
  class GroupFormSimulator :
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
    GroupFormSimulator(RandomGenerator& r, Graph& g,
                       unsigned int v = 0);
    ~GroupFormSimulator() {;}
      
    void initialise();
    bool updateState();
  
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
  
  };

  template <typename RandomGenerator, typename Graph>
  GroupFormSimulator<RandomGenerator, Graph>::
  GroupFormSimulator(RandomGenerator& r, Graph& g,
		     unsigned int v) :
    Simulator<Graph>(g, v), randGen(r, boost::uniform_real<> (0,1)),
    active(), verbose(v)
  {
    this->recorder_options.add_options()
      ("component-dist,t",po::value<double>(),
       "write component distribution in comp directory at arg timesteps")
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
  bool GroupFormSimulator<RandomGenerator, Graph>::
  parse_options(const po::variables_map& vm)
  { 
    bool ret = Simulator<Graph>::parse_options(vm);
    stopComponent = vm["cmin"].as<unsigned int>();
    if (vm.count("comp")) {
      this->statRecorders.push_back(new StatRecorder<Graph>
                                    (new write_comp_dist<Graph>(this->getGraph()),
                                     vm["comp"].as<double>()));
    }
    return ret;
  }


  //----------------------------------------------------------
  /*! \brief Initialise the simulation.
    
  Initialises the graph with the rates for the possible processes and generates 
  the Tree using these rates.
  */
  template <typename RandomGenerator, typename Graph>
  void GroupFormSimulator<RandomGenerator, Graph>::initialise()
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
      if (verbose >=1) {
        std::cout << "Associating trait (";
	for (unsigned int i = 0; i < model->getTraitDim(); ++i) {
	  if (i>0) std::cout << ",";
	  std::cout << myState->getTrait(i);
	}
	std::cout  << ") with vertex " << *vi << std::endl;
      }
    }

    // initially seed groups
    for (unsigned int i = 1; i <= model->getGroups(); ++i) {
      vertex_descriptor v = random_vertex(graph, randGen);
      if (graph[v].state->getState() == 0) {
        graph[v].state->setState(i);
        active.push_back(v);
        if (verbose >=1) {
          std::cout << "Seeding group " << i << " with vertex " << v << std::endl;
        }
      } else {
        --i;
      }
    }
    
  }

  //----------------------------------------------------------
  /*! \brief Perform an update and process one event.
  
  \return true if an event is processed, false if no events can happen or
  something goes wrong 
  */
  template <typename RandomGenerator, typename Graph>
  bool GroupFormSimulator<RandomGenerator, Graph>::updateState()
  {
    Graph& graph = this->getGraph();
    const Models::GroupFormModel<Graph>* model =
      dynamic_cast<const Models::GroupFormModel<Graph>*>(this->getModel());
    
    if (active.size() > 0) {

      if (verbose >=2) {
        std::cout << "Active nodes: " << std::endl;
        for (unsigned int i = 0; i < active.size(); ++i) {
          std::cout << active[i] << std::endl;
        }
      }

      // group formation stage
      if (verbose >=1) {
        std::cout << "Group formation" << std::endl;
      }

      // choose an active node at random
      unsigned int randActive =
        static_cast<unsigned int>((randGen)() * active.size());
      if (verbose >=2) {
        std::cout << "Active node " << active[randActive] << " ("
                  << model->printState(graph[active[randActive]].state)
                  << ") is sending out invitations" << std::endl;
      }
    
      // invite all neighbours to group
      out_edge_iterator oi, oi_end;;
   
      for (tie(oi, oi_end) = boost::out_edges(active[randActive], graph);
           oi != oi_end; ++oi) {
        if (graph[target(*oi, graph)].state->getState() == 0) {
          double randAccept = randGen();
          if (verbose >=2) {
            std::cout << "Invitation to " << target(*oi, graph) << " ("
                      << model->printState(graph[target(*oi, graph)].state)
                      << ")" << std::endl;
          }
          if (randAccept < model->accept(graph[active[randActive]].state,
                                         graph[target(*oi, graph)].state)) {
            // accept invitation
            GroupFormState* myState =
              dynamic_cast<GroupFormState*>(graph[target(*oi, graph)].state);

            myState->setState(graph[active[randActive]].state->getState());
            // add to active nodes
            active.push_back(target(*oi, graph));
            if (verbose >=2) {
              std::cout << " accepted." << std::endl;
            }
          } else {
            if (verbose >=2) {
              std::cout << " refused." << std::endl;
            }
          }
        }
      }

      // remove randActive from active nodes
      active.erase(active.begin() + randActive);

      // update time if group formation phase is finished
      if (active.size() == 0) this->updateTime(1.);
      
    } else {
      // network update stage
      if (verbose >=1) {
        std::cout << "Network update stage, time=" << this->getTime() 
		  << std::endl;
      }

      // go through members group by group
      std::vector<vertex_descriptor> group;
      
      for (unsigned int i = 1; i <= model->getGroups(); ++i) {
        if (verbose >=2) {
          std::cout << "Group " << i << std::endl;
        }
        group.clear();
        
        vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
          if (graph[*vi].state->getState() == i) group.push_back(*vi);
        }

        // loop over all member of group and update ties
        for (unsigned int i = 0; i < group.size(); ++i) {
          if (verbose >=2) {
            std::cout << "Group member " << i << ": " << group[i] << std::endl;
          }
          
          out_edge_iterator oi, oi_end;;

          // find least similar neighbour
          double biggest_distance = 0;
          edge_descriptor weakest_link;
          
          for (tie(oi, oi_end) = boost::out_edges(group[i], graph);
               oi != oi_end; ++oi) {
            double current_distance = 
              model->distance(graph[group[i]].state,
                              graph[target(*oi, graph)].state);
            
            if (current_distance > biggest_distance) {
              weakest_link = *oi;
              biggest_distance = current_distance;
            }
          }
          // find most similar neighbour in group
          double lowest_distance = 1.;
          vertex_descriptor similar_neighbour;

          for (unsigned int j = 0; j < group.size(); ++j) {
            double current_distance = 1.;
            if (i != j) {
              current_distance =
                model->distance(graph[group[i]].state,
                                graph[group[j]].state);
            }
            if (current_distance < lowest_distance &&
                edge(group[i], group[j], graph).second == false) {
              similar_neighbour = group[j];
              lowest_distance = current_distance;
            }
          }

          if (lowest_distance < 1 && biggest_distance > 0) {
            if (verbose >=2) {
              std::cout << "removing link to " << target(weakest_link, graph)
                        << std::endl;
            }
            // remove weakest link
            boost::remove_edge(weakest_link, graph);
            // link to group node
            if (verbose >=2) {
              std::cout << "adding link to " << similar_neighbour
                        << std::endl;
            }
            boost::add_edge(group[i], similar_neighbour, graph);
          }
        }
      }
      // seed groups again
      vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        graph[*vi].state->setState(0);
      }
      for (unsigned int i = 1; i <= model->getGroups(); ++i) {
        vertex_descriptor v = random_vertex(graph, randGen);
        if (graph[v].state->getState() == 0) {
          graph[v].state->setState(i);
          active.push_back(v);
          if (verbose >=2) {
            std::cout << "Seeding group " << i << " with vertex " << v << std::endl;
          }
        } else {
          --i;
        }
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
  void GroupFormSimulator<RandomGenerator, Graph>::print()
  {
    std::cout << std::endl;
  }
  
}

#endif
