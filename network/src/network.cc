/******************************************************************/
// main.cc
// contains the main simulation program
/******************************************************************/
#include <iostream>
#include <string>
#include <map>

#include <boost/graph/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>

#include <math.h>

#include "Vertex.hh"
#include "GillespieSimulator.hh"
#include "Tree.hh"
#include "RandomGenerator.hh"
#include "lattice_generator.hh"
#include "generate_graph.hh"

namespace po = boost::program_options;

struct latticeOptions {

      latticeOptions()
         : sideLength(0),
           dimensions(0)
      {;}
      
      unsigned int sideLength;
      unsigned int dimensions;
      VertexState base;
      std::map<VertexState, unsigned int> init;
};

std::vector<VertexState> possibleStates;

template <typename Graph>
void print_graph_statistics(Graph& g)
{
   std::map<VertexState, unsigned int> stateCounts;
   std::map<std::pair<VertexState, VertexState>, unsigned int> edgeCounts;

   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      stateCounts.insert(std::make_pair(*it,0));
   }
   
   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      for (std::vector<VertexState>::iterator it2 = possibleStates.begin();
        it2 != possibleStates.end(); it2++) {
         edgeCounts.insert(std::make_pair(std::make_pair(*it,*it2),0));
      }
   }

   // count vertices
   
   boost::graph_traits<gillespie_graph>::vertex_iterator vi, vi_end;
   for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      stateCounts[g[*vi].state]++;
   }

   // count edges
   boost::graph_traits<gillespie_graph>::edge_iterator ei, ei_end;
   for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      edgeCounts[std::make_pair(g[source(*ei, g)].state,
                                g[target(*ei, g)].state)]++;
   }

   // print statistics
   std::cout << std::endl;
   
   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      if (stateCounts[*it] > 0) {
         std::cout << *it << ": " << stateCounts[*it] << std::endl;
      }
   }
   
   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      for (std::vector<VertexState>::iterator it2 = possibleStates.begin();
        it2 != possibleStates.end(); it2++) {
         if (edgeCounts[std::make_pair(*it, *it2)] > 0) {
            std::cout << *it << *it2 << ": "
                      << edgeCounts[std::make_pair(*it,*it2)] << std::endl;
         }
      }
   }
}

int main(int argc, char* argv[])
{
   possibleStates.push_back(VertexState(Susceptible, Informed));
   possibleStates.push_back(VertexState(Susceptible, Uninformed));
   possibleStates.push_back(VertexState(Infected, Informed));
   possibleStates.push_back(VertexState(Infected, Uninformed));
   possibleStates.push_back(VertexState(Recovered, Informed));
   possibleStates.push_back(VertexState(Recovered, Uninformed));

   /******************************************************************/
   // read parameters
   /******************************************************************/
   std::string topology;
   unsigned int steps;
   unsigned int outputSteps;
   
   po::options_description main_options
    ("Usage: simulate -p params_file [options]... \n\nAllowed options:");

   po::positional_options_description pd;
   pd.add("model", 1);

   main_options.add_options()
      ("help,h",
       "produce help message")
      ("hh",
       "produce long help message including all topology-related options")
      ("params_file,p",po::value<std::string>(),
       "file containing model parameters")
      ("topology,t", po::value<std::string>(),
       "network topology\n(lattice,random)")
      ("steps,s", po::value<unsigned int>()->default_value(100),
       "number of steps")
      ("output,o", po::value<unsigned int>()->default_value(1),
       "display status every N steps (0 for status display only at beginning and end")
      ;

   po::options_description lattice_options
      ("Lattice options:");
   lattice_options.add_options()
      ("length,l", po::value<unsigned int>(),
       "side length of lattice")
      ("dim,d",  po::value<unsigned int>(),
       "number of dimensions")
      ("base,b", po::value<std::string>()->default_value("Sm"),
       "base state of individuals\n(Sp,Sm,Ip,Im,Rp,Rm)")
      ("random-s", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen informed susceptibles")
      ("random-S", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen uninformed susceptibles")
      ("random-i", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen informed infected")
      ("random-I", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen uninformed infected")
      ("random-r", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen informed recovered")
      ("random-R", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen uninformed recovered")
      ;

   po::options_description all_options;
   all_options.add(main_options).add(lattice_options);
   
   po::options_description partial_options;
   
   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, all_options), vm);
   po::notify(vm);    
   
   if (vm.count("help")) {
      std::cout << main_options << std::endl;
      return 1;
   }
   
   if (vm.count("hh")) {
      std::cout << all_options << std::endl;
      return 1;
   }
   
   Model model;
   if (vm.count("params_file")) {
      if (model.InitFromFile(vm["params_file"].as<std::string>()) > 0) {
         std::cout << "ERROR: could not read params_file" << std::endl;
         std::cout << std::endl;
         std::cout << main_options << std::endl;
         return 1;
      }
   } else {
      std::cout << "ERROR: missing params_file" << std::endl;
      std::cout << std::endl;
      std::cout << main_options << std::endl;
      return 1;
   }
   
   if (vm.count("topology") == 0) {
      std::cout << "ERROR: no topology specified" << std::endl;
      std::cout << std::endl;
      std::cout << main_options << std::endl;
      return 1;
   } else {
      topology = vm["topology"].as<std::string>();
   }

   GillespieSimulator* gSim = new GillespieSimulator();
   gillespie_graph& g = gSim->graph;
   Tree<unsigned int>& t = gSim->tree;

   steps = vm["steps"].as<unsigned int>();
   outputSteps = vm["output"].as<unsigned int>();

   if (topology == "lattice") {
      partial_options.add(main_options).add(lattice_options);
      
      /******************************************************************/
      // read lattice specific parameters
      /******************************************************************/

      latticeOptions opt;
      if (vm.count("length")) {
         opt.sideLength = vm["length"].as<unsigned int>();
      } else {
         std::cout << "ERROR: no length specified" << std::endl;
         std::cout << std::endl;
         std::cout << partial_options << std::endl;
         return 1;
      }
      if (vm.count("dim")) {
         opt.dimensions = vm["dim"].as<unsigned int>();
      } else {
         std::cout << "ERROR: no number of dimensions specified" << std::endl;
         std::cout << std::endl;
         std::cout << partial_options << std::endl;
         return 1;
      }

      for (std::vector<VertexState>::iterator it = possibleStates.begin();
           it != possibleStates.end(); it++) {
         std::stringstream ss;
         ss << "random-" << (*it).getString();
         std::string s(ss.str());
         if (vm.count(s.c_str())) {
            opt.init.insert(std::make_pair(*it,
                                           vm[s.c_str()].as<unsigned int>()));
         }
      }

      opt.base.set(vm["base"].as<std::string>());
      opt.init[opt.base] = 0;

   
      /******************************************************************/
      // generate lattice with desired properties 
      /******************************************************************/
      boost::mt19937 gen(time(0));
      unsigned int N = static_cast<int>(pow(opt.sideLength, opt.dimensions));

      boost::add_vertices(g, N, Vertex(opt.base));
      typedef boost::lattice_iterator<gillespie_graph> lattice_iterator;
      boost::add_edge_structure(g, lattice_iterator(opt.sideLength,
                                                    opt.dimensions, false),
                                lattice_iterator(), Edge(Disease));
      boost::add_edge_structure(g, lattice_iterator(opt.sideLength,
                                                    opt.dimensions, false),
                                lattice_iterator(), Edge(Information));

      /******************************************************************/
      // set vertex states 
      /******************************************************************/
      
      unsigned int initSum = 0;
      for (std::map<VertexState, unsigned int>::iterator it = opt.init.begin();
           it != opt.init.end(); it++) {
         initSum += (*it).second;
      }
      
      if (initSum > N) {
         std::cout << "Error: number of vertices to select randomly"
                   << " higher than number of total vertices" << std::endl;
      }

      boost::graph_traits<gillespie_graph>::vertex_descriptor v;
      for (std::map<VertexState, unsigned int>::iterator it = opt.init.begin();
           it != opt.init.end(); it++) {
         for (unsigned int i=0; i<(*it).second; i++) {
            bool inserted = false;
            while (!inserted) {
               v = boost::random_vertex(g, gen);
               if (g[v].state == opt.base) {
                  g[v].state = (*it).first;
                  inserted = true;
               }
            }
         }
      }
      
      /******************************************************************/
      // initalize GillespieSimulator
      /******************************************************************/
      gSim->initialize(model);
      generateTree(t,g,get(&Vertex::rateSum, g),get(boost::vertex_index, g));
      std::cout << "time elapsed: " << gSim->getTime() << std::endl;
      boost::print_lattice(g, opt.sideLength);
      print_graph_statistics(g);

      /******************************************************************/
      // run simulation
      /******************************************************************/
      unsigned int i=0;
      while (i<steps && gSim->updateState(model)) {
         if ((outputSteps > 0) && ((i+1)%outputSteps == 0)) {
            std::cout << "time elapsed: " << gSim->getTime() << std::endl;
            boost::print_lattice(g, opt.sideLength);
         }
         ++i;
      }
      
      std::cout << "Final status:" << std::endl;
      boost::print_lattice(g, opt.sideLength);
      print_graph_statistics(g);
   } else {
      std::cout << "ERROR: unknown topology: " << topology << std::endl;
      std::cout << std::endl;
      std::cout << main_options << std::endl;
   }      

   return 0;
}
