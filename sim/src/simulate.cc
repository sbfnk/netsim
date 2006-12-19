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
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <math.h>

#include "Vertex.hh"
#include "GillespieSimulator.hh"
#include "Tree.hh"
#include "lattice_generator.hh"
#include "graph_structure.hh"
#include "erdos_renyi_generator2.hh"

namespace po = boost::program_options;

struct latticeOptions {

      latticeOptions()
         : sideLength(0)
      {;}

      unsigned int sideLength;
      unsigned int dimensions;
      bool periodicBoundary;
};

struct rgOptions {
      unsigned int edges;
};

template <typename Graph>
void print_graph_statistics(Graph& g,
                            std::vector<VertexState>& possibleStates,
                            std::vector<EdgeType>& possibleEdgeTypes)
{
   std::map<VertexState, unsigned int> stateCounts;
   std::map<std::pair<VertexState, VertexState>, unsigned int> edgeCounts;

   std::cout << std::endl;
   std::cout << "Vertex count: " << std::endl;

   // generate vertex map
   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      stateCounts.insert(std::make_pair(*it,0));
   }

   // count vertices
   boost::graph_traits<gillespie_graph>::vertex_iterator vi, vi_end;
   for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      stateCounts[g[*vi].state]++;
   }

   // print statistics
   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      if (stateCounts[*it] > 0) {
         std::cout << *it << ": " << stateCounts[*it] << std::endl;
      }
   }
   
   std::cout << std::endl;
   std::cout << "Edge count: " << std::endl;

   // generate edge maps
   for (std::vector<EdgeType>::iterator etIt = possibleEdgeTypes.begin();
        etIt != possibleEdgeTypes.end(); etIt++) {
      edgeCounts.clear();
      std::cout << *etIt << "-type: " << std::endl;

      for (std::vector<VertexState>::iterator it = possibleStates.begin();
           it != possibleStates.end(); it++) {
         for (std::vector<VertexState>::iterator it2 = possibleStates.begin();
              it2 != possibleStates.end(); it2++) {
            edgeCounts.insert(std::make_pair(std::make_pair(*it,*it2),0));
         }
      }

      // count edges
      boost::graph_traits<gillespie_graph>::edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
         if (g[*ei].type == *etIt)
            edgeCounts[std::make_pair(g[source(*ei, g)].state,
                                      g[target(*ei, g)].state)]++;
      }

      for (std::vector<VertexState>::iterator it = possibleStates.begin();
           it != possibleStates.end(); it++) {
         for (std::vector<VertexState>::iterator it2 = it;
              it2 != possibleStates.end(); it2++) {
            unsigned int nEdges;
            if (it == it2) {
               nEdges = edgeCounts[std::make_pair(*it, *it2)];
            } else {
               nEdges = edgeCounts[std::make_pair(*it, *it2)] +
                  edgeCounts[std::make_pair(*it2, *it)];
            }
            if (nEdges > 0) {
               std::cout << *it << *it2 << ": " << nEdges << std::endl;
            }
         }
      }
   }

   std::cout << std::endl;
}

int main(int argc, char* argv[])
{
   Model model;
   std::vector<VertexState> possibleStates = model.getPossibleStates();
   std::vector<EdgeType> possibleEdgeTypes = model.getPossibleEdgeTypes();

   std::map<VertexState, unsigned int> init;
   VertexState base;

   unsigned int lattice_sideLength = 0;

   /******************************************************************/
   // read parameters
   /******************************************************************/
   std::string topology;
   unsigned int N;
   double stop;
   unsigned int outputSteps;
   
   po::options_description main_options
      ("Usage: simulate -p params_file [options]... \n\nAllowed options");

   po::positional_options_description pd;
   pd.add("model", 1);
   
   main_options.add_options()
      ("help,h",
       "produce help message")
      ("longhelp",
       "produce long help message including all topology-related options")
      ("params_file,p",po::value<std::string>(),
       "file containing model parameters")
      ("vertices,n", po::value<unsigned int>(),
       "number of vertices")
      ("topology, t", po::value<std::string>(),
       "network topology\n(lattice,random)")
//       ("i-topology, i", po::value<std::string>(),
//        "information network topology\n(lattice,random)")
      ("stop,s", po::value<double>()->default_value(100.),
       "time after which to stop")
      ("output,o", po::value<unsigned int>()->default_value(1),
       "display status every N steps (0 for display only at start and end")
      ("base,b", po::value<std::string>()->default_value("S"),
       "base state of individuals\n(S,s,I,i,R,r)")
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

   po::options_description lattice_options
      ("Lattice options");
   lattice_options.add_options()
      ("length,l", po::value<unsigned int>(),
       "side length of lattice")
      ("dim",  po::value<unsigned int>()->default_value(2),
       "number of dimensions")
      ("pb", "periodic boundary conditions")
      ;

   po::options_description rg_options
      ("Random graph options");
   rg_options.add_options()
      ("edges,e",  po::value<unsigned int>(),
       "number of edges")
      ;

   po::options_description all_options;
   all_options.add(main_options).add(lattice_options).add(rg_options);
   
   po::options_description partial_options;

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, all_options), vm);
   po::notify(vm);

   if (vm.count("help")) {
      std::cout << main_options << std::endl;
      return 1;
   }

   if (vm.count("longhelp")) {
      std::cout << all_options << std::endl;
      return 1;
   }

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
      topology = vm["topology"].as<std::string>();
   } else {
      std::cout << "ERROR: no topology specified"
                << std::endl;
      std::cout << std::endl;
      std::cout << main_options << std::endl;
      return 1;
   }

   boost::mt19937 gen(time(0));
   GillespieSimulator<boost::mt19937>* gSim =
      new GillespieSimulator<boost::mt19937>(gen);
   gillespie_graph& g = gSim->graph;
   Tree<unsigned int>& t = gSim->tree;

   stop = vm["stop"].as<double>();
   outputSteps = vm["output"].as<unsigned int>();
   
   for (std::vector<VertexState>::iterator it = possibleStates.begin();
        it != possibleStates.end(); it++) {
      std::stringstream ss;
      ss << "random-" << (*it).getString();
      std::string s(ss.str());
      if (vm.count(s.c_str())) {
            unsigned int random = vm[s.c_str()].as<unsigned int>();
            init.insert(std::make_pair(*it,random));
      }
   }
   
   base.set(vm["base"].as<std::string>());
   init[base] = 0;
   
   if (vm.count("vertices")) {
     N = vm["vertices"].as<unsigned int>();
   } else {
     std::cout << "ERROR: no number of vertices specified" << std::endl;
     std::cout << std::endl;
     std::cout << main_options << std::endl;
     return 1;
   }
   
   if (topology == "lattice") {
      partial_options.add(main_options).add(lattice_options);

      /******************************************************************/
      // read lattice specific parameters
      /******************************************************************/

      latticeOptions opt;
      opt.dimensions = vm["dim"].as<unsigned int>();

      opt.sideLength = static_cast<int>(pow(N, 1.0/opt.dimensions));
      if (pow(opt.sideLength, opt.dimensions) != N) { 
        std::cout << "WARNING: lattice resized to " << N << " vertices to "
                  << "make it a square lattice" << std::endl;
      }
      lattice_sideLength = opt.sideLength;

      if (vm.count("pb")) {
         opt.periodicBoundary = true;
      } else {
         opt.periodicBoundary = false;
      }

      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      boost::add_vertices(g, N, Vertex(base));
      typedef boost::lattice_iterator<gillespie_graph> lattice_iterator;

      lattice_iterator li(opt.sideLength,opt.dimensions,opt.periodicBoundary);
      lattice_iterator li_end;
      
      boost::add_edge_structure(g, li, li_end, Edge(Disease));
      boost::add_edge_structure(g, li, li_end, Edge(Information));

   } else if (topology == "random") {

      partial_options.add(main_options).add(rg_options);

      /******************************************************************/
      // read random graph specific parameters
      /******************************************************************/

      rgOptions opt;

      if (vm.count("edges")) {
         opt.edges = vm["edges"].as<unsigned int>();
      } else {
         std::cout << "ERROR: no number of edges specified" << std::endl;
         std::cout << std::endl;
         std::cout << partial_options << std::endl;
         return 1;
      }

      /******************************************************************/
      // generate random graph with desired properties
      /******************************************************************/
      boost::add_vertices(g, N, Vertex(base));
      typedef boost::erdos_renyi_iterator2<boost::mt19937, gillespie_graph>
         rg_iterator;

      double p = (double)opt.edges*2/(double)(N*(N-1));
      rg_iterator ri(gen,N, p);
      rg_iterator ri_end;

      boost::add_edge_structure(g, ri, ri_end, Edge(Disease));
      boost::add_edge_structure(g, ri, ri_end, Edge(Information));

   } else {
      std::cout << "ERROR: unknown topology: " << topology << std::endl;
      std::cout << std::endl;
      std::cout << main_options << std::endl;
      return 1;
   }

   /******************************************************************/
   // set vertex states
   /******************************************************************/
   
   unsigned int initSum = 0;
   for (std::map<VertexState, unsigned int>::iterator it =
           init.begin();
        it != init.end(); it++) {
      initSum += (*it).second;
   }
   
   if (initSum > N) {
      std::cout << "Error: number of vertices to select randomly"
                << " higher than number of total vertices" << std::endl;
   }
   
   boost::graph_traits<gillespie_graph>::vertex_descriptor v;
   for (std::map<VertexState, unsigned int>::iterator it =
           init.begin();
        it != init.end(); it++) {
      for (unsigned int i=0; i<(*it).second; i++) {
         bool inserted = false;
         while (!inserted) {
            v = boost::random_vertex(g, gen);
            if (g[v].state == base) {
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
//    boost::print_lattice(g, lattice_sideLength);
   print_graph_statistics(g, possibleStates, possibleEdgeTypes);
   
   /******************************************************************/
   // run simulation
   /******************************************************************/
   unsigned int steps = 1;
   while (gSim->getTime()<stop && gSim->updateState(model)) {
      if ((outputSteps > 0) && (steps%outputSteps == 0)) {
         std::cout << "time elapsed: " << gSim->getTime() << std::endl;
//          boost::print_lattice(g, lattice_sideLength);
      }
      ++steps;
   }
   
   std::cout << "Final status:" << std::endl;
//    boost::print_lattice(g, lattice_sideLength);
   print_graph_statistics(g, possibleStates, possibleEdgeTypes);

   return 0;
}
