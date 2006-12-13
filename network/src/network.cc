/******************************************************************/
// main.cc
// contains the main simulation program
/******************************************************************/
#include <iostream>
#include <string>
#include <map>

#include <boost/graph/random.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>

#include <math.h>

#include "Vertex.hh"
#include "GillespieSimulator.hh"
#include "Tree.hh"
#include "RandomGenerator.hh"
#include "lattice.hh"

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

int main(int argc, char* argv[])
{
   /******************************************************************/
   // read parameters
   /******************************************************************/
   std::string topology;
   
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
      ("Sp", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen informed susceptibles")
      ("Sm", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen uninformed susceptibles")
      ("Ip", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen informed infected")
      ("Im", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen uninformed infected")
      ("Rp", po::value<unsigned int>()->default_value(0),
       "number of randomly chosen informed recovered")
      ("Rm", po::value<unsigned int>()->default_value(0),
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
   
   Model m;
   if (vm.count("params_file")) {
      if (m.InitFromFile(vm["params_file"].as<std::string>()) > 0) {
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
   unsigned int steps = vm["steps"].as<unsigned int>();

   if (topology == "lattice") {
      partial_options.add(main_options).add(lattice_options);
      
      /******************************************************************/
      // read lattice specific
      /******************************************************************/

      latticeOptions opt;
      if (vm.count("length")) {
         opt.sideLength = vm["length"].as<unsigned int>();
         std::cout << "Side length: " << opt.sideLength << std::endl;
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
      opt.init.insert(std::make_pair(VertexState(Susceptible, Informed),
                                     vm["Sp"].as<unsigned int>()));
      opt.init.insert(std::make_pair(VertexState(Susceptible, Uninformed),
                                     vm["Sm"].as<unsigned int>()));
      opt.init.insert(std::make_pair(VertexState(Infected, Informed),
                                     vm["Ip"].as<unsigned int>()));
      opt.init.insert(std::make_pair(VertexState(Infected, Uninformed),
                                     vm["Im"].as<unsigned int>()));
      opt.init.insert(std::make_pair(VertexState(Recovered, Informed),
                                     vm["Rp"].as<unsigned int>()));
      opt.init.insert(std::make_pair(VertexState(Recovered, Uninformed),
                                     vm["Rm"].as<unsigned int>()));

      std::string baseString = vm["base"].as<std::string>();
      if (baseString.length() > 1) {
         if (baseString[0] == 'S') {
            opt.base.setDisease(Susceptible);
         } else if (baseString[0] == 'I') {
            opt.base.setDisease(Infected);
         } else if (baseString[0] == 'R') {
            opt.base.setDisease(Recovered);
         }
         if (baseString[1] == 'p') {
            opt.base.setInfo(Informed);
         } else if (baseString[1] == 'm') {
            opt.base.setInfo(Uninformed);
         }
      }

      opt.init[opt.base] = 0;

   
      /******************************************************************/
      // generate lattice with desired properties 
      /******************************************************************/
      boost::mt19937 gen;
      boost::generate_lattice(g,opt.dimensions,opt.sideLength);

      /******************************************************************/
      // set vertex states 
      /******************************************************************/
      
      unsigned int N = static_cast<int>(pow(opt.sideLength, opt.dimensions));

      unsigned int initSum = 0;
      for (std::map<VertexState, unsigned int>::iterator it = opt.init.begin();
           it != opt.init.end(); it++) {
         initSum += (*it).second;
      }
      
      if (initSum > N) {
         std::cout << "Error: number of vertices to select randomly"
                   << " higher than number of total vertices" << std::endl;
      }

      std::cout << "Setting base to " << opt.base << std::endl;
      boost::graph_traits<gillespie_graph>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
         g[*vi].state = opt.base;
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
      gSim->initialize(m);
      generateTree(t,g,get(&Vertex::rateSum, g),get(boost::vertex_index, g));
      std::cout << "time elapsed: " << gSim->getTime() << std::endl;
      boost::print_lattice(g, opt.sideLength);

      /******************************************************************/
      // run simulation
      /******************************************************************/
      for (unsigned int i=0; i<steps; i++) {
         gSim->updateState(m);
         std::cout << "time elapsed: " << gSim->getTime() << std::endl;
         boost::print_lattice(g, opt.sideLength);
      }
      
   } else {
      std::cout << "ERROR: unknown topology: " << topology << std::endl;
      std::cout << std::endl;
      std::cout << main_options << std::endl;
   }      

   return 0;
}
