/******************************************************************/
// main.cc
// contains the main simulation program
/******************************************************************/
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

#include <boost/graph/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <math.h>

#include "Vertex.hh"
#include "GillespieSimulator.hh"
#include "Tree.hh"
#include "lattice_generator.hh"
#include "graph_structure.hh"
#include "erdos_renyi_generator2.hh"
#include "visualize_graph.hh"

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

std::string generateFileName(std::string nameBase, unsigned int id)
{
  std::stringstream s;
  s << nameBase << std::setw(3) << std::setfill('0') << id;
  return s.str();
}

int main(int argc, char* argv[])
{
  Model model;
  std::vector<VertexState> possibleStates = model.getPossibleStates();
  std::vector<EdgeType> possibleEdgeTypes = model.getPossibleEdgeTypes();

  std::map<VertexState, unsigned int> init;
  VertexState base;

  unsigned int lattice_sideLength = 0;
  unsigned int outputNum = 0;

  /******************************************************************/
  // read parameters
  /******************************************************************/
  std::string topology;
  unsigned int N;
  double stop;
  unsigned int outputSteps;
   
  po::options_description command_line_options
    ("Usage: simulate -p params_file [options]... \n\nAllowed options");

  command_line_options.add_options()
    ("help,h",
     "produce help message")
    ("config-file,c",po::value<std::string>(),
     "file containing graph confguration")
    ("params-file,p",po::value<std::string>(),
     "file containing model parameters")
    ;
    
  po::options_description main_options;

  main_options.add_options()
    ("vertices,n", po::value<unsigned int>(),
     "number of vertices")
    ("d-topology", po::value<std::string>(),
     "disease network topology\n(lattice,random)")
    ("i-topology", po::value<std::string>(),
     "information network topology\n(lattice,random)")
    ;
  main_options.add_options()
    ("stop,s", po::value<double>()->default_value(100.),
     "time after which to stop")
    ("output,o", po::value<unsigned int>()->default_value(1),
     "display status every N steps (0 for display only at start and end")
    ("base,b", po::value<std::string>()->default_value("S"),
     "base state of individuals\n(S,s,I,i,R,r)")
    ("S+", po::value<unsigned int>()->default_value(0),
     "number of randomly chosen informed susceptibles")
    ("S-", po::value<unsigned int>()->default_value(0),
     "number of randomly chosen uninformed susceptibles")
    ("I+", po::value<unsigned int>()->default_value(0),
     "number of randomly chosen informed infected")
    ("I-", po::value<unsigned int>()->default_value(0),
     "number of randomly chosen uninformed infected")
    ("R+", po::value<unsigned int>()->default_value(0),
     "number of randomly chosen informed recovered")
    ("R-", po::value<unsigned int>()->default_value(0),
     "number of randomly chosen uninformed recovered")
    ;

  // generate topology-specifice options for each type of graph

  // lattice
  std::map<EdgeType, po::options_description> lattice_options;
  for (std::vector<EdgeType>::iterator etIt = possibleEdgeTypes.begin();
       etIt != possibleEdgeTypes.end(); etIt++) {
    std::stringstream s;
    s << *etIt << "-Lattice Options";
    po::options_description lo(s.str().c_str());
    s.str("");
    s << *etIt << "-dim";
    lo.add_options()
      (s.str().c_str(),  po::value<unsigned int>()->default_value(2),
       "number of dimensions");
    s.str("");
    s << *etIt << "-pb";
    lo.add_options()
      (s.str().c_str(), "periodic boundary conditions");
    lattice_options.insert(std::make_pair(*etIt, lo));
  }

  // random graph
  std::map<EdgeType, po::options_description> rg_options;
  for (std::vector<EdgeType>::iterator etIt = possibleEdgeTypes.begin();
       etIt != possibleEdgeTypes.end(); etIt++) {
    std::stringstream s;
    s << *etIt << "-RandomGraph Options";
    po::options_description ro(s.str().c_str());
    s.str("");
    s << *etIt << "-edges";
    ro.add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of edges");
    rg_options.insert(std::make_pair(*etIt, ro));
  }

  po::options_description all_options;
  all_options.add(command_line_options).add(main_options);
  for (std::vector<EdgeType>::iterator etIt = possibleEdgeTypes.begin();
       etIt != possibleEdgeTypes.end(); etIt++) {
    all_options.add(lattice_options[*etIt]);
    all_options.add(rg_options[*etIt]);
  }
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, all_options), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << command_line_options << main_options << std::endl;
    return 1;
  }

  if (vm.count("config-file")) {
    po::options_description config_file_options;
    std::ifstream ifs(vm["config-file"].as<std::string>().c_str());
    config_file_options.add(main_options);
    for (std::vector<EdgeType>::iterator etIt = possibleEdgeTypes.begin();
         etIt != possibleEdgeTypes.end(); etIt++) {
      config_file_options.add(lattice_options[*etIt]);
      config_file_options.add(rg_options[*etIt]);
    }
    po::store(po::parse_config_file(ifs, config_file_options), vm);
  }

  po::notify(vm);

  if (vm.count("params-file")) {
    if (model.InitFromFile(vm["params-file"].as<std::string>()) > 0) {
      std::cerr << "ERROR: could not read params-file" << std::endl;
      return 1;
    }
  } else {
    std::cerr << "ERROR: missing params_file" << std::endl;
    std::cerr << std::endl;
    std::cerr << command_line_options << main_options << std::endl;
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
    std::cerr << "ERROR: no number of vertices specified" << std::endl;
    std::cerr << std::endl;
    std::cerr << command_line_options << main_options << std::endl;
    return 1;
  }

  // generate the two networks
  boost::add_vertices(g, N, Vertex(base));

  for (std::vector<EdgeType>::iterator etIt = possibleEdgeTypes.begin();
       etIt != possibleEdgeTypes.end(); etIt++) {
    po::options_description partial_options;
    partial_options.add(command_line_options).add(main_options);
    
    std::stringstream s;
    s << *etIt << "-topology";
    if (vm.count(s.str())) {
      topology = vm[s.str()].as<std::string>();
    } else {
      std::cerr << "ERROR: no " << s.str() << " specified" << std::endl;
      std::cerr << std::endl;
      std::cerr << command_line_options << main_options << std::endl;
      return 1;
    }

    if (topology == "lattice") {
      partial_options.add(lattice_options[*etIt]);
      
      /******************************************************************/
      // read lattice specific parameters
      /******************************************************************/
      
      latticeOptions opt;
      s.str("");
      s << *etIt << "-dim";
      opt.dimensions = vm[s.str()].as<unsigned int>();
      
      opt.sideLength = static_cast<int>(pow(N, 1.0/opt.dimensions));
      if (pow(opt.sideLength, opt.dimensions) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }
      lattice_sideLength = opt.sideLength;

      s.str("");
      s << *etIt << "-pb";
      if (vm.count(s.str())) {
        opt.periodicBoundary = true;
      } else {
        opt.periodicBoundary = false;
      }

      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      typedef boost::lattice_iterator<gillespie_graph> lattice_iterator;
      
      lattice_iterator li(opt.sideLength,opt.dimensions,opt.periodicBoundary);
      lattice_iterator li_end;
      
      boost::add_edge_structure(g, li, li_end, Edge(*etIt));
      
    } else if (topology == "random") {

      partial_options.add(rg_options[*etIt]);
      
      /******************************************************************/
      // read random graph specific parameters
      /******************************************************************/
      
      rgOptions opt;
      
      s.str("");
      s << *etIt << "-edges";
      if (vm.count(s.str())) {
        opt.edges = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << partial_options << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate random graph with desired properties
      /******************************************************************/
      typedef boost::erdos_renyi_iterator2<boost::mt19937, gillespie_graph>
        rg_iterator;
      
      double p = (double)opt.edges*2/(double)(N*(N-1));
      std::cout << "N: " << N << "p: " << p << std::endl;
      rg_iterator ri(gen,N, p);
      rg_iterator ri_end;
      
      boost::add_edge_structure(g, ri, ri_end, Edge(*etIt));
      
    } else {
      std::cerr << "ERROR: unknown " << s.str() << ": " << topology << std::endl;
      std::cerr << std::endl;
      std::cerr << main_options << std::endl;
      return 1;
    }
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
  write_graph(g, "images/start");
  print_graph_statistics(g, possibleStates, possibleEdgeTypes);
   
  /******************************************************************/
  // run simulation
  /******************************************************************/
  unsigned int steps = 1;
  while (gSim->getTime()<stop && gSim->updateState(model)) {
    if ((outputSteps > 0) && (steps%outputSteps == 0)) {
      std::cout << "time elapsed: " << gSim->getTime() << std::endl;
      write_graph(g, generateFileName("images/frame", outputNum));
      ++outputNum;
    }
    ++steps;
  }
   
  std::cout << "Final status:" << std::endl;
  write_graph(g, "images/end");
  print_graph_statistics(g, possibleStates, possibleEdgeTypes);

  return 0;
}
