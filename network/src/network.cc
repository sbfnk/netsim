/******************************************************************/
// main.cc
// contains the main simulation program
/******************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

#include <sys/time.h>

#include <boost/graph/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>

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

struct swOptions {
  unsigned int neighbours;
  double rewiringProb;
};

struct sfOptions {
  double alpha;
  double beta;
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

  std::map<VertexState, unsigned int> init;
  VertexState base;

  unsigned int outputNum = 0;

  /******************************************************************/
  // read parameters
  /******************************************************************/
  std::string topology;
  unsigned int N;
  double stop;

  double outputData = 0.;
  double outputGraphviz = 0.;
  std::string graphDir = "images";

  std::string outputFileName = "";
  
  std::ofstream* outputFile = 0;

  bool verbose = false;
   
  po::options_description command_line_options
    ("Usage: simulate -p params_file [options]... \n\nAllowed options");

  command_line_options.add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message including all options")
    ("verbose,v",
     "produce verbose output")
    ("params-file,p",po::value<std::string>(),
     "file containing graph parameters")
    ("model-file,m",po::value<std::string>(),
     "file containing model parameters")
    ;
    
  po::options_description main_options;

  main_options.add_options()
    ("tmax", po::value<double>()->default_value(100.),
     "time after which to stop")
    ("output", po::value<double>()->default_value(0.),
     "write output data at arg timesteps")
    ("graphviz,g", po::value<double>()->default_value(1.),
     "create graphviz output in the images directory at arg timesteps")
    ("graph-dir", po::value<std::string>()->default_value(graphDir),
     "set ouput dir for graphs")
    ("write-file,f", po::value<std::string>(),
     "output data to file (.sim.dat will be appended)")
    ;

  // declare hidden option for suppression of graph output --
  // needed for do_all script so that graphviz output
  // is not generated at each run
  po::options_description hidden_option;
  hidden_option.add_options()
    ("no-graph",
     "do not produce graphviz output no matter what the other settings")
    ;
  
  po::options_description graph_options;
  
  graph_options.add_options()
    ("d-topology", po::value<std::string>(),
     "disease network topology\n((tri-)lattice,random,small-world,scale-free,complete)")
    ("i-topology", po::value<std::string>(),
     "information network topology\n((tri-)lattice,random,small-world,scale-free,complete,copy)")
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
  for (std::vector<EdgeType>::const_iterator etIt =
         model.getPossibleEdgeTypes().begin();
       etIt != model.getPossibleEdgeTypes().end(); etIt++) {
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
  for (std::vector<EdgeType>::const_iterator etIt =
         model.getPossibleEdgeTypes().begin();
       etIt != model.getPossibleEdgeTypes().end(); etIt++) {
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

  // small-world graph
  std::map<EdgeType, po::options_description> sw_options;
  for (std::vector<EdgeType>::const_iterator etIt =
         model.getPossibleEdgeTypes().begin();
       etIt != model.getPossibleEdgeTypes().end(); etIt++) {
    std::stringstream s;
    s << *etIt << "-SmallWorld Options";
    po::options_description swo(s.str().c_str());
    s.str("");
    s << *etIt << "-neighbours";
    swo.add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of neighbours of each node");
    s.str("");
    s << *etIt << "-rewiring-prob";
    swo.add_options()
      (s.str().c_str(), po::value<double>(),
       "rewiring probability");
    sw_options.insert(std::make_pair(*etIt, swo));
  }

  // scale-free graph
  std::map<EdgeType, po::options_description> sf_options;
  for (std::vector<EdgeType>::const_iterator etIt =
         model.getPossibleEdgeTypes().begin();
       etIt != model.getPossibleEdgeTypes().end(); etIt++) {
    std::stringstream s;
    s << *etIt << "-ScaleFree Options";
    po::options_description sfo(s.str().c_str());
    s.str("");
    s << *etIt << "-alpha";
    sfo.add_options()
      (s.str().c_str(), po::value<double>(),
       "alpha (index of power law)");
    s.str("");
    s << *etIt << "-beta";
    sfo.add_options()
      (s.str().c_str(), po::value<double>(),
       "beta (multiplicative factor of power law)");
    sf_options.insert(std::make_pair(*etIt, sfo));
  }

  po::options_description model_options
    ("Model parameters");
  
  model_options.add_options()
    ("vertices,N", po::value<unsigned int>(),
     "number of vertices")
    ("beta--", po::value<double>(),
     "disease transmission rate uninformed->uninformed")
    ("beta+-", po::value<double>(),
     "disease transmission rate informed->uninformed")
    ("beta-+", po::value<double>(),
     "disease transmission rate uninformed->informed")
    ("beta++", po::value<double>(),
     "disease transmission rate informed->informed")
    ("gamma-", po::value<double>(),
     "recovery rate of uninformed")
    ("gamma+", po::value<double>(),
     "recovery rate of informed")
    ("delta-", po::value<double>(),
     "loss of immunity rate of uninformed")
    ("delta+", po::value<double>(),
     "loss of immunity rate of informed")
    ("alpha", po::value<double>(),
     "information transmission rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("omega", po::value<double>(),
     "local information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ;

  // read options from command line
  po::options_description visible_options;
  visible_options.add(command_line_options).add(main_options).add(graph_options).
    add(model_options);
  
  po::options_description all_options;
  all_options.add(visible_options).add(hidden_option);
  
  for (std::vector<EdgeType>::const_iterator etIt =
         model.getPossibleEdgeTypes().begin();
       etIt != model.getPossibleEdgeTypes().end(); etIt++) {
    all_options.add(lattice_options[*etIt]);
    all_options.add(rg_options[*etIt]);
    all_options.add(sw_options[*etIt]);
    all_options.add(sf_options[*etIt]);
  }
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all_options).
            allow_unregistered().run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << command_line_options << main_options << std::endl;
    return 1;
  }

  if (vm.count("longhelp")) {
    std::cout << visible_options << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  if (vm.count("params-file")) {
    std::ifstream ifs(vm["params-file"].as<std::string>().c_str());
    try {
      po::store(po::parse_config_file(ifs, all_options), vm);
    }
    catch (std::exception& e) {
      std::cerr << "Error parsing params file: " << e.what() << std::endl;
      return 1;
    }
  }

  if (vm.count("model-file")) {
    std::ifstream ifs(vm["model-file"].as<std::string>().c_str());
    try {
      po::store(po::parse_config_file(ifs, all_options), vm);
    }
    catch (std::exception& e) {
      std::cerr << "Error parsing model file: " << e.what() << std::endl;
      return 1;
    }
  }

  po::notify(vm);

  model.Init(vm);

  unsigned int seed;
  struct timeval tv;

  gettimeofday(&tv, 0);
  seed = tv.tv_sec + tv.tv_usec;
  
  boost::mt19937 gen(seed);
  GillespieSimulator<boost::mt19937>* gSim =
    new GillespieSimulator<boost::mt19937>(gen);
  gillespie_graph& graph = gSim->graph;
  Tree<unsigned int>& tree = gSim->tree;

  stop = vm["tmax"].as<double>();
  outputData = vm["output"].as<double>();
  if (vm.count("graphviz")) {
    outputGraphviz = vm["graphviz"].as<double>();
  }
  if (vm.count("graph-dir")) {
    graphDir = vm["graph-dir"].as<std::string>();
  }
  // no-graph overrides other graph options
  if (vm.count("no-graph")) {
    outputGraphviz = 0.;
  }

  if (vm.count("write-file")) {
    outputFileName = (vm["write-file"].as<std::string>())+".sim.dat";
    try {
      outputFile = new std::ofstream();
      outputFile->open(outputFileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
    }
  }
   
  for (std::vector<VertexState>::const_iterator it =
         model.getPossibleStates().begin();
       it != model.getPossibleStates().end(); it++) {
    std::stringstream ss;
    ss << (*it).getString();
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

  // generate vertices
  boost::add_vertices(graph, N, Vertex(base));

  // generate edges
  for (std::vector<EdgeType>::const_iterator etIt =
         model.getPossibleEdgeTypes().begin();
       etIt != model.getPossibleEdgeTypes().end(); etIt++) {

    onetype_graph temp_graph;
    boost::add_vertices(temp_graph, N);
    
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
      
      boost::add_edge_structure(temp_graph, li, li_end, Edge(*etIt));
      
    } else if (topology == "tri-lattice") {

      /******************************************************************/
      // read tri-lattice specific parameters
      /******************************************************************/
      
      latticeOptions opt;
      s.str("");

      opt.sideLength = static_cast<int>(sqrt(N));
      if (pow(opt.sideLength, 2) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }

      s << *etIt << "-pb";
      if (vm.count(s.str())) {
        opt.periodicBoundary = true;
      } else {
        opt.periodicBoundary = false;
      }

      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      typedef boost::tri_lattice_iterator<gillespie_graph> tri_lattice_iterator;
      
      tri_lattice_iterator tli(opt.sideLength,opt.periodicBoundary);
      tri_lattice_iterator tli_end;
      
      boost::add_edge_structure(temp_graph, tli, tli_end, Edge(*etIt));

    } else if (topology == "random") {

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
        std::cerr << lattice_options[*etIt] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate random graph with desired properties
      /******************************************************************/
      typedef boost::erdos_renyi_iterator2<boost::mt19937, gillespie_graph>
        rg_iterator;
      
      double p = (double)opt.edges*2/(double)(N*(N-1));
      rg_iterator ri(gen,N, p);
      rg_iterator ri_end;
      
      boost::add_edge_structure(temp_graph, ri, ri_end, Edge(*etIt));
      
    } else if (topology == "small-world") {

      /******************************************************************/
      // read small-world graph specific parameters
      /******************************************************************/
      
      swOptions opt;
      
      s.str("");
      s << *etIt << "-neighbours";
      if (vm.count(s.str())) {
        opt.neighbours = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of neighbours specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << sw_options[*etIt] << std::endl;
        return 1;
      }
      s.str("");
      s << *etIt << "-rewiring-prob";
      if (vm.count(s.str())) {
        opt.rewiringProb = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no rewiring probability" << std::endl;
        std::cerr << std::endl;
        std::cerr << sw_options[*etIt] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate small-world graph with desired properties
      /******************************************************************/
      typedef boost::small_world_iterator<boost::mt19937, gillespie_graph>
        sw_iterator;

      sw_iterator swi(gen,N,opt.neighbours, opt.rewiringProb);
      sw_iterator swi_end;
      
      boost::add_edge_structure(temp_graph, swi, swi_end, Edge(*etIt));
    } else if (topology == "scale-free") {

      /******************************************************************/
      // read scale-free graph specific parameters
      /******************************************************************/
      
      sfOptions opt;
      
      s.str("");
      s << *etIt << "-alpha";
      if (vm.count(s.str())) {
        opt.alpha = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no alpha specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << sf_options[*etIt] << std::endl;
        return 1;
      }
      s.str("");
      s << *etIt << "-beta";
      if (vm.count(s.str())) {
        opt.beta = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no beta specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << sf_options[*etIt] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate scale-free graph with desired properties
      /******************************************************************/
      typedef boost::plod_iterator<boost::mt19937, gillespie_graph>
        sf_iterator;
      
      sf_iterator sfi(gen,N,opt.alpha,opt.beta);
      sf_iterator sfi_end;
      
      boost::add_edge_structure(temp_graph, sfi, sfi_end, Edge(*etIt));
    } else if (topology == "complete") {
      boost::graph_traits<gillespie_graph>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; vi++) {
        boost::graph_traits<gillespie_graph>::vertex_iterator vi2;
        for (vi2 = vi+1; vi2 != vi_end; vi2++) {
          add_edge(*vi, *vi2, Edge(*etIt), temp_graph);
        }
      }
    } else if (topology == "copy") {
      boost::graph_traits<gillespie_graph>::edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(graph); ei != ei_end; ei++) {
        if (graph[*ei].type == *(model.getPossibleEdgeTypes().begin())) {
          add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
                   Edge(*etIt), temp_graph);
        }
      }
    } else {
      std::cerr << "ERROR: unknown " << s.str() << ": " << topology
                << std::endl;
      std::cerr << std::endl;
      std::cerr << main_options << std::endl;
      return 1;
    }

    // copy edges to main graph
    boost::graph_traits<gillespie_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(temp_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
               Edge(temp_graph[*ei].type), graph);
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
    std::cerr << "Error: number of vertices to select randomly"
              << " higher than number of total vertices" << std::endl;
  }
   
  boost::graph_traits<gillespie_graph>::vertex_descriptor v;
  for (std::map<VertexState, unsigned int>::iterator it =
         init.begin();
       it != init.end(); it++) {
    for (unsigned int i=0; i<(*it).second; i++) {
      bool inserted = false;
      while (!inserted) {
        v = boost::random_vertex(graph, gen);
        if (graph[v].state == base) {
          graph[v].state = (*it).first;
          inserted = true;
        }
      }
    }
  }
   
  /******************************************************************/
  // initalize GillespieSimulator
  /******************************************************************/
  gSim->initialize(model);
  generateTree(tree,graph,get(&Vertex::rateSum, graph),
               get(boost::vertex_index, graph));
  if (verbose) std::cout << "time elapsed: " << gSim->getTime() << std::endl;
  if (outputGraphviz > 0) write_graph(graph, (graphDir + "/start"),-1);
  std::string lastLine = "";
  if (outputFile) {
    lastLine = write_graph_data(graph, gSim->getTime(), *outputFile,
                                model.getPossibleStates(),
                                model.getPossibleEdgeTypes());
  }
  if (verbose) print_graph_statistics(graph, model.getPossibleStates(),
                                      model.getPossibleEdgeTypes());
   
  /******************************************************************/
  // run simulation
  /******************************************************************/
  double nextGraphStep = outputGraphviz;
  double nextDataStep = outputData;
  
  unsigned int steps = 0;

  while (gSim->updateState(model) && gSim->getTime()<stop) {
    if (verbose && steps%100 == 0) {
      std::cout << "time elapsed: " << gSim->getTime() << std::endl;
    }

    if ((outputGraphviz > 0) && (gSim->getTime() > nextGraphStep)) {
      write_graph(graph, generateFileName((graphDir +"/frame"),outputNum),
                  gSim->getTime());
      do {
        nextGraphStep += outputGraphviz;
      } while (gSim->getTime() > nextGraphStep);
      ++outputNum;
    }
    if (outputFile && gSim->getTime() > nextDataStep) {
      lastLine = 
        write_graph_data(graph, gSim->getTime(), *outputFile,
                         model.getPossibleStates(), model.getPossibleEdgeTypes());
      if (outputData > 0) {
        do {
          nextDataStep += outputData;
        } while (gSim->getTime() > nextDataStep);
      }
    }
    ++steps;
  }
   
  if (verbose) std::cout << "Final status:" << std::endl;
  if (outputGraphviz > 0) write_graph(graph, (graphDir + "/end"), gSim->getTime());
  if (outputFile) {
    *outputFile << stop << '\t' << lastLine;
    outputFile->close();
  }
  
  if (verbose) print_graph_statistics(graph, model.getPossibleStates(),
                                      model.getPossibleEdgeTypes());

  return 0;
}
