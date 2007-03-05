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

#include "lattice_generator.hh"
#include "graph_structure.hh"
#include "erdos_renyi_generator2.hh"
#include "visualize_graph.hh"

#include "GillespieSimulator.hh"
#include "ChrisSimulator.hh"
#include "Vertex.hh"
#include "InfoSIRS.hh"

namespace po = boost::program_options;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> dualtype_graph;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> onetype_graph;

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
  InfoSIRS model;

  std::vector <unsigned int> init;

  unsigned int outputNum = 0;

  /******************************************************************/
  // read parameters
  /******************************************************************/
  std::string topology;
  
  unsigned int N;
  double stop;

  double outputData = 0.;
  int outputGraphviz = 0;
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
    ("sim", po::value<std::string>()->default_value("Gillespie"),
     "simulator to use (Gillespie, Chris)")
    ("tmax", po::value<double>()->default_value(100.),
     "time after which to stop")
    ("output", po::value<double>()->default_value(0.),
     "write output data at arg timesteps")
    ("graphviz,g", po::value<int>()->default_value(0),
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
    ("vertices,N", po::value<unsigned int>(),
     "number of vertices")
    ("d-topology", po::value<std::string>(),
     "disease network topology\n((tri-)lattice,random,small-world,scale-free,complete)")
    ("i-topology", po::value<std::string>(),
     "information network topology\n((tri-)lattice,random,small-world,scale-free,complete,copy)")
    ("base,b", po::value<std::string>()->default_value
     (model.getVertexStates().begin()->getText()),
     "base state of individuals")
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
  std::vector<po::options_description*> lattice_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i] << "-Lattice Options";
    po::options_description* lo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i] << "-dim";
    lo->add_options()
      (s.str().c_str(),  po::value<unsigned int>()->default_value(2),
       "number of dimensions");
    s.str("");
    s << model.getEdgeTypes()[i] << "-pb";
    lo->add_options()
      (s.str().c_str(), "periodic boundary conditions");
    lattice_options.push_back(lo);
  }

  // random graph
  std::vector<po::options_description*> rg_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i] << "-RandomGraph Options";
    po::options_description* ro
      = new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i] << "-edges";
    ro->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of edges");
    rg_options.push_back(ro);
  }

  // small-world graph
  std::vector<po::options_description*> sw_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i] << "-SmallWorld Options";
    po::options_description* swo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i] << "-neighbours";
    swo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of neighbours of each node");
    s.str("");
    s << model.getEdgeTypes()[i] << "-rewiring-prob";
    swo->add_options()
      (s.str().c_str(), po::value<double>(),
       "rewiring probability");
    sw_options.push_back(swo);
  }

  // scale-free graph
  std::vector<po::options_description*> sf_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i] << "-ScaleFree Options";
    po::options_description* sfo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i] << "-alpha";
    sfo->add_options()
      (s.str().c_str(), po::value<double>(),
       "alpha (index of power law)");
    s.str("");
    s << model.getEdgeTypes()[i] << "-beta";
    sfo->add_options()
      (s.str().c_str(), po::value<double>(),
       "beta (multiplicative factor of power law)");
    sf_options.push_back(sfo);
  }

  po::options_description model_options = model.getOptions();

  // read options from command line
  po::options_description visible_options;
  visible_options.add(command_line_options).add(main_options).add(graph_options).
    add(model_options);
  
  po::options_description all_options;
  all_options.add(visible_options).add(hidden_option);
  
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    all_options.add(*(lattice_options[i]));
    all_options.add(*(rg_options[i]));
    all_options.add(*(sw_options[i]));
    all_options.add(*(sf_options[i]));
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

  /******************************************************************/
  // initialize model
  /******************************************************************/

  model.Init(vm);

  
  /******************************************************************/
  // create simulator
  /******************************************************************/

  unsigned int seed;
  struct timeval tv;
  gettimeofday(&tv, 0);
  seed = tv.tv_sec + tv.tv_usec;
  boost::mt19937 gen(seed);
  dualtype_graph graph;
  Simulator* sim;

  if (vm.count("sim")) {
    std::string simType = vm["sim"].as<std::string>();
    if (simType == "Gillespie") {
      sim = new GillespieSimulator<boost::mt19937, dualtype_graph>
        (gen, graph, model);
    } else if (simType == "Chris") {
      sim = new ChrisSimulator<boost::mt19937, dualtype_graph>
        (gen, graph, model);
    } else {
      std::cerr << "Error: unknown simulator: " << simType << std::endl;
      return 1;
    }
  } else {
    std::cerr << "Error: no simulator specified" << std::endl;
    return 1;
  }
  
  /******************************************************************/
  // read simulation paramters
  /******************************************************************/

  stop = vm["tmax"].as<double>();
  outputData = vm["output"].as<double>();
  if (vm.count("graphviz")) {
    outputGraphviz = vm["graphviz"].as<int>();
  }
  if (vm.count("graph-dir")) {
    graphDir = vm["graph-dir"].as<std::string>();
  }
  // no-graph overrides other graph options
  if (vm.count("no-graph")) {
    outputGraphviz = -1;
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
   
  /******************************************************************/
  // set initial vertex states
  /******************************************************************/

  for (unsigned int i = 0; i < model.getVertexStates().size(); i++) {
    std::stringstream ss;
    ss << model.getVertexStates()[i].getText();
    std::string s(ss.str());
    unsigned int random = 0;
    if (vm.count(s.c_str())) {
      random = vm[s.c_str()].as<unsigned int>();
    }
    init.push_back(random);
  }


  std::string baseString = vm["base"].as<std::string>();
  unsigned int baseState = 0;
  while (baseState < model.getVertexStates().size() &&
         (model.getVertexStates()[baseState].getText() != baseString)) {
    baseState++;
  }
  if (baseState < model.getVertexStates().size()) {
    init[baseState] = 0;
  } else {
    std::cerr << "ERROR: no unknown base state: " << baseString << std::endl;
    std::cerr << std::endl;
    std::cerr << command_line_options << main_options << std::endl;
    return 1;
  }
   
  if (vm.count("vertices")) {
    N = vm["vertices"].as<unsigned int>();
  } else {
    std::cerr << "ERROR: no number of vertices specified" << std::endl;
    std::cerr << std::endl;
    std::cerr << command_line_options << main_options << std::endl;
    return 1;
  }

  /******************************************************************/
  // generate vertices
  /******************************************************************/

  boost::add_vertices(graph, N, Vertex(baseState));
  unsigned int initSum = 0;
  for (std::vector<unsigned int>::iterator it = init.begin();
       it != init.end(); it++) {
    initSum += (*it);
  }
   
  if (initSum > N) {
    std::cerr << "Error: number of vertices to select randomly"
              << " higher than number of total vertices" << std::endl;
  }
   
  boost::graph_traits<dualtype_graph>::vertex_descriptor v;
  for (unsigned int i=0; i<model.getVertexStates().size(); i++) {
    for (unsigned int j=0; j<init[i]; j++) {
      bool inserted = false;
      while (!inserted) {
        v = boost::random_vertex(graph, gen);
        if (graph[v].state == baseState) {
          graph[v].state = i;
          inserted = true;
        }
      }
    }
  }
   

  /******************************************************************/
  // generate edges
  /******************************************************************/

  for (unsigned int i = 0; i < model.getEdgeTypes().size() ; i++) {

    onetype_graph temp_graph;
    boost::add_vertices(temp_graph, N);
    
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-topology";
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
      s << model.getEdgeTypes()[i] << "-dim";
      opt.dimensions = vm[s.str()].as<unsigned int>();
      
      opt.sideLength = static_cast<int>(pow(N, 1.0/opt.dimensions));
      if (pow(opt.sideLength, opt.dimensions) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }

      s.str("");
      s << model.getEdgeTypes()[i] << "-pb";
      if (vm.count(s.str())) {
        opt.periodicBoundary = true;
      } else {
        opt.periodicBoundary = false;
      }

      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      typedef boost::lattice_iterator<dualtype_graph> lattice_iterator;
      
      lattice_iterator li(opt.sideLength,opt.dimensions,opt.periodicBoundary);
      lattice_iterator li_end;
      
      boost::add_edge_structure(temp_graph, li, li_end, Edge(i));
      
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

      s << model.getEdgeTypes()[i] << "-pb";
      if (vm.count(s.str())) {
        opt.periodicBoundary = true;
      } else {
        opt.periodicBoundary = false;
      }

      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      typedef boost::tri_lattice_iterator<dualtype_graph> tri_lattice_iterator;
      
      tri_lattice_iterator tli(opt.sideLength,opt.periodicBoundary);
      tri_lattice_iterator tli_end;
      
      boost::add_edge_structure(temp_graph, tli, tli_end, Edge(i));

    } else if (topology == "random") {

      /******************************************************************/
      // read random graph specific parameters
      /******************************************************************/
      
      rgOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i] << "-edges";
      if (vm.count(s.str())) {
        opt.edges = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << lattice_options[i] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate random graph with desired properties
      /******************************************************************/
      typedef boost::erdos_renyi_iterator2<boost::mt19937, dualtype_graph>
        rg_iterator;
      
      double p = (double)opt.edges*2/(double)(N*(N-1));
      rg_iterator ri(gen,N, p);
      rg_iterator ri_end;
      
      boost::add_edge_structure(temp_graph, ri, ri_end, Edge(i));
      
    } else if (topology == "small-world") {

      /******************************************************************/
      // read small-world graph specific parameters
      /******************************************************************/
      
      swOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i] << "-neighbours";
      if (vm.count(s.str())) {
        opt.neighbours = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of neighbours specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << sw_options[i] << std::endl;
        return 1;
      }
      s.str("");
      s << model.getEdgeTypes()[i] << "-rewiring-prob";
      if (vm.count(s.str())) {
        opt.rewiringProb = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no rewiring probability" << std::endl;
        std::cerr << std::endl;
        std::cerr << sw_options[i] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate small-world graph with desired properties
      /******************************************************************/
      typedef boost::small_world_iterator<boost::mt19937, dualtype_graph>
        sw_iterator;

      sw_iterator swi(gen,N,opt.neighbours, opt.rewiringProb);
      sw_iterator swi_end;
      
      boost::add_edge_structure(temp_graph, swi, swi_end, Edge(i));
    } else if (topology == "scale-free") {

      /******************************************************************/
      // read scale-free graph specific parameters
      /******************************************************************/
      
      sfOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i] << "-alpha";
      if (vm.count(s.str())) {
        opt.alpha = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no alpha specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << sf_options[i] << std::endl;
        return 1;
      }
      s.str("");
      s << model.getEdgeTypes()[i] << "-beta";
      if (vm.count(s.str())) {
        opt.beta = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no beta specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << sf_options[i] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate scale-free graph with desired properties
      /******************************************************************/
      typedef boost::plod_iterator<boost::mt19937, dualtype_graph>
        sf_iterator;
      
      sf_iterator sfi(gen,N,opt.alpha,opt.beta);
      sf_iterator sfi_end;
      
      boost::add_edge_structure(temp_graph, sfi, sfi_end, Edge(i));
    } else if (topology == "complete") {
      boost::graph_traits<dualtype_graph>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; vi++) {
        boost::graph_traits<dualtype_graph>::vertex_iterator vi2;
        for (vi2 = vi+1; vi2 != vi_end; vi2++) {
          add_edge(*vi, *vi2, Edge(i), temp_graph);
        }
      }
    } else if (topology == "copy") {
      boost::graph_traits<dualtype_graph>::edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(graph); ei != ei_end; ei++) {
        if (graph[*ei].type == 0) {
          add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
                   Edge(i), temp_graph);
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
    boost::graph_traits<dualtype_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(temp_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
               Edge(temp_graph[*ei].type), graph);
    }
  }

  /******************************************************************/
  // initalize Simulator
  /******************************************************************/
  sim->initialize();
  if (verbose) std::cout << "time elapsed: " << sim->getTime() << std::endl;
  if (outputGraphviz >= 0) write_graph(graph,model,(graphDir + "/start"),-1);
  std::string lastLine = "";
  if (outputFile) {
    lastLine = write_graph_data(graph, model, sim->getTime(), *outputFile);
  }
  if (verbose) print_graph_statistics(graph, model);
   
  /******************************************************************/
  // run simulation
  /******************************************************************/
  double nextDataStep = outputData;
  
  unsigned int steps = 0;

//   sim->print();

  while (sim->updateState() && sim->getTime()<stop) {
//     sim->print();
    
    if (verbose && steps%100 == 0) {
      std::cout << "time elapsed: " << sim->getTime() << std::endl;
    }

    if ((outputGraphviz > 0) && (outputGraphviz % steps == 0)) {
      write_graph(graph, model,
                  generateFileName((graphDir +"/frame"),outputNum),
                  sim->getTime());
      ++outputNum;
    }
    if (outputFile && sim->getTime() > nextDataStep) {
      lastLine = 
        write_graph_data(graph, model, sim->getTime(), *outputFile);
      if (outputData > 0) {
        do {
          nextDataStep += outputData;
        } while (sim->getTime() > nextDataStep);
      }
    }
    ++steps;
  }
   
  if (verbose) std::cout << "Final status:" << std::endl;
  if (outputGraphviz >= 0) {
    write_graph(graph, model, (graphDir + "/end"), sim->getTime());
  }
  
  if (outputFile) {
    *outputFile << stop << '\t' << lastLine;
    outputFile->close();
  }
  
  if (verbose) print_graph_statistics(graph, model);

  return 0;
}
