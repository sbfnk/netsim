/*! \file simulate.cc
  \brief The main simulation program.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

#include <sys/time.h>
#include <sys/stat.h>

#include <boost/graph/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>
// #include <boost/graph/erdos_renyi_generator.hpp>
#include <math.h>

#include "graph_structure.hh"
#include "graph_statistics.hh"
#include "lattice_generator.hh"
#include "erdos_renyi_generator2.hh"
#include "albert_barabasi_generator.hh"
#include "visualise_graph.hh"
#include "cluster_coeffs.hh"
#include "assortativity.hh"

#include "GillespieSimulator.hh"
#include "ChrisSimulator.hh"
#include "Vertex.hh"

// models
#include "InfoSIRS.hh"
#include "ProtectiveSIRS.hh"
#include "VaccinationSIRS.hh"

namespace po = boost::program_options;

//! Classes and functions of the boost libraries.
namespace boost {}

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

struct rrgOptions {

  rrgOptions()
    : degree(0), jointDegree(0)
  {;}
  
  unsigned int degree;
  unsigned int jointDegree;
};

struct swOptions {
  unsigned int neighbours;
  double rewiringProb;
};

struct plodOptions {
  double alpha;
  double beta;
};

struct abOptions 
{
  unsigned int new_edges;
};

struct readFileOptions {
  std::string fileName;
  bool getStates;
};

std::string generateFileName(std::string nameBase, unsigned int id)
{
  std::stringstream s;
  s << nameBase << std::setw(3) << std::setfill('0') << id;
  return s.str();
}

int main(int argc, char* argv[])
{
  Model* model;
  
  std::vector <unsigned int> init;
  
  /******************************************************************/
  // read parameters
  /******************************************************************/
  std::string topology;
  
  unsigned int N = 0;
  double stopTime = 0;
  unsigned int stopRecoveries;
  unsigned int stopInfections;
  unsigned int stopInformations;

  double outputData = 0.;
  int outputGraphviz = 0;
  std::string graphDir = "images";
  std::string outputFileName = "";
  std::ofstream* outputFile = 0;

  unsigned int verbose = 0;
  bool printStats = false;

  unsigned int numSims = 1;
  
  std::string readGraph = ""; // default is to generate graph.
  bool generateIC = true; // default is to generate i.c.

  bool pairs = true;

  po::options_description command_line_options
    ("\nUsage: simulate -p params_file [options]... \n\nMain options");

  command_line_options.add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message including all options")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("print-stats",
     "print stats at the end of run")
    ("params-file,p",po::value<std::string>(),
     "file containing graph parameters")
    ("model-file,m",po::value<std::string>(),
     "file containing model parameters")
    ;

  po::options_description sim_options
    ("\nSimulation options");

  sim_options.add_options()
    ("sim", po::value<std::string>()->default_value("Gillespie"),
     "simulator to use (Gillespie, Chris)")
    ("usemodel", po::value<std::string>()->default_value("InfoSIRS"),
     "model to use (InfoSIRS, ProtectiveSIRS)")
    ("tmax", po::value<double>()->default_value(stopTime),
     "time after which to stop\n(use tmax=0 to just generate graph and exit)")
    ("rmax", po::value<unsigned int>()->default_value(0),
     "number of recoveries after which to stop (if >0)")
    ("imax", po::value<unsigned int>()->default_value(0),
     "number of infections after which to stop (if >0)")
    ("pmax", po::value<unsigned int>()->default_value(0),
     "number of informations after which to stop (if >0)")
    ("output", po::value<double>()->default_value(0.),
     "write output data at arg timesteps")
    ("graphviz,g", po::value<int>()->default_value(0),
     "create graphviz output in the images directory at arg timesteps")
    ("nsims", po::value<unsigned int>()->default_value(1),
     "number of simulation runs to produce (on a given graph)")
    ("graph-dir", po::value<std::string>()->default_value(graphDir),
     "set ouput dir for graphs")
    ("write-file,f", po::value<std::string>(),
     "output data to file (.sim.dat will be appended)")
    ("cluster-coeff",
     "write clustering coefficients to baseName.cluster file")    
    ("write-Js",
     "write adjacency matrices to files  baseName.Jd/i")        
    ;
  
  po::options_description temp_options;
  temp_options.add(command_line_options).add(sim_options);
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(temp_options).
              allow_unregistered().run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  po::notify(vm);

  if (vm.count("verbose")) {
    verbose = 1;
  }
  if (vm.count("very-verbose")) {
    verbose = 2;
  }
  if (vm.count("print-stats")) {
    printStats = true;
  }

  if (vm.count("help")) {
    std::cout << command_line_options << sim_options << std::endl;
    return 0;
  }

  if (vm.count("usemodel")) {
    std::string modelString = vm["usemodel"].as<std::string>();
    if (modelString == "InfoSIRS") {
      model = new Models::InfoSIRS(verbose);
    } else if (modelString == "ProtectiveSIRS") {
      model = new Models::ProtectiveSIRS(verbose);
    } else if (modelString == "Vaccination") {
      model = new Models::VaccinationSIRS(verbose);
    } else {
      std::cerr << "Error: unknown model: " << modelString << std::endl;
      return 1;
    }
  } else {
    std::cerr << "Error: no model specified" << std::endl;
    return 1;
  }
    

  // declare hidden option for suppression of graph output --
  // needed for do_all script so that graphviz output
  // is not generated at each run
  po::options_description hidden_options
    ("\nAdditional options");
  
  hidden_options.add_options()
    ("no-graph",
     "do not produce graphviz output no matter what the other settings")
    ("no-degree-dist",
     "do not write degree distribution to baseName.degree file")
    ("no-pairs",
     "do not consider pairs")
    ;
  
  po::options_description graph_options
    ("\nGraph options");
  
  graph_options.add_options()
    ("vertices,N", po::value<unsigned int>(),
     "number of vertices");
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    graph_options.add_options()
      ((it->getText() + "-topology").c_str(), po::value<std::string>(),
       (it->getText() + "-network topology\n((tri-)lattice,random,random-regular,"+
        "small-world,plod,albert-barabasi,complete,read,null)").c_str());
  }
  graph_options.add_options()
    ("base,b", po::value<std::string>()->default_value
     (model->getVertexStates().begin()->getText()),
     "base state of individuals");
  
  po::options_description assortativity_options
    ("\nAssortativity options");
  
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    for (std::vector<Label>::const_iterator it2 = model->getEdgeTypes().begin();
         it2 != model->getEdgeTypes().end(); it2++) {
      assortativity_options.add_options()
        ((it->getText() + it2->getText() + "-assortativity").c_str(),
         po::value<double>(),
         ("desired assortativity between "+it->getText()+
          "- and "+it2->getText()+"-edges").c_str());
    }
  }
  
  for (std::vector<Label>::const_iterator it =
         model->getVertexStates().begin();
       it != model->getVertexStates().end(); it++) {
    graph_options.add_options()
      (it->getText().c_str(), po::value<unsigned int>()->default_value(0),
       ("number of randomly chosen " + it->getText()).c_str());
  }
  
  // generate topology-specifice options for each type of graph
  
  // lattice
  std::vector<po::options_description*> lattice_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-Lattice options";
    po::options_description* lo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-dim";
    lo->add_options()
      (s.str().c_str(),  po::value<unsigned int>()->default_value(2),
       "number of dimensions");
    s.str("");
    s << it->getText() << "-pb";
    lo->add_options()
      (s.str().c_str(), "periodic boundary conditions");
    lattice_options.push_back(lo);
  }

  // random graph
  std::vector<po::options_description*> rg_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-RandomGraph options";
    po::options_description* ro
      = new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-edges";
    ro->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of edges");
    rg_options.push_back(ro);
  }

  // random regular graph
  std::vector<po::options_description*> rrg_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-RandomRegularGraph options";
    po::options_description* rrgo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-degree";
    rrgo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "degree of the graph G(n,d)");
    s.str("");
    s << it->getText() << "-joint-degree";
    rrgo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "degree of the joint core of two random regular graphs. If not given, the overlapping is random.");
    rrg_options.push_back(rrgo);
  }
  
  // small-world graph
  std::vector<po::options_description*> sw_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-SmallWorld options";
    po::options_description* swo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-neighbours";
    swo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of neighbours of each node");
    s.str("");
    s << it->getText() << "-rewiring-prob";
    swo->add_options()
      (s.str().c_str(), po::value<double>(),
       "rewiring probability");
    sw_options.push_back(swo);
  }

  // plod graph
  std::vector<po::options_description*> plod_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-Power-Law-Out-Degree options";
    po::options_description* plodo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-alpha";
    plodo->add_options()
      (s.str().c_str(), po::value<double>(),
       "alpha (index of power law)");
    s.str("");
    s << it->getText() << "-beta";
    plodo->add_options()
      (s.str().c_str(), po::value<double>(),
       "beta (multiplicative factor of power law)");
    plod_options.push_back(plodo);
  }

  // Albert-Barabasi graph
  std::vector<po::options_description*> ab_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-Albert-Barabasi options";
    po::options_description* abo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-newedges";
    abo->add_options()
      (s.str().c_str(), po::value<unsigned int>()->default_value(1),
       "number of edges to add per vertex");
    ab_options.push_back(abo);
  }

  // read graph
  std::vector<po::options_description*> readFile_options;
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-Read options";
    po::options_description* rfo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << it->getText() << "-file";
    rfo->add_options()
      (s.str().c_str(), po::value<std::string>(),
       "name of graph file to read");
    s.str("");
    s << it->getText() << "-getstates";
    rfo->add_options()
      (s.str().c_str(), 
       "get vertex states as well (and do not randomize initial conditions)");
    readFile_options.push_back(rfo);
  }

  // rewiring options
  po::options_description rewiring_options("Rewiring options");
  
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    std::stringstream s;
    s << it->getText() << "-rewire";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to rewire");
    s.str("");
    s << it->getText() << "-remove";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to remove");
    s.str("");
    s << it->getText() << "-add";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to add");
  }

  po::options_description model_options = model->getOptions();

  // read options from command line
  po::options_description visible_options;
  visible_options.add(command_line_options).add(sim_options).add(graph_options).
    add(rewiring_options).add(assortativity_options).add(model_options);
  
  po::options_description all_options;
  all_options.add(visible_options).add(hidden_options);
  
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    all_options.add(*(lattice_options[it->getId()]));
    all_options.add(*(rg_options[it->getId()]));
    all_options.add(*(rrg_options[it->getId()]));
    all_options.add(*(sw_options[it->getId()]));
    all_options.add(*(plod_options[it->getId()]));
    all_options.add(*(ab_options[it->getId()]));
    all_options.add(*(readFile_options[it->getId()]));
  }

  po::parsed_options parsed = po::command_line_parser(argc, argv).
    options(all_options).allow_unregistered().run();
  try {
    po::store(parsed, vm);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  po::notify(vm);

  std::vector<std::string> unregistered =
    po::collect_unrecognized(parsed.options, po::exclude_positional);
  for (std::vector<std::string>::iterator it = unregistered.begin();
       it != unregistered.end(); it++) {
    std::cerr << "WARNING: ignoring unknown option " << *it << std::endl;
  }
 
  if (vm.count("longhelp")) {
    std::cout << all_options << std::endl;
    return 0;
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

  if (vm.count("vertices")) {
    N = vm["vertices"].as<unsigned int>();
  }

  /******************************************************************/
  // initialize random generator
  /******************************************************************/

  unsigned int seed;
  struct timeval tv;
  gettimeofday(&tv, 0);
  seed = tv.tv_sec + tv.tv_usec;
  boost::mt19937 gen(seed);

  /******************************************************************/
  // create graph variable
  /******************************************************************/
  dualtype_graph graph, saved_graph;
  void (*draw_function)(const dualtype_graph&, const Model&,
                        std::string, double) = 0;
            
  /******************************************************************/
  // read graph from file or generate it
  /******************************************************************/
  
  // adding N vertices to graph
  boost::add_vertices(graph, N);
  
  /******************************************************************/
  // generate edges
  /******************************************************************/

  std::vector<Label> edgeTypes = model->getEdgeTypes();
  
  for (unsigned int i = 0; i < edgeTypes.size(); i++) {

    onetype_graph temp_graph;
    boost::add_vertices(temp_graph, N);
    
    std::stringstream s;
    s << edgeTypes[i].getText() << "-topology";
    if (vm.count(s.str())) {
      topology = vm[s.str()].as<std::string>();
    } else {
      std::cerr << "ERROR: no " << s.str() << " specified" << std::endl;
      std::cerr << command_line_options << graph_options << std::endl;
      return 1;
    }
    
    if (topology == "lattice") {

      draw_function = &draw_lattice;
      
      /******************************************************************/
      // read lattice specific parameters
      /******************************************************************/
      
      latticeOptions opt;
      s.str("");
      s << edgeTypes[i].getText() << "-dim";
      opt.dimensions = vm[s.str()].as<unsigned int>();
      
      opt.sideLength = static_cast<int>(pow(N, 1.0/opt.dimensions));
      if (pow(opt.sideLength, opt.dimensions) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }
      s.str("");
      s << edgeTypes[i].getText() << "-pb";
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
      
      boost::add_edge_structure(temp_graph, li, li_end, Edge(edgeTypes[i].getId()));
      
    } else if (topology == "tri-lattice") {
      
      draw_function = &draw_lattice;
      
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
      
      s << edgeTypes[i].getText() << "-pb";
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
      
      boost::add_edge_structure(temp_graph, tli, tli_end, Edge(edgeTypes[i].getId()));
      
    } else if (topology == "random") {

      /******************************************************************/
      // read random graph specific parameters
      /******************************************************************/
      
      rgOptions opt;
      
      s.str("");
      s << edgeTypes[i].getText() << "-edges";
      if (vm.count(s.str())) {
        opt.edges = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges specified" << std::endl;
        std::cerr << *rg_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate random graph with desired properties
      /******************************************************************/
      typedef boost::erdos_renyi_iterator2<boost::mt19937, dualtype_graph>
        rg_iterator;
      
      //         typedef boost::erdos_renyi_iterator<boost::mt19937, dualtype_graph>
      //           rg_iterator;
      
      double p = (double)opt.edges*2.0/((double)N*(double)(N-1));
      
      // fix for E = p * N * N
      //         p = p * ((double)(N-1) / (double)N);
      
      //         std::cout << "p=" << p << std::endl;
      //         std::cout << p*(double)N*(double)N << std::endl;
      //         std::cout << (int)(p*N*N)/2 << std::endl;
      
      rg_iterator ri(gen, N, p);
      rg_iterator ri_end;
      
      boost::add_edge_structure(temp_graph, ri, ri_end, Edge(edgeTypes[i].getId()));

    } else if (topology == "random-regular") {
      
      /******************************************************************/
      // read random regular graph specific parameters
      /******************************************************************/
      
      rrgOptions opt;
      
      s.str("");
      s << edgeTypes[i].getText() << "-degree";
      if (vm.count(s.str())) {
        opt.degree = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: Graph degree not spcified" << std::endl;
        std::cerr << *rrg_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      s.str("");
      s << edgeTypes[i].getText() << "-joint-degree";
      if (vm.count(s.str())) {
        opt.jointDegree = vm[s.str()].as<unsigned int>();
      } 
      
      /******************************************************************/
      // generate random regular graph with desired properties
      /******************************************************************/
      
      bool success = false;
      unsigned int count = 0;

      typedef std::vector<std::pair<unsigned int, unsigned int> > GraphEdges;
      GraphEdges rrg_edges;
      boost::uniform_01<boost::mt19937, double> uni_gen(gen);

      if (opt.jointDegree > 0) {
        // generate joint graph
        rrg_edges.clear();         
        ++count;
        
        while (!success) {
          success = boost::random_regular_graph(graph, rrg_edges,
                                                opt.jointDegree, N, uni_gen);
          if (success) {            
            for (GraphEdges::iterator it = rrg_edges.begin();
                 it != rrg_edges.end(); ++it){
              boost::add_edge((*it).first, (*it).second,
                              Edge(edgeTypes[i].getId()), temp_graph);
            }
          }
        }
        
        if (verbose) {
          std::cout << "Random regular graph of degree " << opt.jointDegree
                    << " was generated in " << count << " trials\n";
        }

        // copy graph to all edgetypes
        for (unsigned int j = 0; j < edgeTypes.size(); j++) {
          boost::copy_graph(temp_graph, graph, Edge(edgeTypes[j].getId()));
        }
        temp_graph.clear();
      }

      if (opt.degree > 0) {
        // generate additional graph, excluding existing edges
        success = 0;
        count = 0;
        rrg_edges.clear();
        while (!success) {
          success = boost::random_regular_graph(graph, rrg_edges,
                                                opt.degree, N, uni_gen);
          if (success) {            
            for (GraphEdges::iterator it = rrg_edges.begin();
                 it != rrg_edges.end(); ++it){
              boost::add_edge((*it).first, (*it).second,
                              Edge(edgeTypes[i].getId()), temp_graph);
            }
          }
        }
        
        if (verbose) std::cout << "Random regular graph of degree " << opt.degree
                               << " was generated in " << count << " trials\n";
      }
      
    } else if (topology == "small-world") {
      
      draw_function = &draw_ring;

      /******************************************************************/
      // read small-world graph specific parameters
      /******************************************************************/
      
      swOptions opt;
      
      s.str("");
      s << edgeTypes[i].getText() << "-neighbours";
      if (vm.count(s.str())) {
        opt.neighbours = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of neighbours specified" << std::endl;
        std::cerr << *sw_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      s.str("");
      s << edgeTypes[i].getText() << "-rewiring-prob";
      if (vm.count(s.str())) {
        opt.rewiringProb = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no rewiring probability" << std::endl;
        std::cerr << *sw_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate small-world graph with desired properties
      /******************************************************************/
      typedef boost::small_world_iterator<boost::mt19937, dualtype_graph>
        sw_iterator;
      
      sw_iterator swi(gen,N,opt.neighbours, opt.rewiringProb);
      sw_iterator swi_end;
      
      boost::add_edge_structure(temp_graph, swi, swi_end, Edge(edgeTypes[i].getId()));
      
    } else if (topology == "plod") {
      
      /******************************************************************/
      // read plod graph specific parameters
      /******************************************************************/
      
      plodOptions opt;
      
      s.str("");
      s << edgeTypes[i].getText() << "-alpha";
      if (vm.count(s.str())) {
        opt.alpha = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no alpha specified" << std::endl;
        std::cerr << *plod_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      s.str("");
      s << edgeTypes[i].getText() << "-beta";
      if (vm.count(s.str())) {
        opt.beta = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no beta specified" << std::endl;
        std::cerr << *plod_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate plod graph with desired properties
      /******************************************************************/
      typedef boost::plod_iterator<boost::mt19937, dualtype_graph>
        plod_iterator;
      
      plod_iterator plodi(gen,N,opt.alpha,opt.beta);
      plod_iterator plodi_end;
      
      boost::add_edge_structure(temp_graph, plodi, plodi_end, Edge(edgeTypes[i].getId()));

    } else if (topology == "albert-barabasi") {
      
      /******************************************************************/
      // read Albert-Barabasi graph specific parameters
      /******************************************************************/
      
      abOptions opt;
      
      s.str("");
      s << edgeTypes[i].getText() << "-newedges";
      if (vm.count(s.str())) {
        opt.new_edges = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges to add per vertex specified" << std::endl;
        std::cerr << *ab_options[edgeTypes[i].getId()] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate Albert-Barabasi with desired properties
      /******************************************************************/
      typedef boost::albert_barabasi_iterator<boost::mt19937, dualtype_graph>
        albert_barabasi_iterator;
      
      albert_barabasi_iterator ab(gen,N,opt.new_edges);
      albert_barabasi_iterator ab_end;
      
      boost::add_edge_structure(temp_graph, ab, ab_end, Edge(edgeTypes[i].getId()));

    } else if (topology == "complete") {

      boost::graph_traits<dualtype_graph>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; vi++) {
        boost::graph_traits<dualtype_graph>::vertex_iterator vi2;
        for (vi2 = vi+1; vi2 != vi_end; vi2++) {
          add_edge(*vi, *vi2, Edge(edgeTypes[i].getId()), temp_graph);
        }
      }
      
    } else if (topology == "copy") {
      
      /******************************************************************/
      // test if we have graph to copy from
      /******************************************************************/

      if (num_edges(graph) > 0) {
        boost::copy_graph(graph, temp_graph, Edge(edgeTypes[i].getId()));
      } else {
        // push back for later
        edgeTypes.push_back(edgeTypes[i]);
        if (edgeTypes.size() == 2*model->getEdgeTypes().size()) {
          std::cerr << "ERROR: no edges to copy from, " << edgeTypes[i]
                    << "-graph to remain empty" << std::endl;
          std::cerr << command_line_options << graph_options << std::endl;
          return 1;
        }
      }

    } else if (topology == "read") {
        
      /******************************************************************/
      // read file graph specific parameters
      /******************************************************************/
      readFileOptions opt;

      s.str("");
      s << edgeTypes[i].getText() << "-file";
      if (vm.count(s.str())) {
        if (readGraph.size() == 0) {
          readGraph = vm[s.str()].as<std::string>();
        }
        opt.fileName = vm[s.str()].as<std::string>();
      } else {
        if (readGraph.size() == 0) {
          std::cerr << "ERROR: no file name specified" << std::endl;
          std::cerr << *readFile_options[edgeTypes[i].getId()] << std::endl;
          return 1;
        } else {
          opt.fileName = readGraph;
        }
      }

      s.str("");
      s << edgeTypes[i].getText() << "-getstates";
      if (vm.count(s.str()) && generateIC == true) {
        opt.getStates = true;
        generateIC = false;
      } else {
        opt.getStates = false;
      }
      
      // reading graph structure and initial state from file
      if (read_graph(graph, *model, opt.fileName, edgeTypes[i].getId(),
                     opt.getStates, verbose) == 0) {
        
        // update number of vertices
        N = num_vertices(graph);
        if (verbose) {
          std::cout << "graph file " << readGraph << " was read ok\n";
        }
      } else {
        std::cerr << "ERROR: something wrong in read graph from "
                  << readGraph << std::endl;
        return 1;
      }
    } else if (topology != "null") {
      std::cerr << "ERROR: unknown " << s.str() << ": " << topology
                << std::endl;
      std::cerr << command_line_options << graph_options << std::endl;
    }

    // do graph modifications if desired
    if (vm.count(edgeTypes[i].getText()+"-remove")) {
      double removeFraction = vm[edgeTypes[i].getText()+"-remove"].as<double>();
      if (removeFraction > 0) {
        boost::removeEdges(temp_graph, gen, removeFraction);
      }
    }
    if (vm.count(edgeTypes[i].getText()+"-rewire")) {
      double rewireFraction = vm[edgeTypes[i].getText()+"-rewire"].as<double>();
      if (rewireFraction > 0) {
        int rewire_result = -1;
        unsigned int count = 0;
        onetype_graph rewire_graph;
        while (rewire_result < 0) {
          rewire_graph = temp_graph;
          rewire_result = boost::rewireEdges(rewire_graph, gen,
                                             rewireFraction, verbose);
          ++count;
        }
        temp_graph = rewire_graph;
        if (verbose) {
          std::cout << "graph rewired after " << count << " trials\n";
        }
      }
    }
    if (vm.count(edgeTypes[i].getText()+"-add")) {
      double addFraction = vm[edgeTypes[i].getText()+"-add"].as<double>();
      boost::addEdges(temp_graph, gen, addFraction);
    }
    
    // copy edges to main graph
    boost::graph_traits<dualtype_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(temp_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
               temp_graph[*ei], graph);
    }         
  }
  
  // checking graph  
  if (num_vertices(graph) == 0) {
    std::cerr << "ERROR: no vertices" << std::endl;
    std::cerr << command_line_options << graph_options << std::endl;
    return 1;
  }

  if (draw_function == 0) draw_function = &write_graph;

  /******************************************************************/
  // create given assortativity if desired
  /******************************************************************/
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    for (std::vector<Label>::const_iterator it2 = model->getEdgeTypes().begin();
         it2 != model->getEdgeTypes().end(); it2++) {
      std::string s = it->getText()+it2->getText()+"-assortativity";
      if (vm.count(s.c_str())) {
        double ass = (vm[s.c_str()].as<double>());
        rewire_assortatively(graph, gen, ass,
                             Edge(it->getId()), Edge(it2->getId()), verbose);
      }
    }
  }
  
  /******************************************************************/
  // mark parallel edges
  /******************************************************************/

  if (pairs) {
    // mark parallel edges
    unsigned int parallel_edges = mark_parallel_edges(graph);
    if (verbose) std::cout << "No. of parallel edges is: " << parallel_edges
                           << std::endl;
  }

  // set base file name
  std::string baseFileName;
  if (vm.count("write-file")) {
    baseFileName = vm["write-file"].as<std::string>();
    outputFileName = baseFileName+".sim.dat";
    
    try {
      outputFile = new std::ofstream();
      outputFile->open(outputFileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
    }
    
    if (outputGraphviz >=0) {
      
      std::string outputGraphName =
        baseFileName+".graph";
      
      // write graph
      write_graph(graph, *model, outputGraphName, -1);
    }
    
    // calculate degree distribution
    std::string degreeFileName = baseFileName+".degree";
    bool status = write_degree(graph, *model, degreeFileName);
    if (verbose)
      if (!status) {
        std::cout << "degree file " << degreeFileName << " was written ok\n";
      } else {
        std::cout << "ERROR: something wrong in writing degree file "
                  << degreeFileName << std::endl;
      }

    // create sparse adjacency matrices and clustering coefficients
    if (vm.count("cluster-coeff")) {
      bool writeJs = false;
      if (vm.count("write-Js")) writeJs = true;
      
      bool status = boost::cluster_coeff(graph, baseFileName, writeJs);
      if (verbose)
        if (!status) {
          std::cout << "cluster coeff file was written ok\n";
        } else {
          std::cout << "ERROR: something wrong in writing cluster coeff file "
                    << std::endl;
        }
    }
  }

  numSims = vm["nsims"].as<unsigned int>();
  
  unsigned int extLength = 0;
  
  if (numSims > 1) {
    std::stringstream ext;
    ext << numSims;
    extLength = ext.str().length();
    std::cout << "Running " << numSims << " simulations" << std::endl;
  }

  for (unsigned int nSim = 1; nSim <= numSims; nSim++) {

    if (printStats) {
      std::cout << "----- " << "run #" << nSim << std::endl;
    } else if (numSims > 0) {
      std::cout << ".";
      std::cout.flush();
    }
    
    std::stringstream fileName;
    fileName << baseFileName << std::setfill('0') << std::setw(extLength)
             << nSim;

    stopTime = vm["tmax"].as<double>();

    /******************************************************************/
    // generate new initial state
    /******************************************************************/
   
    if (generateIC) { // generate new initial state
      
      /******************************************************************/
      // set initial vertex states
      /******************************************************************/
      
      // how many random vertices of each state are to be initialized
      // over the background of the base state
      for (unsigned int i = 0; i < model->getVertexStates().size(); i++) {
        std::stringstream ss;
        ss << model->getVertexStates()[i].getText();
        std::string s(ss.str());
        unsigned int random = 0;
        if (vm.count(s.c_str())) {
          random = vm[s.c_str()].as<unsigned int>();
        }
        init.push_back(random);
      }
      
      
      std::string baseString = vm["base"].as<std::string>();
      unsigned int baseState = 0;
      
      // assign baseState
      while (baseState < model->getVertexStates().size() &&
             (model->getVertexStates()[baseState].getText() != baseString)) {
        baseState++;
      }
      
      // set random vertices of baseState to zero
      if (baseState < model->getVertexStates().size()) {
        init[baseState] = 0;
      } else {
        std::cerr << "ERROR: no unknown base state: " << baseString << std::endl;
        std::cerr << command_line_options << sim_options << std::endl;
        return 1;
      }              
      
      /******************************************************************/
      // generate vertices' state
      /******************************************************************/
      
      // add N vertices in state baseState to graph
      boost::graph_traits<dualtype_graph>::vertex_iterator vi, vi_end;      
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi)
        graph[*vi].state = baseState;
      
      // sum over vector init to make sure the sum is less than N
      unsigned int initSum = 0;
      for (std::vector<unsigned int>::iterator it = init.begin();
           it != init.end(); it++) {
        initSum += (*it);
      }   
      if (initSum > N) {
        std::cerr << "Error: number of vertices to select randomly"
                  << " higher than number of total vertices" << std::endl;
      }
      
      // inserting init[i] vertices of type i
      boost::graph_traits<dualtype_graph>::vertex_descriptor v;
      for (unsigned int i=0; i<model->getVertexStates().size(); i++) {
        for (unsigned int j=0; j<init[i]; j++) {
          bool inserted = false;
          while (!inserted) {
            v = boost::random_vertex(graph, gen);
            if (graph[v].state == baseState) {
              graph[v].state = i;
              inserted = true;
            }
          }
          if (verbose >= 2) {
            std::cout << "Vertex #" << v << " is assigned state " 
                      << model->getVertexStates()[i] << std::endl;
          }
        }
      }
      
    } else { // if generateIC is not set
      if (nSim == 1) {
        // save graph states
        saved_graph = graph;
      } else {
        // recover graph states
        graph = saved_graph;
      }
    }
    
    /******************************************************************/
    // initialize model
    /******************************************************************/
    
    model->Init(vm);
    if (verbose >=1) model->Print();
    
    /******************************************************************/
    // create simulator
    /******************************************************************/
    
    Simulator* sim;
    
    if (vm.count("sim")) {
      std::string simType = vm["sim"].as<std::string>();
      if (simType == "Gillespie") {
        sim = new Simulators::GillespieSimulator<boost::mt19937, dualtype_graph>
          (gen, graph, *model, verbose);
      } else if (simType == "Chris") {
        sim = new Simulators::ChrisSimulator<boost::mt19937, dualtype_graph>
          (gen, graph, *model, verbose);
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
    
    stopRecoveries = vm["rmax"].as<unsigned int>();
    stopInfections = vm["imax"].as<unsigned int>();
    stopInformations = vm["pmax"].as<unsigned int>();
    outputData = vm["output"].as<double>();
        
    // graph directory
    if (vm.count("graph-dir")) {
      graphDir = vm["graph-dir"].as<std::string>();
    }
    std::stringstream graphDirName(graphDir, std::ios::in | std::ios::out |
                                   std::ios::ate);
    graphDirName << "/run" << std::setfill('0') << std::setw(extLength) << nSim;
    
    // timesteps after which to write graphviz output
    if (vm.count("graphviz")) {
      outputGraphviz = vm["graphviz"].as<int>();
    }
    // no-graph overrides other graph options
    if (vm.count("no-graph")) {
      outputGraphviz = -1;
    }
    // do not consider pairs
    if (vm.count("no-pairs")) {
      pairs = false;
    }
    
    /******************************************************************/
    // open output file
    /******************************************************************/
    if (vm.count("write-file")) {
      outputFileName = fileName.str()+".sim.dat";
      
      try {
        outputFile = new std::ofstream();
        outputFile->open(outputFileName.c_str(), std::ios::out);
      }
      catch (std::exception &e) {
        std::cerr << "Unable to open output file: " << e.what() << std::endl;
      }
    }
    
    /******************************************************************/
    // initialize Simulator
    /******************************************************************/
    sim->initialize();
    
    // print time
    if (verbose)
      std::cout << "time elapsed: " << sim->getTime() << std::endl;
    
    // GraphViz output
    if (outputGraphviz >= 0) {
      // create graph directories
      if (nSim == 1) mkdir(graphDir.c_str(), 0755);
      mkdir(graphDirName.str().c_str(), 0755);
      draw_function(graph, *model, (graphDirName.str() + "/frame000"), -1);
    }
    
    // prints data to outputFile
    std::string lastLine = "";
    if (outputFile) {
      lastLine = write_graph_data(graph, *model, sim->getTime(), *outputFile,
                                  pairs);
    }
    if (verbose) print_graph_statistics(graph, *model, pairs);
    
    /******************************************************************/
    // run simulation
    /******************************************************************/
    double nextDataStep = outputData;
    
    unsigned int steps = 0;
    unsigned int outputNum = 1;
    
    bool doSim = !(stopTime == 0 && stopInfections == 0 &&
                   stopRecoveries == 0 && stopInformations == 0);
    
    while ((stopTime == 0 || sim->getTime()<stopTime) &&
           (stopInfections == 0 || (sim->getNumInfections() < stopInfections &&
                                    sim->getNumInfections()+1 > sim->getNumRecoveries())) &&
           (stopInformations == 0 || (sim->getNumInformations() < stopInformations && 
                                      sim->getNumInformations()+1 > sim->getNumForgettings())) &&
           doSim && sim->updateState()) {
      
      if (verbose >= 2) {
        print_graph_statistics(graph, *model, pairs);
      }
      
      if (verbose && steps%100 == 0) {
        std::cout << "time elapsed: " << sim->getTime() << std::endl;
      }
      
      if ((outputGraphviz > 0) && (steps % outputGraphviz == 0)) {
        draw_function(graph, *model,
                      generateFileName((graphDirName.str() +"/frame"),
                                       outputNum),
                      sim->getTime());
        ++outputNum;
      }
      if (outputFile && sim->getTime() > nextDataStep) {
        lastLine = 
          write_graph_data(graph, *model, sim->getTime(), *outputFile, pairs);
        if (outputData > 0) {
          do {
            nextDataStep += outputData;
          } while (sim->getTime() > nextDataStep);
        }
      }
      ++steps;
    }
    
    if (stopTime == 0) stopTime = sim->getTime();
    
    if (verbose) std::cout << "Final status (" << sim->getTime() << "): " 
                           << std::endl;
    if (doSim && outputGraphviz >= 0) {
      draw_function(graph, *model,
                    generateFileName((graphDirName.str()+"/frame"),
                                     outputNum),sim->getTime());
    }
    
    if (outputFile) {
      if (sim->getTime() < stopTime) {
        lastLine = 
          write_graph_data(graph, *model, sim->getTime(), *outputFile, pairs);
      }
      *outputFile << stopTime << '\t' << lastLine;
      outputFile->close();
      delete outputFile;
//       std::ofstream statsFile;
//       std::string statsFileName = fileName.str()+".stats";
//       statsFile.open(statsFileName.c_str(), std::ios::out);
//       statsFile << "Cumulative number of infections: "
//                 << sim->getNumInfections() << std::endl;
//       statsFile.close();
    }
    
    if (verbose || printStats) {
      std::cout << "Cumulative number of infections: " << sim->getNumInfections() 
                << std::endl;
      std::cout << "Cumulative number of informations: " << sim->getNumInformations() 
                << std::endl;
    }
    
    if (verbose) print_graph_statistics(graph, *model, pairs);
    
    // free memory
    delete sim;
    
  }

  if (vm.count("write-file")) {
    
    std::ofstream gpFile;
    std::string gpFileName = baseFileName+".gp";
    
    try {
      gpFile.open(gpFileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "... unable to open gnuplot output file " 
                << gpFileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return 1; 
    }
    
    gpFile << "### model parameters generated by simulate" << std::endl;
    gpFile << "N=" << num_vertices(graph) << std::endl;
    gpFile << "Tmax=" << stopTime << std::endl;
    gpFile << "### end of model parameters" << std::endl;
    
    try {
      gpFile.close();
    }
    catch (std::exception &e) {
      std::cerr << "... unable to close gnuplot output file " 
                << gpFileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return 1; 
    }
  }
  
  // free memory
  for (std::vector<Label>::const_iterator it = model->getEdgeTypes().begin();
       it != model->getEdgeTypes().end(); it++) {
    delete lattice_options[it->getId()];
    delete rg_options[it->getId()];
    delete rrg_options[it->getId()];
    delete sw_options[it->getId()];
    delete plod_options[it->getId()];
    delete ab_options[it->getId()];
    delete readFile_options[it->getId()];
  }
  
  delete model;
  
  std::cout << std::endl;
  
  return 0;
  
}
