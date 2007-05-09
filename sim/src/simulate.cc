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
// #include <boost/graph/erdos_renyi_generator.hpp>
#include <math.h>

#include "lattice_generator.hh"
#include "graph_structure.hh"
#include "erdos_renyi_generator2.hh"
#include "albert_barabasi_generator.hh"
#include "visualize_graph.hh"
#include "cluster_coeffs.hh"

#include "GillespieSimulator.hh"
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

struct rrgOptions {
  unsigned int d;
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
  InfoSIRS model;
  
  std::vector <unsigned int> init;
  
  /******************************************************************/
  // read parameters
  /******************************************************************/
  std::string topology;
  
  unsigned int N = 0;
  double stop;

  double outputData = 0.;
  int outputGraphviz = 0;
  std::string graphDir = "images";
  std::string outputFileName = "";
  std::ofstream* outputFile = 0;

  bool verbose = false;
  
  std::string readGraph = ""; // default is to generate graph.
  bool generateIC = true; // default is to generate i.c.
            
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
     "time after which to stop\n(use tmax=0 to just generate graph and exit)")
    ("output", po::value<double>()->default_value(0.),
     "write output data at arg timesteps")
    ("graphviz,g", po::value<int>()->default_value(0),
     "create graphviz output in the images directory at arg timesteps")
    ("graph-dir", po::value<std::string>()->default_value(graphDir),
     "set ouput dir for graphs")
    ("write-file,f", po::value<std::string>(),
     "output data to file (.sim.dat will be appended)")
    ("generate-ode-ic-file", po::value<std::string>(),
     "generate init file for ode solver and stop")
    ("cluster-coeff",
     "write clustering coefficients to baseName.cluster file")    
    ("write-Js",
     "write adjacency matrices to files  baseName.Jd/i")        
    ;
  
  // declare hidden option for suppression of graph output --
  // needed for do_all script so that graphviz output
  // is not generated at each run
  po::options_description hidden_option;
  hidden_option.add_options()
    ("no-graph",
     "do not produce graphviz output no matter what the other settings")
    ("no-degree-dist",
     "do not write degree distribution to baseName.degree file")
    ;
  
  po::options_description graph_options;
  graph_options.add_options()
    ("vertices,N", po::value<unsigned int>(),
     "number of vertices")
    ("d-topology", po::value<std::string>(),
     "disease network topology\n((tri-)lattice,random,random-regular,small-world,plod,albert-barabasi,complete,read,null)")
    ("i-topology", po::value<std::string>(),
     "information network topology\n((tri-)lattice,random,random-regular,small-world,plod,albert-barabasi,complete,read,null)")
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
    s << model.getEdgeTypes()[i].getText() << "-Lattice Options";
    po::options_description* lo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-dim";
    lo->add_options()
      (s.str().c_str(),  po::value<unsigned int>()->default_value(2),
       "number of dimensions");
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-pb";
    lo->add_options()
      (s.str().c_str(), "periodic boundary conditions");
    lattice_options.push_back(lo);
  }

  // random graph
  std::vector<po::options_description*> rg_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-RandomGraph Options";
    po::options_description* ro
      = new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-edges";
    ro->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of edges");
    rg_options.push_back(ro);
  }

  // random regular graph
  std::vector<po::options_description*> rrg_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-RandomRegularGraph Options";
    po::options_description* rrgo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-degree";
    rrgo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "degree of the graph G(n,d)");
    rrg_options.push_back(rrgo);
  }
  
  // small-world graph
  std::vector<po::options_description*> sw_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-SmallWorld Options";
    po::options_description* swo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-neighbours";
    swo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of neighbours of each node");
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-rewiring-prob";
    swo->add_options()
      (s.str().c_str(), po::value<double>(),
       "rewiring probability");
    sw_options.push_back(swo);
  }

  // plod graph
  std::vector<po::options_description*> plod_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-Power-Law-Out-Degree Options";
    po::options_description* plodo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-alpha";
    plodo->add_options()
      (s.str().c_str(), po::value<double>(),
       "alpha (index of power law)");
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-beta";
    plodo->add_options()
      (s.str().c_str(), po::value<double>(),
       "beta (multiplicative factor of power law)");
    plod_options.push_back(plodo);
  }

  // Albert-Barabasi graph
  std::vector<po::options_description*> ab_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-Albert-Barabasi Options";
    po::options_description* abo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-newedges";
    abo->add_options()
      (s.str().c_str(), po::value<unsigned int>()->default_value(1),
       "number of edges to add per vertex");
    ab_options.push_back(abo);
  }

  std::vector<po::options_description*> readFile_options;
  for (unsigned int i = 0; i < model.getEdgeTypes().size(); i++) {
    std::stringstream s;
    s << model.getEdgeTypes()[i].getText() << "-Read Options";
    po::options_description* rfo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-file";
    rfo->add_options()
      (s.str().c_str(), po::value<std::string>(),
       "name of graph file to read");
    s.str("");
    s << model.getEdgeTypes()[i].getText() << "-getstates";
    rfo->add_options()
      (s.str().c_str(), 
       "get vertex states as well (and do not randomize initial conditions)");
    readFile_options.push_back(rfo);
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
    all_options.add(*(rrg_options[i]));
    all_options.add(*(sw_options[i]));
    all_options.add(*(plod_options[i]));
    all_options.add(*(ab_options[i]));
    all_options.add(*(readFile_options[i]));
  }
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all_options).
            allow_unregistered().run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << command_line_options << main_options << std::endl;
    return 0;
  }

  if (vm.count("longhelp")) {
    std::cout << all_options << std::endl;
    return 0;
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
  
  if (vm.count("vertices")) {
    N = vm["vertices"].as<unsigned int>();
  }

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
  boost::uniform_01<boost::mt19937, double> uni_gen(gen); // for rrg 
  dualtype_graph graph;
  Simulator* sim;

  if (vm.count("sim")) {
    std::string simType = vm["sim"].as<std::string>();
    if (simType == "Gillespie") {
      sim = new GillespieSimulator<boost::mt19937, dualtype_graph>
        (gen, graph, model);
//     } else if (simType == "Chris") {
//       sim = new ChrisSimulator<boost::mt19937, dualtype_graph>
//         (gen, graph, model);
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
  
  /******************************************************************/
  // read graph from file or generate it
  /******************************************************************/
  
  // adding N vertices to graph
  boost::add_vertices(graph, N);
  
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
      s << model.getEdgeTypes()[i].getText() << "-dim";
      opt.dimensions = vm[s.str()].as<unsigned int>();
      
      opt.sideLength = static_cast<int>(pow(N, 1.0/opt.dimensions));
      if (pow(opt.sideLength, opt.dimensions) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-pb";
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
      
      s << model.getEdgeTypes()[i].getText() << "-pb";
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
      s << model.getEdgeTypes()[i].getText() << "-edges";
      if (vm.count(s.str())) {
        opt.edges = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << *rg_options[i] << std::endl;
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
      
      boost::add_edge_structure(temp_graph, ri, ri_end, Edge(i));

    } else if (topology == "random-regular") {
      
      /******************************************************************/
      // read random regular graph specific parameters
      /******************************************************************/
      
      rrgOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-degree";
      if (vm.count(s.str())) {
        opt.d = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: Graph degree not spcified" << std::endl;
        std::cerr << std::endl;
        std::cerr << rrg_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate random regular graph with desired properties
      /******************************************************************/
      
      bool success = 0;
      unsigned int count = 0;
      typedef std::vector<std::pair<unsigned int, unsigned int> > GraphEdges;
      GraphEdges rrg_edges;
      
      while (!success) {
        success = boost::random_regular_graph(rrg_edges, opt.d, N, uni_gen);
        
        if (success) {            
          for (GraphEdges::iterator it = rrg_edges.begin(); it != rrg_edges.end(); ++it) 
            boost::add_edge((*it).first, (*it).second, Edge(i), temp_graph);            
        } 

        rrg_edges.clear();         
        ++count;
      }
      
      if (verbose) std::cout << "Random regular graph of degree " << opt.d
                             << " was generated in " << count << " trials\n";
      
    } else if (topology == "small-world") {
      
      /******************************************************************/
      // read small-world graph specific parameters
      /******************************************************************/
      
      swOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-neighbours";
      if (vm.count(s.str())) {
        opt.neighbours = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of neighbours specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << *sw_options[i] << std::endl;
        return 1;
      }
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-rewiring-prob";
      if (vm.count(s.str())) {
        opt.rewiringProb = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no rewiring probability" << std::endl;
        std::cerr << std::endl;
        std::cerr << *sw_options[i] << std::endl;
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
      
    } else if (topology == "plod") {
      
      /******************************************************************/
      // read plod graph specific parameters
      /******************************************************************/
      
      plodOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-alpha";
      if (vm.count(s.str())) {
        opt.alpha = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no alpha specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << *plod_options[i] << std::endl;
        return 1;
      }
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-beta";
      if (vm.count(s.str())) {
        opt.beta = vm[s.str()].as<double>();
      } else {
        std::cerr << "ERROR: no beta specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << *plod_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate plod graph with desired properties
      /******************************************************************/
      typedef boost::plod_iterator<boost::mt19937, dualtype_graph>
        plod_iterator;
      
      plod_iterator plodi(gen,N,opt.alpha,opt.beta);
      plod_iterator plodi_end;
      
      boost::add_edge_structure(temp_graph, plodi, plodi_end, Edge(i));

    } else if (topology == "albert-barabasi") {
      
      /******************************************************************/
      // read Albert-Barabasi graph specific parameters
      /******************************************************************/
      
      abOptions opt;
      
      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-newedges";
      if (vm.count(s.str())) {
        opt.new_edges = vm[s.str()].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges to add per vertex specified" << std::endl;
        std::cerr << std::endl;
        std::cerr << *ab_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate Albert-Barabasi with desired properties
      /******************************************************************/
      typedef boost::albert_barabasi_iterator<boost::mt19937, dualtype_graph>
        albert_barabasi_iterator;
      
      albert_barabasi_iterator ab(gen,N,opt.new_edges);
      albert_barabasi_iterator ab_end;
      
      boost::add_edge_structure(temp_graph, ab, ab_end, Edge(i));

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
      
    } else if (topology == "read") {

      readFileOptions opt;

      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-file";
      if (vm.count(s.str())) {
        if (readGraph.size() == 0) {
          readGraph = vm[s.str()].as<std::string>();
        }
        opt.fileName = vm[s.str()].as<std::string>();
      } else {
        if (readGraph.size() == 0) {
          std::cerr << "ERROR: no file name specified" << std::endl;
          std::cerr << std::endl;
          std::cerr << *readFile_options[i] << std::endl;
          return 1;
        } else {
          opt.fileName = readGraph;
        }
      }

      s.str("");
      s << model.getEdgeTypes()[i].getText() << "-getstates";
      if (vm.count(s.str()) && generateIC == true) {
        opt.getStates = true;
        generateIC = false;
      } else {
        opt.getStates = false;
      }
      
      // reading graph structure and initial state from file
      if (read_graph(graph, model, opt.fileName, i, opt.getStates, verbose) == 0) {
        
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
      std::cerr << std::endl;
      std::cerr << command_line_options << graph_options << std::endl;
      return 1;
    }
    
    // copy edges to main graph
    boost::graph_traits<dualtype_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(temp_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
               Edge(temp_graph[*ei].type), graph);
    }         
  }

  if (num_vertices(graph) == 0) {
    std::cerr << "ERROR: no vertices" << std::endl;
    std::cerr << std::endl;
    std::cerr << command_line_options << graph_options << std::endl;
    return 1;
  }
    
  // mark parallel edges
  unsigned int parallel_edges = mark_parallel_edges(graph);
  if (verbose) std::cout << "No. of parallel edges is: " << parallel_edges
                         << std::endl;

  // create sparse adjacency matrices and clustering coefficients
  if (vm.count("cluster-coeff")) {
    std::string baseFileName = (vm["write-file"].as<std::string>());
    
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

  
  /******************************************************************/
  // generate new initial state
  /******************************************************************/
   
  if (generateIC) { // generate new initial state
      
    /******************************************************************/
    // set initial vertex states
    /******************************************************************/
      
    // how many random vertices of each state are to be initialized
    // over the background of the base state
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
      
    // assign baseState
    while (baseState < model.getVertexStates().size() &&
           (model.getVertexStates()[baseState].getText() != baseString)) {
      baseState++;
    }
      
    // set random vertices of baseState to zero
    if (baseState < model.getVertexStates().size()) {
      init[baseState] = 0;
    } else {
      std::cerr << "ERROR: no unknown base state: " << baseString << std::endl;
      std::cerr << std::endl;
      std::cerr << command_line_options << main_options << std::endl;
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
      
  } // end generateIC

  /******************************************************************/
  // check for generate-ode-ic-file
  /******************************************************************/

  if (vm.count("generate-ode-ic-file")) {
    std::ofstream odeIcFile;
    std::string odeIcFileName =
      (vm["generate-ode-ic-file"].as<std::string>());

    // open ode ic file for writing
    try {
      odeIcFile.open(odeIcFileName.c_str(), std::ios::out);
    }
    catch (std::exception& e) {
      std::cerr << "ERROR:  unable to open file for ode ic file: "
                << e.what() << std::endl;
      return 1;
    }
      
    // write ode ic file
    bool status = write_ode_ic_file(graph, model, odeIcFile);
      
    // print message
    if (verbose) 
      if (!status) {
        std::cout << "ode ic file " << odeIcFileName << " written ok\n";
      } else {
        std::cout << "ERROR: something wrong in writing ode ic file "
                  << odeIcFileName << std::endl;
        return 1;
      }

    // close file
    odeIcFile.close();
      
    // exit
    return 0;   
  }
   
  /******************************************************************/
  // open output file
  /******************************************************************/
  if (vm.count("write-file")) {
    outputFileName = (vm["write-file"].as<std::string>())+".sim.dat";

    try {
      outputFile = new std::ofstream();
      outputFile->open(outputFileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
    }

    if (outputGraphviz >=0) {

      std::string outputGraphName =
        (vm["write-file"].as<std::string>())+".graph";

      // write graph
      write_graph(graph, model, outputGraphName, -1);
    }

    // calculate degree distribution
    if (!vm.count("no-degree-dist")) {      
      std::string degreeFileName = (vm["write-file"].as<std::string>())+".degree";
      bool status = write_degree(graph, model, degreeFileName);
      if (verbose)
        if (!status) {
          std::cout << "degree file " << degreeFileName << " was written ok\n";
        } else {
          std::cout << "ERROR: something wrong in writing degree file "
                    << degreeFileName << std::endl;
        }
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
  if (outputGraphviz >= 0)
    write_graph(graph, model, (graphDir + "/frame000"), -1);

  // prints data to outputFile
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
  unsigned int outputNum = 1;

  while (sim->updateState() && sim->getTime()<stop) {
    
    if (verbose && steps%100 == 0) {
      std::cout << "time elapsed: " << sim->getTime() << std::endl;
    }

    if ((outputGraphviz > 0) && (steps % outputGraphviz == 0)) {
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
  if (stop > 0 && outputGraphviz >= 0) {
      write_graph(graph, model,
                  generateFileName((graphDir +"/frame"),outputNum),
                  sim->getTime());
  }
  
  if (outputFile) {
    *outputFile << stop << '\t' << lastLine;
    outputFile->close();
    delete outputFile;
  }
  
  if (verbose) print_graph_statistics(graph, model);

  // free memory
  delete sim;
  
  return 0;
}
