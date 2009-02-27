/*! \file network.cc
  \brief The main graph program.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <algorithm>

#include <sys/time.h>
#include <sys/stat.h>

#include <boost/graph/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <math.h>

#include "graph_io.hh"
#include "graph_structure.hh"
#include "graph_statistics.hh"
#include "lattice_generator.hh"
#include "tree_generator.hh"
#include "erdos_renyi_generator2.hh"
#include "albert_barabasi_generator.hh"
#include "community_generator.hh"
#include "cluster_coeffs.hh"
#include "assortativity.hh"
#include "degree_overlap.hh"
#include "path_length.hh"
#include "community_structure.hh"

namespace po = boost::program_options;

//! Classes and functions of the boost libraries.
namespace boost {
//! Special classes and functions of the boost libraries.
  namespace detail {}
}

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::no_property, Edge> multitype_graph;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                              boost::no_property, Edge> onetype_graph;

struct latticeOptions {
  latticeOptions() : sideLength(0) {}
  unsigned int sideLength;
  unsigned int dimensions;
  bool periodicBoundary;
};

struct treeOptions {
  unsigned int branches;
};

struct rgOptions {
  unsigned int edges;
};

struct rrgOptions {
  rrgOptions() : degree(0), jointDegree(0) {}
  unsigned int degree;
  unsigned int jointDegree;
};

struct swOptions {
  unsigned int degree;
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
};

struct communityOptions {
  double delta, pl, pr, pd;
  unsigned int iterations;
};

int main(int argc, char* argv[])
{
  /******************************************************************/
  // read parameters
  /******************************************************************/
  std::string topology;

  unsigned int nEdgeTypes = 2;
  unsigned int N = 0;
  unsigned int verbose = 0;

  std::string readGraph = ""; // default is to generate graph.
  bool readAll = false;
  std::string edgeLabels = "di";

  bool allStats = false;

  po::options_description main_options
    ("Usage: network -p params_file [options]... \n\nMain options");

  main_options.add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message including all options")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ;

  po::options_description edgetype_options
    ("Edge types");

  edgetype_options.add_options()
    ("ntypes,n",po::value<unsigned int>()->default_value(nEdgeTypes),
     "number of edge types")
    ("labels,l",po::value<std::string>()->default_value(edgeLabels),
     "labels of the edge types")
    ;

  po::options_description output_options
    ("Output options");
  
  output_options.add_options()
    ("output-file,o", po::value<std::string>(),
     "output graph to file (.graph will be appended)")
    ("split,s", 
     "split graph output input")
    ;

  po::options_description temp_options;
  temp_options.add(main_options).add(edgetype_options).add(output_options);
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

  if (vm.count("help")) {
    std::cout << temp_options << std::endl;
    return 0;
  }

  if (vm.count("labels")) {
    edgeLabels = vm["labels"].as<std::string>();
  }

  if (vm.count("ntypes")) {
    nEdgeTypes = vm["ntypes"].as<unsigned int>();
  }
  
  if (edgeLabels.size() < nEdgeTypes) {
    edgeLabels.clear();
    for (unsigned int i = 0; i < nEdgeTypes; ++i) {
      std::stringstream ss;
      ss << (i+1);
      edgeLabels.append(ss.str());
    }
  } 

  po::options_description statistics_options
    ("Graph statistics");
  
  statistics_options.add_options()
    ("all-stats",
     "print all statistics")
    ("degree-dist",
     "compute degree distribution")
    ("pairs",
     "count pairs")
    ("triples",
     "count triples")
    ("local-cluster-coeff",
     "compute averaged local clustering coefficients")
    ("global-cluster-coeff",
     "compute global clustering coefficients (from network)")
    ("assortativity",
     "compute assortativities")
    ("degree-overlap",
     "compute degree overlap")
    ("path-length",
     "compute shortest path lengths")
    ("print-degrees",
     "print the degrees of vertices")
    ("cluster-coeff",
     "write clustering coefficients to baseName.cluster file (from adjacency matrices)")
    ("write-Js",
     "write adjacency matrices to files  baseName.Jd/i")        
    ("community",
     "determine community structure")        
    ("modularity", 
     "determine modularity")        
    ("components", 
     "determine component distribution")        
    ;
  
  po::options_description graph_options
    ("Graph options");
  
  graph_options.add_options()
    ("vertices,N", po::value<unsigned int>(),
     "number of vertices")
    ("read-file,r", po::value<std::string>(),
     "read full graph from file (ignores topology options")
    ;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::string prefix("");
    if (nEdgeTypes > 1) {
      prefix = std::string(1,edgeLabels[i]) + "-"; 
    }
    graph_options.add_options()
      ((prefix + "topology").c_str(),
       po::value<std::string>(),
       (prefix + "network topology\n((tri-)lattice, "+
        "tree, random, random-regular, small-world, plod, albert-barabasi, "+
        "community, complete, read, null)").c_str());
  }
  
  po::options_description assortativity_options
    ("Assortativity options");
  
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    for (unsigned int j = 0; j < nEdgeTypes; ++j) {
      for (unsigned int k = j; k < nEdgeTypes; ++k) {
        std::string prefix("");
        std::string desc("");
        if (nEdgeTypes > 1) {
          prefix =
            std::string(1,edgeLabels[i]) + "-" +
            std::string(1,edgeLabels[j]) +
            std::string(1,edgeLabels[k]) + "-";
          desc =
            " between "+ std::string(1,edgeLabels[j]) + "- and " +
            std::string(1,edgeLabels[k]) + "- along " +
            std::string(1,edgeLabels[i]) +"-edges";
        }
        assortativity_options.add_options()
          ((prefix + "assortativity").c_str(),
           po::value<double>(),
           ("degree assortativity" + desc).c_str());
      }
    }
  }

  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    for (unsigned int j = i+1; j < nEdgeTypes; ++j) {
      assortativity_options.add_options()
        ((std::string(1,edgeLabels[i]) + std::string(1,edgeLabels[j]) +
          "-degree-overlap").c_str(), po::value<double>(),
         ("degree overlap between " + std::string(1,edgeLabels[i]) + "- and " +
          std::string(1,edgeLabels[j]) + "-edges").c_str());
    }
  }

  // generate topology-specific options for each type of graph
  
  // lattice
  std::vector<po::options_description*> lattice_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "lattice options";
    po::options_description* lo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "dim";
    lo->add_options()
      (s.str().c_str(),  po::value<unsigned int>()->default_value(2),
       "number of dimensions");
    s.str("");
    s << prefix.str() << "pb";
    lo->add_options()
      (s.str().c_str(), "periodic boundary conditions");
    lattice_options.push_back(lo);
  }

  // tree
  std::vector<po::options_description*> tree_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "tree options";
    po::options_description* to =
      new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "branches";
    to->add_options()
      (s.str().c_str(),  po::value<unsigned int>(),
       "number of branches to create per vertex");
    tree_options.push_back(to);
  }

  // random graph
  std::vector<po::options_description*> rg_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "-random graph options";
    po::options_description* ro
      = new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "edges";
    ro->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "number of edges");
    s.str("");
    s << prefix.str() << "avg-degree";
    ro->add_options()
      (s.str().c_str(), po::value<double>(),
       "average degree");
    rg_options.push_back(ro);
  }

  // random regular graph
  std::vector<po::options_description*> rrg_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "random regular graph options";
    po::options_description* rrgo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "node-degree";
    rrgo->add_options()
      (s.str().c_str(), po::value<double>(),
       "degree of the graph G(n,d)");
    s.str("");
    s << prefix.str() << "joint-degree";
    rrgo->add_options()
      (s.str().c_str(), po::value<unsigned int>(),
       "degree of the joint core of two random regular graphs. If not given, the overlapping is random.");
    rrg_options.push_back(rrgo);
  }
  
  // small-world graph
  std::vector<po::options_description*> sw_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "small world options";
    po::options_description* swo
      = new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "node-degree";
    swo->add_options()
      (s.str().c_str(), po::value<double>(),
       "degree of each node");
    sw_options.push_back(swo);
  }

  // plod graph
  std::vector<po::options_description*> plod_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "power law out degree options";
    po::options_description* plodo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "alpha";
    plodo->add_options()
      (s.str().c_str(), po::value<double>(),
       "alpha (index of power law)");
    s.str("");
    s << prefix.str() << "beta";
    plodo->add_options()
      (s.str().c_str(), po::value<double>(),
       "beta (multiplicative factor of power law)");
    plod_options.push_back(plodo);
  }

  // Albert-Barabasi graph
  std::vector<po::options_description*> ab_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "albert barabasi options";
    po::options_description* abo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "newedges";
    abo->add_options()
      (s.str().c_str(), po::value<unsigned int>()->default_value(1),
       "number of edges to add per vertex");
    ab_options.push_back(abo);
  }

  // Community graph
  std::vector<po::options_description*> com_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "community options";
    po::options_description* co =
      new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "delta";
    co->add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.5),
       "amount to increase link weight by if strenghtened");
    s.str("");
    s << prefix.str() << "pl";
    co->add_options()
      (s.str().c_str(), po::value<double>()->default_value(5e-4),
       "probability for local attachment");
    s.str("");
    s << prefix.str() << "pr";
    co->add_options()
      (s.str().c_str(), po::value<double>()->default_value(5e-4),
       "probability for creation of random link");
    s.str("");
    s << prefix.str() << "pd";
    co->add_options()
      (s.str().c_str(), po::value<double>()->default_value(1e-3),
       "probabiliy for node deletion");
    s.str("");
    s << prefix.str() << "iter";
    co->add_options()
      (s.str().c_str(), po::value<unsigned int>()->default_value(25),
       "number of iterations");
    com_options.push_back(co);
  }

  // read graph
  std::vector<po::options_description*> readFile_options;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "read options";
    po::options_description* rfo =
      new po::options_description(s.str().c_str());
    s.str("");
    s << prefix.str() << "file";
    rfo->add_options()
      (s.str().c_str(), po::value<std::string>(),
       "name of graph file to read");
    s.str("");
    readFile_options.push_back(rfo);
  }

  // rewiring options
  po::options_description rewiring_options("Rewiring options");
  
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "rewire";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to rewire");
    s.str("");
    s << prefix.str() << "rewire-clustered";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to rewire and cluster");
    s.str("");
    s << prefix.str() << "remove";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to remove");
    s.str("");
    s << prefix.str() << "add";
    rewiring_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.),
       "fraction of edges to add");
  }

  // rewiring options
  po::options_description additional_options("Additional options");
  
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    std::stringstream s, prefix;
    if (nEdgeTypes > 1) {
      prefix << edgeLabels[i] << "-";
    }
    s << prefix.str() << "randomise";
    additional_options.add_options()
      (s.str().c_str(), po::value<double>()->default_value(0.), 
       ("randomise fraction of vertices before adding "+prefix.str()+
        "edges").c_str());
  }

  // read options from command line
  po::options_description all_options;
  all_options.add(main_options).add(edgetype_options).add(output_options).
    add(graph_options).add(rewiring_options).add(additional_options).
    add(assortativity_options).add(statistics_options);
                                   
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    all_options.add(*(lattice_options[i]));
    all_options.add(*(tree_options[i]));
    all_options.add(*(rg_options[i]));
    all_options.add(*(rrg_options[i]));
    all_options.add(*(sw_options[i]));
    all_options.add(*(plod_options[i]));
    all_options.add(*(ab_options[i]));
    all_options.add(*(com_options[i]));
    all_options.add(*(readFile_options[i]));
  }

  std::vector<std::string> unregistered;
  
  try {
    po::parsed_options parsed = po::command_line_parser(argc, argv).
      options(all_options).allow_unregistered().run();
    po::store(parsed, vm);
    unregistered =
      po::collect_unrecognized(parsed.options, po::exclude_positional);
  }
  catch (std::exception& e) {
    std::cerr << "Error parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }

  po::notify(vm);
  for (std::vector<std::string>::iterator it = unregistered.begin();
       it != unregistered.end(); it++) {
    std::cerr << "WARNING: ignoring unknown option " << *it << std::endl;
  }
  
  if (vm.count("longhelp")) {
    std::cout << all_options << std::endl;
    return 0;
  }

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
  multitype_graph graph;
            
  /******************************************************************/
  // read graph from file or generate it
  /******************************************************************/
  
  // adding N vertices to graph
  boost::add_vertices(graph, N);
  std::vector<std::string>* vertexOptions = 0;

  /******************************************************************/
  // generate edges
  /******************************************************************/

  std::vector<unsigned int> edgeTypes;
  for (unsigned int i = 0; i < nEdgeTypes; ++i) edgeTypes.push_back(i);

  if (vm.count("read-file")) {
    readGraph = vm["read-file"].as<std::string>();
    readAll = true;
  }

  for (unsigned int i = 0; i < edgeTypes.size(); ++i) {

    onetype_graph temp_graph;
    boost::add_vertices(temp_graph, N);

    std::string currentEdgeLabel("");
    
    if (nEdgeTypes > 1) {
      currentEdgeLabel = std::string(1, edgeLabels[edgeTypes[i]]) + "-";
    }
    if (verbose >= 2) {
      std::cout << "Creating " << currentEdgeLabel << "graph" << std::endl;
    }
    std::string optStr = currentEdgeLabel + "topology";
    if (readAll) {
      topology = "read";
    } else if (vm.count(optStr)) {
      topology = vm[optStr].as<std::string>();
    } else {
      std::cerr << "ERROR: no " << optStr << " specified" << std::endl;
      std::cerr << main_options << graph_options << std::endl;
      return 1;
    }
    
    if (topology == "lattice") {

      /******************************************************************/
      // read lattice specific parameters
      /******************************************************************/
      
      latticeOptions opt;
      optStr = currentEdgeLabel + "dim";
      opt.dimensions = vm[optStr].as<unsigned int>();
      opt.sideLength = static_cast<int>(pow(N, 1.0/opt.dimensions));
      if (pow(opt.sideLength, opt.dimensions) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }
      optStr = currentEdgeLabel + "pb";
      if (vm.count(optStr)) {
        opt.periodicBoundary = true;
      } else {
        opt.periodicBoundary = false;
      }
      
      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      typedef boost::lattice_iterator<onetype_graph> lattice_iterator;
      
      lattice_iterator li(opt.sideLength,opt.dimensions,opt.periodicBoundary);
      lattice_iterator li_end;
      
      boost::add_edge_structure(temp_graph, li, li_end, Edge(i));
      
    } else if (topology == "tri-lattice") {
      
      /******************************************************************/
      // read tri-lattice specific parameters
      /******************************************************************/
      
      latticeOptions opt;
      
      opt.sideLength = static_cast<int>(sqrt(N));
      if (pow(opt.sideLength, 2) != N) { 
        std::cerr << "ERROR: cannot generate square lattice out of " << N
                  << " vertices." << std::endl;
        return 1;
      }
      
      optStr = std::string(1, edgeLabels[i]) + std::string("-pb");
      if (vm.count(optStr)) {
        opt.periodicBoundary = true;
      } else {
        opt.periodicBoundary = false;
      }
      
      /******************************************************************/
      // generate lattice with desired properties
      /******************************************************************/
      typedef boost::tri_lattice_iterator<onetype_graph> tri_lattice_iterator;
      
      tri_lattice_iterator tli(opt.sideLength,opt.periodicBoundary);
      tri_lattice_iterator tli_end;
      
      boost::add_edge_structure(temp_graph, tli, tli_end, Edge(i));

    } else if (topology == "tree") {
      
      /******************************************************************/
      // read tree specific parameters
      /******************************************************************/
      
      treeOptions opt;

      optStr = currentEdgeLabel + "branches";
      if (vm.count(optStr)) {
        opt.branches = vm[optStr].as<unsigned int>();
      } else {
        std::cerr << "ERROR: number of branches not specified" << std::endl;
        std::cerr << *tree_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate tree with desired properties
      /******************************************************************/
      typedef boost::tree_iterator<onetype_graph> tree_iterator;
      
      tree_iterator ti(N, opt.branches);
      tree_iterator ti_end;
      
      boost::add_edge_structure(temp_graph, ti, ti_end, Edge(i));
      
    } else if (topology == "random") {

      /******************************************************************/
      // read random graph specific parameters
      /******************************************************************/
      
      rgOptions opt;

      optStr = currentEdgeLabel + "edges";
      std::string optStr2 = currentEdgeLabel + "avg-degree";
      if (vm.count(optStr)) {
        opt.edges = vm[optStr].as<unsigned int>();
        if (vm.count(optStr2)) {
          std::cerr << "WARNING: " << optStr << " and " << optStr2
                     << " exclude each other, ignoring " << optStr2
                     << std::endl;
        }
      } else {
        if (vm.count(optStr2)) {
          opt.edges =
            static_cast<unsigned int>(vm[optStr2].as<double>() * N / 2.);
        } else {
          std::cerr << "ERROR: neither " << optStr << " nor "
                    << optStr2 << " specified." << std::endl;
          std::cerr << *rg_options[i] << std::endl;
          return 1;
        }
      }
      
      /******************************************************************/
      // generate random graph with desired properties
      /******************************************************************/
//       typedef boost::erdos_renyi_iterator2<boost::mt19937, dualtype_graph>
//         rg_iterator;
      
      typedef boost::erdos_renyi_iterator<boost::mt19937, onetype_graph>
        rg_iterator;
      
//       double p = (double)opt.edges*2.0/((double)N*(double)(N-1));
      double p = (double)opt.edges*2.0/((double)N*(double)N);
      if (p> 1) p=1;
      
      rg_iterator ri(gen, N, p);
      rg_iterator ri_end;
      
      boost::add_edge_structure(temp_graph, ri, ri_end, Edge(i));

    } else if (topology == "random-regular") {
      
      /******************************************************************/
      // read random regular graph specific parameters
      /******************************************************************/
      
      rrgOptions opt;

      optStr = currentEdgeLabel + "node-degree";
      if (vm.count(optStr)) {
        opt.degree = static_cast<unsigned int>(vm[optStr].as<double>());
      } 
      optStr = currentEdgeLabel + "joint-degree";
      if (vm.count(optStr)) {
        opt.jointDegree = vm[optStr].as<unsigned int>();
      }
      if (opt.degree + opt.jointDegree == 0) {
        std::cerr << "WARNING: Neither degree nor joint-degree specified, "
                  << " for random-regular " << edgeLabels[i] << "-graph."
                  << std::endl;
        std::cerr << *rrg_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate random regular graph with desired properties
      /******************************************************************/
      
      bool success = false;
      unsigned int count = 0;

      typedef std::vector<std::pair<unsigned int, unsigned int> > GraphEdges;
      GraphEdges rrg_edges;

      if (opt.jointDegree > 0) {
        // generate joint graph
        rrg_edges.clear();         
        ++count;
        
        while (!success) {
          success = boost::random_regular_graph(graph, rrg_edges,
                                                opt.jointDegree, N, gen);
          if (success) {            
            for (GraphEdges::iterator it = rrg_edges.begin();
                 it != rrg_edges.end(); ++it){
              boost::add_edge((*it).first, (*it).second,
                              Edge(i), temp_graph);
            }
          }
        }
        
        if (verbose) {
          std::cout << "Random regular graph of degree " << opt.jointDegree
                    << " was generated in " << count << " trials\n";
        }

        // copy graph to all edgetypes
        for (unsigned int j = 0; j < edgeTypes.size(); j++) {
          boost::copy_edges(temp_graph, graph, Edge(j));
        }
        temp_graph.clear();
      }

      if (opt.degree > 0) {
        // generate additional graph, excluding existing edges
        success = 0;
        count = 0;
        rrg_edges.clear();
        ++count;
        while (!success) {
          success = boost::random_regular_graph(graph, rrg_edges,
                                                opt.degree, N, gen);
          if (success) {            
            for (GraphEdges::iterator it = rrg_edges.begin();
                 it != rrg_edges.end(); ++it){
              boost::add_edge((*it).first, (*it).second,
                              Edge(i), temp_graph);
            }
          }
        }
        
        if (verbose) {
          std::cout << "Random regular graph of degree " << opt.degree
                    << " was generated in " << count << " trials\n";
        }
      }
      
    } else if (topology == "small-world") {
      
      /******************************************************************/
      // read small-world graph specific parameters
      /******************************************************************/
      
      swOptions opt;
      
      optStr = currentEdgeLabel + "node-degree";
      if (vm.count(optStr)) {
        opt.degree = static_cast<unsigned int>(vm[optStr].as<double>());
      } else {
        std::cerr << "ERROR: no degree specified" << std::endl;
        std::cerr << *sw_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate small-world graph with desired properties
      /******************************************************************/
      typedef boost::small_world_iterator<boost::mt19937, onetype_graph>
        sw_iterator;
      
      sw_iterator swi(gen,N,opt.degree,0);
      sw_iterator swi_end;
      
      boost::add_edge_structure(temp_graph, swi, swi_end, Edge(i));
      
    } else if (topology == "plod") {
      
      /******************************************************************/
      // read plod graph specific parameters
      /******************************************************************/
      
      plodOptions opt;
      
      optStr = currentEdgeLabel + "alpha";
      if (vm.count(optStr)) {
        opt.alpha = vm[optStr].as<double>();
      } else {
        std::cerr << "ERROR: no alpha specified" << std::endl;
        std::cerr << *plod_options[i] << std::endl;
        return 1;
      }
      optStr = currentEdgeLabel + "beta";
      if (vm.count(optStr)) {
        opt.beta = vm[optStr].as<double>();
      } else {
        std::cerr << "ERROR: no beta specified" << std::endl;
        std::cerr << *plod_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate plod graph with desired properties
      /******************************************************************/
      typedef boost::plod_iterator<boost::mt19937, onetype_graph>
        plod_iterator;
      
      plod_iterator plodi(gen,N,opt.alpha,opt.beta);
      plod_iterator plodi_end;
      
      boost::add_edge_structure(temp_graph, plodi, plodi_end, Edge(i));

    } else if (topology == "albert-barabasi") {
      
      /******************************************************************/
      // read Albert-Barabasi graph specific parameters
      /******************************************************************/
      
      abOptions opt;
      
      optStr = currentEdgeLabel + "newedges";
      if (vm.count(optStr)) {
        opt.new_edges = vm[optStr].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no number of edges to add per vertex specified"
                  << std::endl;
        std::cerr << *ab_options[i] << std::endl;
        return 1;
      }
      
      /******************************************************************/
      // generate Albert-Barabasi graph with desired properties
      /******************************************************************/
      typedef boost::albert_barabasi_iterator<boost::mt19937, onetype_graph>
        albert_barabasi_iterator;
      
      albert_barabasi_iterator ab(gen,N,opt.new_edges);
      albert_barabasi_iterator ab_end;
      
      boost::add_edge_structure(temp_graph, ab, ab_end, Edge(i));

    } else if (topology == "community") {
      
      /******************************************************************/
      // read community graph specific parameters
      /******************************************************************/

      communityOptions opt;

      optStr = currentEdgeLabel + "delta";
      if (vm.count(optStr)) {
        opt.delta = vm[optStr].as<double>();
      } else {
        std::cerr << "ERROR: no " << optStr << " specified"
                  << std::endl;
        std::cerr << *com_options[i] << std::endl;
        return 1;
      }
      optStr = currentEdgeLabel + "pl";
      if (vm.count(optStr)) {
        opt.pl = vm[optStr].as<double>();
      } else {
        std::cerr << "ERROR: no " << optStr << " specified"
                  << std::endl;
        std::cerr << *com_options[i] << std::endl;
        return 1;
      }
      optStr = currentEdgeLabel + "pr";
      if (vm.count(optStr)) {
        opt.pr = vm[optStr].as<double>();
      } else {
        std::cerr << "ERROR: no " << optStr << " specified"
                  << std::endl;
        std::cerr << *com_options[i] << std::endl;
        return 1;
      }
      optStr = currentEdgeLabel + "pd";
      if (vm.count(optStr)) {
        opt.pd = vm[optStr].as<double>();
      } else {
        std::cerr << "ERROR: no " << optStr << " specified"
                  << std::endl;
        std::cerr << *com_options[i] << std::endl;
        return 1;
      }
      optStr = currentEdgeLabel + "iter";
      if (vm.count(optStr)) {
        opt.iterations = vm[optStr].as<unsigned int>();
      } else {
        std::cerr << "ERROR: no " << optStr << " specified"
                  << std::endl;
        std::cerr << *com_options[i] << std::endl;
        return 1;
      }

      /******************************************************************/
      // generate community graph with desired properties
      /******************************************************************/
      typedef boost::community_iterator<boost::mt19937, onetype_graph>
        community_iterator;
      
      community_iterator com(gen, &temp_graph, opt.delta, opt.pl, opt.pr, opt.pd,
                             opt.iterations);
      community_iterator com_end;
      
      boost::add_edge_structure(temp_graph, com, com_end, Edge(i));

    } else if (topology == "complete") {

      boost::graph_traits<multitype_graph>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; vi++) {
        boost::graph_traits<multitype_graph>::vertex_iterator vi2;
        for (vi2 = vi+1; vi2 != vi_end; vi2++) {
          add_edge(*vi, *vi2, Edge(i), temp_graph);
        }
      }
      
    } else if (topology == "copy") {
      
      /******************************************************************/
      // test if we have graph to copy from
      /******************************************************************/

      if (num_edges(graph) > 0) {
        boost::copy_edges(graph, temp_graph, Edge(i));
      } else {
        if (edgeTypes.size() == (2*nEdgeTypes-1) || i == (edgeTypes.size()-1)) {
          std::cerr << "ERROR: no edges to copy from, " << edgeTypes[i]
                    << "-graph to remain empty" << std::endl;
          std::cerr << main_options << graph_options << std::endl;
          return 1;
        } else {
          // push back for later
          edgeTypes.push_back(edgeTypes[i]);
	}
      }
    } else if (topology == "read") {
      
      /******************************************************************/
      // read file graph specific parameters
      /******************************************************************/
      readFileOptions opt;

      optStr = currentEdgeLabel + "file";
      if (vm.count(optStr)) {
        if (readGraph.size() == 0) {
          readGraph = vm[optStr].as<std::string>();
        }
        opt.fileName = vm[optStr].as<std::string>();
      } else {
        if (readGraph.size() == 0) {
          std::cerr << "ERROR: no file name specified" << std::endl;
          std::cerr << *readFile_options[i] << std::endl;
          return 1;
        } else {
          opt.fileName = readGraph;
        }
      }

      // reading graph structure
      if (vertexOptions == 0) vertexOptions = new std::vector<std::string>(N);
      int read_result = read_graph(temp_graph, opt.fileName, i,
                                   vertexOptions);
      if (read_result > 0) {
        // update number of vertices
        N = num_vertices(graph);
        if (verbose) {
          std::cout << currentEdgeLabel << "graph file " << opt.fileName 
	            << " was read ok" << std::endl;
        }
      } else if (read_result == 0) {
        std::cerr << "WARNING: no " << edgeLabels[i] << "-edges in "
                  << opt.fileName << std::endl;
      } else {
        std::cerr << "ERROR: something wrong in reading "
	          << edgeLabels[i] << "-graph from "
                  << opt.fileName << std::endl;
        return 1;
      }
    } else if (topology == "null") {
    } else {
      std::cerr << "ERROR: unknown " << optStr << ": " << topology
                << std::endl;
      std::cerr << main_options << graph_options << std::endl;
      return 1;
    }

    // do graph modifications if desired
    if (num_edges(temp_graph) > 0 &&
        vm.count(currentEdgeLabel+"remove")) {
      double removeFraction =
        vm[currentEdgeLabel+"remove"].as<double>();
      if (removeFraction > 0) {
        boost::removeEdges(temp_graph, gen, removeFraction);
      }
    }
    if (num_edges(temp_graph) > 0 &&
        vm.count(currentEdgeLabel+"rewire")) {
      double rewireFraction =
        vm[currentEdgeLabel+"rewire"].as<double>();
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
    if (num_edges(temp_graph) > 0 &&
        vm.count(currentEdgeLabel+"rewire-clustered")) {
      double clusteredRewireFraction =
        vm[currentEdgeLabel+"rewire-clustered"].as<double>();
      if (clusteredRewireFraction > 0) {
//        onetype_graph rewire_graph = temp_graph;
        boost::rewireClustered(temp_graph, temp_graph, gen,
                               clusteredRewireFraction, verbose);
//        temp_graph = rewire_graph;
        if (verbose) {
          std::cout << "graph rewired\n";
        }
      }
    }
    if (num_edges(temp_graph) > 0 && vm.count(currentEdgeLabel+"add")) {
      double addFraction = vm[currentEdgeLabel+"add"].as<double>();
      boost::addEdges(temp_graph, gen, addFraction);
    }
    
    if (vm.count(currentEdgeLabel+"randomise")) {
      double randFraction = vm[currentEdgeLabel+"randomise"].as<double>();
      randomise_vertices(temp_graph, gen, randFraction);
    }

    // copy edges to main graph
    boost::copy_edges(temp_graph, graph);
  }

  // checking graph  
  if (num_vertices(graph) == 0) {
    std::cerr << "ERROR: no vertices" << std::endl;
    std::cerr << main_options << graph_options << std::endl;
    return 1;
  }

  /******************************************************************/
  // create given assortativity if desired
  /******************************************************************/
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    for (unsigned int j = 0; j < nEdgeTypes; ++j) {
      for (unsigned int k = j; k < nEdgeTypes; ++k) {
        std::string s("");
        if (nEdgeTypes > 1) {
          s = std::string(1,edgeLabels[i]) + "-" +
            std::string(1,edgeLabels[j]) + std::string(1,edgeLabels[k]) +
            "-";
        }
        s = s + "assortativity";
        if (vm.count(s.c_str())) {
          double ass = (vm[s.c_str()].as<double>());
          rewire_assortatively(graph, gen, ass, i, j, k, verbose);
        }
      }
    }
  }
  
  /******************************************************************/
  // create given degree overlap if desired
  /******************************************************************/

  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    for (unsigned int j = i+1; j < nEdgeTypes; ++j) {
      std::string s = std::string(1,edgeLabels[i]) +
        std::string(1,edgeLabels[j]) + "-degree-overlap";
      if (vm.count(s.c_str())) {
        double deg_overlap = (vm[s.c_str()].as<double>());
        rewire_degree_overlap(graph, gen, deg_overlap, i, j, verbose);
      }
    }
  }

  /********* graph creation done **********************************/

  if (verbose) {
    std::cout << "\nGraph created." << std::endl;
  }

  setup_edge_index_map(graph);

  /******************************************************************/
  // Mark parallel edges
  /******************************************************************/
  unsigned int parallel_edges = 0;
  if (nEdgeTypes > 1) {
    parallel_edges = mark_parallel_edges(graph);
    if (verbose >= 2) {
      std::cout << "No. of parallel edges is: " << parallel_edges
                << std::endl;
    }
  }

  /******************************************************************/
  // Write graph files
  /******************************************************************/
  std::string baseFileName("");
  if (vm.count("output-file")) {
    baseFileName = vm["output-file"].as<std::string>();
  } else if (readAll && (vm.count("split") || !verbose)) {
    baseFileName = readGraph.substr(0,readGraph.rfind("."));
  }

  std::string ext;
  if (baseFileName.size() > 0) {
    ext = ".graph";
    if (vm.count("split")) {
      // split graph into several files
      onetype_graph split_graph;
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        split_graph.clear();
        copy_edges_type(graph, split_graph, Edge(i));
        std::string outputGraphName = baseFileName + "_" +
          std::string(1, edgeLabels[i]) + ext;
        write_graph(split_graph, outputGraphName, false,
                    vertexOptions);
      }
    } else if (!readAll) {
      // write whole graph in one file
      std::string outputGraphName =
        baseFileName+".graph";
      write_graph(graph, outputGraphName, true,
                  vertexOptions);
    }
  
    // create sparse adjacency matrices and clustering coefficients
    if (vm.count("cluster-coeff")) {
      bool writeJs = false;
      if (vm.count("write-Js")) writeJs = true;
      
      bool status = boost::cluster_coeff(graph, baseFileName, writeJs);
      if (verbose) {
        if (!status) {
          std::cout << "cluster coeff file was written ok\n";
        } else {
          std::cout << "ERROR: something wrong in writing cluster coeff file "
                    << std::endl;
        }
      }
    }
  }

  // calculate community structure
  if (vm.count("community") || vm.count("modularity")) {
    multitype_graph saved_graph = graph;
    double mod =
      boost::community_structure(saved_graph, graph, 1., verbose,
                                 vm.count("community"),
                                 vm.count("modularity") || verbose); 
    if (vm.count("community") && baseFileName.length() > 0) {
      write_graph(graph, baseFileName+".comm"+ext,
                  true, vertexOptions);
    }
    if (vm.count("modularity")) {
      std::stringstream output;
      output << "  M = " << mod << std::endl;
      if (baseFileName.length() == 0 || verbose) {
        std::cout << "\nModularity:" << std::endl;
        std::cout << output.str();
      }
      if (baseFileName.length() > 0) {
        std::string modFileName = baseFileName + ".stat.mod";
        std::ofstream modFile(modFileName.c_str(), std::ios::out);
        modFile << output.str();
        modFile.close();
      }
    }
    graph = saved_graph;
  }

  if (vm.count("components") && baseFileName.length() > 0) {
    write_component_dist(graph, baseFileName+".stat.comp");
  }

  if (vm.count("all-stats")) {
    allStats = true;
//    if (verbose == 0) verbose = 1;
  }

  /******************************************************************/
  // count vertices
  /******************************************************************/
  if (vm.count("pairs") || vm.count("triples") || verbose || allStats) {
    std::stringstream output;
    if (baseFileName.length() == 0 || verbose) {
      output << "\nVertex count:" << std::endl;
    }
    output << "  N: " << num_vertices(graph) << std::endl;;
    
    /******************************************************************/
    // count pairs
    /******************************************************************/
    if (vm.count("pairs") || verbose || allStats) {
      
      if (baseFileName.length() == 0 || verbose) {
        output << "\nPair count:" << std::endl;
      }
      std::vector<unsigned int> pairs = count_pairs(graph, nEdgeTypes);
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        std::string prefix;
        if (nEdgeTypes > 1) {
          prefix = std::string(1,edgeLabels[i]) + "-";
        }
        output << "  " << prefix << "pairs: " << pairs[i]
               << " (avg degree: " << pairs[i]*2./num_vertices(graph)
	       << ")" << std::endl;
      }
      if (nEdgeTypes > 1) {
        output << " " << "parallel: " << parallel_edges << std::endl;
      }
    }
    
    /******************************************************************/
    // count triples
    /******************************************************************/
    if (vm.count("triples") || allStats) {
      if (baseFileName.length() == 0 || verbose) {
        output << "\nTriple count:" << std::endl;
      }
      boost::multi_array<unsigned int, 2> triples =
        count_triples(graph, nEdgeTypes);
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        for (unsigned int j = i; j < nEdgeTypes; ++j) {
          std::string prefix;
          if (nEdgeTypes > 1) {
            prefix = std::string(1,edgeLabels[i]) + 
	             std::string(1,edgeLabels[j]) + "-";
          }
          output << "  " << prefix << "triples: "
                 << triples[i][j] << std::endl;
        }
      }
    }

    if (baseFileName.length() == 0 || verbose) {
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string countFileName = baseFileName+".stat.counts";
      std::ofstream countFile(countFileName.c_str(), std::ios::out);
      countFile << output.str();
      countFile.close();
    }
  }

  /******************************************************************/
  // calculate degree distribution
  /******************************************************************/
  if (vm.count("degree-dist") || allStats) {
    std::stringstream output;

    boost::multi_array<unsigned int, 2> dd;
    unsigned int max_degree = degree_dist(graph, nEdgeTypes, dd);

    std::stringstream s;
    s << num_vertices(graph);
    int len = std::max(6, static_cast<int>(s.str().length()));

    output << "# " << "degree" << "\t";
    if (nEdgeTypes > 1) {
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        output << edgeLabels[i] << " only" << "\t";
      }
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        output << edgeLabels[i] << " all " << "\t";
      }
      output << "parall" << "\t" << " all  " << std::endl;
    } else {
      output << " nodes" << std::endl;
    }
    
    for (unsigned int i = 0; i < max_degree+1; ++i) {

      output << "  " << std::setw(len) << i << "\t";
      if (nEdgeTypes > 1) {
        for (unsigned int j = 0; j < nEdgeTypes*2+2; ++j) {
          output << std::setw(len) << dd[j][i] << "\t";
        }
      } else {
        output << std::setw(len) << dd[0][i];
      }
      output << std::endl;
    }

    if (baseFileName.length() == 0 || verbose) {
      std::cout << "\nDegree distribution:" << std::endl;
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string degreeFileName = baseFileName+".stat.degree";
      std::ofstream degreeFile(degreeFileName.c_str(), std::ios::out);
      degreeFile << output.str();
      degreeFile.close();
    }
  }

  /******************************************************************/
  // calculate local clustering coefficients
  /******************************************************************/
  if (vm.count("local-cluster-coeff")) {
    std::stringstream output;
    
    boost::multi_array<double, 3> lcc =
      local_cluster_coeffs(graph, nEdgeTypes);

    for (unsigned int i = 0; i< nEdgeTypes; ++i) {
      for (unsigned int j = i; j< nEdgeTypes; ++j) {
        for (unsigned int k = 0; k< nEdgeTypes; ++k) {
          std::string prefix;
          if (nEdgeTypes > 1) {
            prefix = std::string(1,edgeLabels[i]) + 
	             std::string(1,edgeLabels[j]) + 
		     std::string(1,edgeLabels[k]);
          }
          output << "  C" << prefix << " = " << lcc[i][j][k] << std::endl;
        }
      }
    }

    if (baseFileName.length() == 0 || verbose) {
      std::cout << "\nLocal clustering coeffients:" << std::endl;
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string lccFileName = baseFileName + ".stat.lcc";
      std::ofstream lccFile(lccFileName.c_str(), std::ios::out);
      lccFile << output.str();
      lccFile.close();
    }
  }

  /******************************************************************/
  // calculate global clustering coefficients
  /******************************************************************/
  if (vm.count("global-cluster-coeff") || allStats) {
    std::stringstream output;
    
    boost::multi_array<double, 3> gcc =
      global_cluster_coeffs(graph, nEdgeTypes);

    for (unsigned int i = 0; i< nEdgeTypes; ++i) {
      for (unsigned int j = i; j< nEdgeTypes; ++j) {
        for (unsigned int k = 0; k< nEdgeTypes; ++k) {
          std::string prefix;
          if (nEdgeTypes > 1) {
            prefix = std::string(1,edgeLabels[i]) + 
	             std::string(1,edgeLabels[j]) + 
		     std::string(1,edgeLabels[k]);
          }
          output << "  C" << prefix << " = " << gcc[i][j][k] << std::endl;
        }
      }
    }

    if (baseFileName.length() == 0 || verbose) {
      std::cout << "\nGlobal clustering coeffients:" << std::endl;
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string gccFileName = baseFileName + ".stat.gcc";
      std::ofstream gccFile(gccFileName.c_str(), std::ios::out);
      gccFile << output.str();
      gccFile.close();
    }
  }

  /******************************************************************/
  // calculate average path lengths
  /******************************************************************/
//  if (vm.count("path-length") || allStats) {
  if (vm.count("path-length")) {
    std::stringstream output;
    
    for (unsigned int i = 0; i < nEdgeTypes; ++i) {
      std::string prefix;
      if (nEdgeTypes > 1) {
        prefix = std::string(1,edgeLabels[i]) + "-";
      }
      output << "  on " << prefix << "edges: "
             << boost::avg_shortest_path_length(graph, Edge(i)) << std::endl;
      for (unsigned int j = 0; j < nEdgeTypes; ++j) {
        if (i != j) {
          output << "  " << std::string(1,edgeLabels[j]) << "-neighbours on "
                 << std::string(1,edgeLabels[i]) << "-edges: "
                 << boost::avg_nb_shortest_path_length(graph, Edge(j), Edge(i))
                 << std::endl;
        }
      }
    }

    if (baseFileName.length() == 0 || verbose) {
      std::cout << "\nAverage shortest path lengths: " << std::endl;
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string splFileName = baseFileName + ".stat.spl";
      std::ofstream splFile(splFileName.c_str(), std::ios::out);
      splFile << output.str();
      splFile.close();
    }
  }

  /******************************************************************/
  // calculate assortativities 
  /******************************************************************/
  if (vm.count("assortativity") || allStats) {
    std::stringstream output;
    
    for (unsigned int i = 0; i < nEdgeTypes; ++i) {
      for (unsigned int j = 0; j < nEdgeTypes; ++j) {
        for (unsigned int k = j; k < nEdgeTypes; ++k) {
          std::string prefix;
          if (nEdgeTypes > 1) {
            prefix = std::string(1,edgeLabels[i]) + "-" + 
	             std::string(1,edgeLabels[j]) +
                     std::string(1,edgeLabels[k]) + "-";
          }
          output << "  " << prefix << "assortativity: "
                 << boost::assortativity(graph, i, j, k)
                 << std::endl;
        }
      }
    }
    if (baseFileName.length() == 0 || verbose) {
      std::cout << "\nAssortatitivies: " << std::endl;
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string assFileName = baseFileName + ".stat.ass";
      std::ofstream assFile(assFileName.c_str(), std::ios::out);
      assFile << output.str();
      assFile.close();
    }
  }

  /******************************************************************/
  // calculate degree overlaps 
  /******************************************************************/
  if (vm.count("degree-overlap") || (allStats && nEdgeTypes > 1)) {
    std::stringstream output;
    
    for (unsigned int i = 0; i < nEdgeTypes; ++i) {
      for (unsigned int j = i+1; j < nEdgeTypes; ++j) {
        output << "  " << std::string(1,edgeLabels[i])
               << std::string(1,edgeLabels[j])
               << "-degree-overlap: "
               << boost::degree_overlap(graph, i, j)
               << std::endl;
      }
    }

    if ((baseFileName.length() == 0 || verbose) && nEdgeTypes > 1) {
      std::cout << "\nDegree overlap: " << std::endl;
      std::cout << output.str();
    }
    if (baseFileName.length() > 0) {
      std::string dolFileName = baseFileName + ".stat.dol";
      std::ofstream dolFile(dolFileName.c_str(), std::ios::out);
      dolFile << output.str();
      dolFile.close();
    }
  }

  /******************************************************************/
  // print vertex degrees
  /******************************************************************/
  if (vm.count("print-degrees")) {
    print_degrees(graph, nEdgeTypes);
  }
  
  /******************************************************************/
  // free memory
  /******************************************************************/
  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    delete lattice_options[i];
    delete tree_options[i];
    delete rg_options[i];
    delete rrg_options[i];
    delete sw_options[i];
    delete plod_options[i];
    delete ab_options[i];
    delete com_options[i];
    delete readFile_options[i];
  }
  
  return 0;
  
}
