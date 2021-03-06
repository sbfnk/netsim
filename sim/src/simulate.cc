/*! \file simulate.cc
  \brief The main simulation program.

  This is the main program for running network simulations of disease spread or
  network evolution. It reads one or more graphs from a .dot file, sets or reads
  the initial states of nodes, and runs a selected model according to the
  parameters given to that model.

  The program relies on an object of a Simulator class to run the
  simulation. The object of the Simulator class in turn contains a Model object
  which specifies which model exactly is to be run.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>

#include <boost/filesystem.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>

#include "network/include/graph_structure.hh"
#include "network/include/graph_io.hh"
#include "network/include/nearest_infected.hh"
#include "network/include/Edge.hh"

#include "sim_statistics.hh"

#include "RewireSimulator.hh"
#include "EpiSimulator.hh"
#include "EpiRewireSimulator.hh"
#include "ChrisSimulator.hh"
#include "NetEvoSimulator.hh"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

// define two types of graph: multitype_graph for graphs with multiple types of
// edges, and onetype_graph for graphs with just one type
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> multitype_graph;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> onetype_graph;

int main(int argc, char* argv[])
{
  const Model<multitype_graph>* model = 0; //!< the model to be run
  Simulator<multitype_graph>* sim = 0; //!< the simulator to use

  std::vector <unsigned int> init; //!< initial network state

  /******************************************************************/
  // parameters
  /******************************************************************/
  unsigned int N = 0; //!< number of edges

  char* dataDirEnv = getenv("DATADIR");
  std::string outputDir = "";
  if (dataDirEnv) {
    mkdir(dataDirEnv, 0755);
    outputDir = dataDirEnv;
  }

  unsigned int verbose = 0; //!< Level of verbosity

  unsigned int numSims = 1; //!< Number of simulation runs

  bool clearIC = true; //!< Clear initial conditions
  bool keepIC = false; //!< Keep initial network state equal between runs
  bool randomIC = true; //!< Generate random initial conditions

  bool doIO = false; //!< Write something to disk?

  /******************************************************************/
  // read main command line options
  /******************************************************************/
  po::options_description main_options
    ("\nUsage: simulate [options]... \n\nMain options");

  main_options.add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message including all options")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("sim", po::value<std::string>()->default_value("EpiGillespie"),
     "simulator to use (EpiGillespie, Chris, Rewire, NetEvo)")
    ("output-dir,o",po::value<std::string>(),
     "output directory for graph and data output (as subdir of $DATADIR)")
    ("nsims", po::value<unsigned int>()->default_value(1),
     "number of simulation runs to produce (on a given graph)")
    ;

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(main_options).
              allow_unregistered().run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR parsing command line parameters: " << e.what()
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
    std::cout << main_options << std::endl;
    return 0;
  }

  numSims = vm["nsims"].as<unsigned int>();

  if (numSims < 1) {
    std::cerr << "nsims smaller than one, not running any simulations" << std::endl;
    return 0;
  }

  /******************************************************************/
  // create graph variable
  /******************************************************************/
  multitype_graph graph, saved_graph;

  /******************************************************************/
  // initialise random generator
  /******************************************************************/

  unsigned int seed;
  struct timeval tv;
  gettimeofday(&tv, 0);
  seed = tv.tv_sec + tv.tv_usec;
  boost::mt19937 gen(seed);

  /******************************************************************/
  // create simulator
  /******************************************************************/

  if (vm.count("sim")) {
    std::string simType = vm["sim"].as<std::string>();
    if (simType == "EpiGillespie") {
      sim = new Simulators::EpiSimulator<boost::mt19937, multitype_graph>
        (gen, graph, verbose);
    } else if (simType == "EpiRewire") {
      sim = new Simulators::EpiRewireSimulator<boost::mt19937, multitype_graph>
        (gen, graph, verbose);
    } else if (simType == "Chris") {
      sim = new Simulators::ChrisSimulator<boost::mt19937, multitype_graph>
        (gen, graph, verbose);
    } else if (simType == "Rewire") {
      sim = new Simulators::RewireSimulator<boost::mt19937, multitype_graph>
        (gen, graph, verbose);
    } else if (simType == "NetEvo") {
      sim = new Simulators::NetEvoSimulator<boost::mt19937, multitype_graph>
        (gen, graph, verbose);
    } else {
      std::cerr << "ERROR: unknown simulator: " << simType << std::endl;
      return 1;
    }
  } else {
    std::cerr << "ERROR: no simulator specified" << std::endl;
    return 1;
  }

  // read additional simulation options
  po::options_description temp_options;
  temp_options.add(main_options).add(sim->getOptions());
  try {
    po::store(po::command_line_parser(argc, argv).options(temp_options).
              allow_unregistered().run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }
  po::notify(vm);
  sim->parse_options(vm);

  /******************************************************************/
  // read model options
  /******************************************************************/
  model = sim->getModel(); // model is specified in Simulator options
  if (model == 0) {
    std::cerr << "ERROR: no model" << std::endl;
    return 1;
  }

  temp_options.add(model->getOptions());

  try {
    po::store(po::command_line_parser(argc, argv).options(temp_options).
              allow_unregistered().run(), vm);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }

  // initialise model
  sim->InitModel(vm);

  std::vector<Label> edgeTypes = model->getEdgeTypes(); //!< available edge
                                                        //!types
  unsigned int nEdgeTypes = edgeTypes.size(); //!< number of edge types


  /******************************************************************/
  // read network
  /******************************************************************/
  po::options_description graph_options
    ("\nGraph options");

  graph_options.add_options()
    ("file,f", po::value<std::string>(),
     "file to read from (one file for all edge types)");

  if (nEdgeTypes > 1) {
    for (unsigned int i = 0; i < nEdgeTypes; ++i) {
      graph_options.add_options()
        ((edgeTypes[i].getText()+"-file").c_str(), po::value<std::string>(),
         ("file to read "+edgeTypes[i].getText()+"-network from").c_str());
    }
  }

  /******************************************************************/
  // determine initial conditions (read from file or randomise)
  /******************************************************************/
  po::options_description ic_options
    ("\nInitial conditions");

  ic_options.add_options()
    ("init,i", po::value<std::string>(),
     "file to get initial conditions from")
    ("init-lattice",
     "init file is png lattice rather than graphviz")
    ("same-ic",
     "start with the same initial conditions for each run")
    ("base,b", po::value<std::string>()->default_value
     (model->getVertexStates().begin()->getText()),
     "base state of individuals")
    ;

  for (std::vector<Label>::const_iterator it =
         model->getVertexStates().begin();
       it != model->getVertexStates().end(); it++) {
    ic_options.add_options()
      (it->getText().c_str(), po::value<unsigned int>()->default_value(0),
       ("number of randomly chosen " + it->getText()).c_str());
  }

  // read options from command line
  po::options_description all_options;
  all_options.add(temp_options).add(ic_options).add(graph_options);

  std::vector<std::string> unregistered;

  try {
    po::parsed_options parsed = po::command_line_parser(argc, argv).
      options(all_options).allow_unregistered().run();
    po::store(parsed, vm);
    unregistered =
      po::collect_unrecognized(parsed.options, po::exclude_positional);
  }
  catch (std::exception& e) {
    std::cerr << "ERROR parsing command line parameters: " << e.what()
              << std::endl;
    return 1;
  }

  po::notify(vm);
  for (std::vector<std::string>::iterator it = unregistered.begin();
       it != unregistered.end(); it++) {
    std::cerr << "WARNING: ignoring unknown option " << *it << std::endl;
  }

  // print long help with all options if requested
  if (vm.count("longhelp")) {
    std::cout << all_options << std::endl;
    return 0;
  }

  po::notify(vm);

  /******************************************************************/
  // collect filenames for the different edge types
  /******************************************************************/

  std::vector<std::string> fileNames; //!< The filenames of the .dot files for
                                      //!each edge type

  if (vm.count("file")) {
    // read all from one file
    for (unsigned int i = 0; i < edgeTypes.size(); i++) {
      fileNames.push_back
        (vm["file"].as<std::string>());
    }
  } else if (edgeTypes.size() > 1) {
    bool graphError = true;
    for (unsigned int i = 0; i < edgeTypes.size(); i++) {
      std::stringstream s;
      s.str("");
      s << edgeTypes[i].getText() << "-file";
      if (vm.count(s.str())) {
        fileNames.push_back(vm[s.str()].as<std::string>());
        graphError = false;
      } else {
        fileNames.push_back("");
        std::cerr << "WARNING: no " << edgeTypes[i].getText()
                  << "-graph file specified" << std::endl;
      }
      if (graphError) {
        std::cerr << "ERROR: no graphs" << std::endl;
        return 1;
      }
    }
  } else if (vm.count("init")) {
    for (unsigned int i = 0; i < edgeTypes.size(); i++) {
      fileNames.push_back
        (vm["init"].as<std::string>());
    }
  } else {
    std::cerr << "ERROR: no graph file specified" << std::endl;
    return 1;
  }

  /******************************************************************/
  // create network from the files for each edge type
  /******************************************************************/

  for (unsigned int i = 0; i < edgeTypes.size(); i++) {

    onetype_graph temp_graph;

    // reading graph structure and initial state from file
    if (fileNames[i].size() > 0) {
      int edgesRead = read_graph(temp_graph, fileNames[i], i);
      if (edgesRead > 0) {

        // update number of vertices
        if (verbose > 1) {
          std::cout << "Read " << edgesRead << " edges from graph file "
                    << fileNames[i] << std::endl;
        }
      } else if (edgesRead == 0) {
        std::cerr << "WARNING: no " << model->getEdgeTypes()[i] << "-edges to "
                  << "read from " << fileNames[i] << std::endl;
      } else {
        std::cerr << "ERROR: Could not read from graph file "
                  << fileNames[i] << std::endl;
        return 1;
      }
      if (num_vertices(temp_graph) > num_vertices(graph)) {
        boost::add_vertices(graph,
                            num_vertices(temp_graph) - num_vertices(graph));
        N = num_vertices(graph);
      }
    }

    // copy edges to main graph
    boost::graph_traits<multitype_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(temp_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, temp_graph), target(*ei, temp_graph),
               temp_graph[*ei], graph);
    }
  }

  // checking graph
  if (num_vertices(graph) == 0) {
    std::cerr << "ERROR: no vertices" << std::endl;
    std::cerr << main_options << graph_options << std::endl;
    return 1;
  }

  /******************************************************************/
  // mark parallel edges
  /******************************************************************/

  // mark parallel edges (to be used later for counting parallel edges)
  unsigned int parallel_edges = mark_parallel_edges(graph);
  if (verbose > 1 && nEdgeTypes > 1) {
    std::cout << "No. of parallel edges is: " << parallel_edges
              << std::endl;
  }
  setup_edge_index_map(graph);

  doIO = (sim->doIO());

  unsigned int extLength = 0;

  if (numSims > 1) {
    std::stringstream ext;
    ext << numSims;
    extLength = ext.str().length();
    std::cout << "Running " << numSims << " simulations" << std::endl;
  }

  /******************************************************************/
  // set initial conditions
  /******************************************************************/
  std::string baseString = vm["base"].as<std::string>();
  unsigned int baseState = 0;

  // assign baseState
  while (baseState < model->getVertexStates().size() &&
         (model->getVertexStates()[baseState].getText() != baseString)) {
    baseState++;
  }

  // check if baseState known
  if (baseState >= model->getVertexStates().size()) {
    std::cerr << "WARNING: unknown base state: " << baseString << std::endl;
    std::cerr << "Setting base state to "
              << model->getVertexStates()[0].getText() << std::endl;
    baseState = 0;
  }

  if (vm.count("init")) { // read initial state from file

    std::string icFileName = vm["init"].as<std::string>();
    int verticesRead;
    if (vm.count("init-lattice")) {
      // verticesRead = read_initial_lattice(graph, icFileName, *model);
      std::cerr << "ERROR: cannot read initial lattices at the moment "
                << "because of clang/pngwriter incompatibilities."
                << std::endl;
      return 1;
    } else {
      verticesRead = read_initial_graph(graph, icFileName, *model, baseState,
                                        verbose);
    }
    if (verticesRead < 0) {
      std::cerr << "ERROR: could not read from file " << icFileName
                << std::endl;
      return 1;
    } else if (static_cast<unsigned int>(verticesRead) <
               num_vertices(graph)) {
      std::cerr << "ERROR: found only " << verticesRead << " vertex "
                << "states in " << icFileName << std::endl;
      return 1;
    } else if (verbose > 1) {
      std::cout << "Read " << verticesRead << " initial states from "
                << icFileName << std::endl;
    }
    clearIC = false;
  }
  if (vm.count("same-ic")) {
    keepIC = true;
  }

  if (vm.count("output-dir")) {
    outputDir += "/"+vm["output-dir"].as<std::string>();
  }

  /******************************************************************/
  // read simulation paramters
  /******************************************************************/

  if (doIO) {
    // remove existing data
    if (fs::exists(outputDir)) {
      fs::directory_iterator iter(outputDir), end_iter;
      for (; iter != end_iter; ++iter) {
        if (iter->path().filename().string().substr(0,3) == "run") {
          if (fs::is_directory(*iter)) {
            fs::remove_all(outputDir+"/"+(iter->path().filename().string()));
          }
        }
      }
    } else {
      mkdir(outputDir.c_str(), 0755);
    }
    std::ofstream cmdLineFile;
    std::string cmdLineFileName = outputDir + "/cmd_line.txt";
    cmdLineFile.open(cmdLineFileName.c_str(), std::ios::out);
    cmdLineFile << argv[0] << " ";
    for (int i = 1; i < argc; i++) {
      cmdLineFile << argv[i] << " ";
    }
    cmdLineFile << std::endl;
    cmdLineFile.close();
  }

  /******************************************************************/
  // run the simulations
  /******************************************************************/

  for (unsigned int nSim = 1; nSim <= numSims; nSim++) {

    if (verbose >= 1) {
      std::cout << "----- " << "run #" << nSim << std::endl;
    } else if (numSims > 0 && !verbose) {
      std::cout << ".";
      std::cout.flush();
    }

    if (nSim > 1) {
      // recover graph
      graph.clear();
      copy_graph(saved_graph, graph);
    }

    /******************************************************************/
    // set initial vertex states
    /******************************************************************/

    init.clear();
    if (randomIC) {
      // how many random vertices of each state are to be initialised
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
    }
    // // set random vertices of baseState to zero
    // init[baseState] = 0;

    /******************************************************************/
    // generate vertices' state
    /******************************************************************/

    if (clearIC) {
      // set the initial state of all vertices to the base state
      boost::graph_traits<multitype_graph>::vertex_iterator vi, vi_end;
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        graph[*vi].setStatePtr(model->newState());
        graph[*vi].state->setState(baseState);
      }
    }

    // sum over vector init to make sure the sum is less than N
    unsigned int initSum = 0;
    for (std::vector<unsigned int>::iterator it = init.begin();
         it != init.end(); it++) {
      initSum += (*it);
    }
    if (initSum > N) {
      std::cerr << "WARNING: number of vertices to select randomly ("
                << initSum << ") higher than number of total vertices ("
                << N << ")" << std::endl;
    }

    std::vector <bool> touched(N, false);
    // initialise init[i] vertices of type i
    boost::graph_traits<multitype_graph>::vertex_descriptor v;
    for (unsigned int i = 0; i < model->getVertexStates().size(); i++) {
      for (unsigned int j = 0; j < init[i]; j++) {
        bool inserted = false;
        while (!inserted) {
          v = boost::random_vertex(graph, gen);
          if (touched[v] == false) {
            graph[v].state->setState(i);
            touched[v] = true;
            inserted = true;
          }
        }
        if (verbose >= 2) {
          std::cout << "Vertex #" << v << " is assigned state "
                    << model->getVertexStates()[i] << std::endl;
        }
      }
    }


    if (nSim == 1) {
      // save graph states
      copy_graph(graph, saved_graph);
    }

    if (keepIC) {
      clearIC = false;
      randomIC = false;
    }

    std::stringstream runStr(std::ios::in | std::ios::out | std::ios::ate);
    runStr << "run" << std::setfill('0') << std::setw(extLength) << nSim;

    std::string runOutputDir = outputDir + "/" + runStr.str();

    if (doIO) mkdir(runOutputDir.c_str(), 0755);

    sim->setDir(runOutputDir);

    /******************************************************************/
    // write model parameters
    /******************************************************************/

    std::ofstream paramFile;
    std::string paramFileName = runOutputDir + "/" + runStr.str() + ".prm";

    if (doIO) {
      try {
        paramFile.open(paramFileName.c_str(), std::ios::out);
      }
      catch (std::exception &e) {
        std::cerr << "... unable to open parameter file "
                  << paramFileName << std::endl;
        std::cerr << "... Standard exception: " << e.what() << std::endl;
        return 1;
      }

      for (unsigned int i = 0; i < edgeTypes.size(); i++) {
        paramFile << edgeTypes[i].getText() << "-graph: " << fileNames[i]
                  << std::endl;
      }
      paramFile << std::endl << *model << std::endl;

      try {
        paramFile.close();
      }
      catch (std::exception &e) {
        std::cerr << "... unable to close parameter file "
                  << paramFileName << std::endl;
        std::cerr << "... Standard exception: " << e.what() << std::endl;
        return 1;
      }
    }

    /******************************************************************/
    // initialise Simulator
    /******************************************************************/
    if (!sim->initialise()) {
      std::cerr << "ERROR: Could not intialise simulator" << std::endl;
      return 1;
    }

    /******************************************************************/
    // run simulation
    /******************************************************************/
    while ((sim->extraStopCondition() == false) && sim->updateState()) {
      sim->updateStats();
    }

    if (doIO) {
      if (sim->getTime() > 0) {
        std::ofstream gpFile;
        std::string gpFileName = runOutputDir+"/"+runStr.str()+".gp";

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
        gpFile << "Tmax=" << sim->getTime() << std::endl;
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
    }

  }

  // free memory
  delete sim;

  std::cout << std::endl;

  return 0;

}
