/*! \file simulate.cc
  \brief The main simulation program.
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
#include "network/include/Edge.hh"

#include "sim_statistics.hh"

#include "GillespieSimulator.hh"
#include "ChrisSimulator.hh"
#include "Vertex.hh"

// models
#include "InfoSIRS.hh"
#include "DimInfoSIRS.hh"
#include "VaccinationSIRS.hh"
#include "ProtectiveSIRS.hh"
#include "SingleSIRS.hh"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> multitype_graph;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                              Vertex, Edge> onetype_graph;

std::string generateFileName(std::string nameBase, unsigned int id)
{
  std::stringstream s;
  s << nameBase << std::setw(6) << std::setfill('0') << id;
  return s.str();
}

int main(int argc, char* argv[])
{
  Model* model = 0;
  Simulator* sim = 0;

  std::vector <unsigned int> init;
  
  /******************************************************************/
  // read parameters
  /******************************************************************/
  unsigned int N = 0;
  double stopTime = 0;
  unsigned int stopInfections;
  unsigned int stopInformations;

  double outputData = 0.;
  double outputGraphviz = -1.;
  double outputDist = -1.;

  char* dataDir = getenv("DATADIR");
  std::string outputDir;
  if (dataDir) {
    mkdir(dataDir, 0755);
    outputDir = dataDir;
  } else {
    outputDir = "output";
  }

  std::string outputFileName = "";
  std::ofstream* outputFile = 0;

  unsigned int verbose = 0;
  bool printStats = false;

  unsigned int numSims = 1;

  bool generateIC = true; // default is to generate initial conditions
  bool keepIC = false; // default is not to keep initial conditions fro each run
  
  bool allFromOne;

  bool pairs = false;
  bool triples = false;

  po::options_description command_line_options
    ("\nUsage: simulate [options]... \n\nMain options");

  command_line_options.add_options()
    ("help,h",
     "produce help message")
    ("longhelp,H",
     "produce long help message including all options")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("params-file,p",po::value<std::string>(),
     "file containing model parameters")
    ;

  po::options_description sim_options
    ("\nSimulation options");

  sim_options.add_options()
    ("sim", po::value<std::string>()->default_value("Gillespie"),
     "simulator to use (Gillespie, Chris)")
    ("usemodel,m", po::value<std::string>()->default_value("InfoSIRS"),
     "model to use (InfoSIRS, DimInfoSIRS, VaccinationSIRS, ProtectiveSIRS, SingleSIRS)")
    ("tmax", po::value<double>()->default_value(stopTime),
     "time after which to stop\n(use tmax=0 to run until extinction)")
    ("imax", po::value<unsigned int>()->default_value(0),
     "number of infections after which to stop (if >0)")
    ("pmax", po::value<unsigned int>()->default_value(0),
     "number of informations after which to stop (if >0)")
    ("data,d", po::value<double>()->default_value(outputData),
     "write output data at arg timesteps (0 for no data output)")
    ("graphviz,g", po::value<double>()->default_value(outputGraphviz),
     "create graphviz output in the images directory at arg timesteps")
    ("info-dist,i", po::value<double>()->default_value(outputDist),
     "create information distribution in the dist directory at arg timesteps")
    ("lattice,l", 
     "paint output as pixelised lattices")
    ("output-dir",po::value<std::string>()->default_value(outputDir),
     "output directory for graph and data output")
    ("nsims", po::value<unsigned int>()->default_value(1),
     "number of simulation runs to produce (on a given graph)")
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
  if (vm.count("help")) {
    std::cout << command_line_options << sim_options << std::endl;
    return 0;
  }

  if (vm.count("usemodel")) {
    std::string modelString = vm["usemodel"].as<std::string>();
    if (modelString == "InfoSIRS") {
      model = new Models::InfoSIRS(verbose);
    } else if (modelString == "DimInfoSIRS") {
      model = new Models::DimInfoSIRS(verbose); 
    } else if (modelString == "VaccinationSIRS") {
      model = new Models::VaccinationSIRS(verbose); 
    } else if (modelString == "ProtectiveSIRS") {
      model = new Models::ProtectiveSIRS(verbose); 
    } else if (modelString == "SingleSIRS") {
      model = new Models::SingleSIRS(verbose); 
    } else {
      std::cerr << "Error: unknown model: " << modelString << std::endl;
      return 1;
    }
  } else {
    std::cerr << "Error: no model specified" << std::endl;
    return 1;
  }

  std::vector<Label> edgeTypes = model->getEdgeTypes();
  unsigned int nEdgeTypes = edgeTypes.size();

  po::options_description graph_options
    ("\nGraph options");

  graph_options.add_options()
    ("file,f", po::value<std::string>(),
     "file to read from (one file for all edge types)");

  for (unsigned int i = 0; i < nEdgeTypes; ++i) {
    graph_options.add_options()
      ((edgeTypes[i].getText()+"-file").c_str(), po::value<std::string>(),
       ("file to read "+edgeTypes[i].getText()+"-network from").c_str());
  }
  
  po::options_description statistics_options
    ("\nStatistics");
  
  statistics_options.add_options()
    ("print-stats",
     "print stats at the end of run")
    ("pairs",
     "count pairs")
    ("triples",
     "count triples")
    ;
  
  po::options_description ic_options
    ("\nInitial conditions");

  ic_options.add_options()
    ("init,i", po::value<std::string>(),
     "graphviz file to get initial conditions from")
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
  
  po::options_description model_options = model->getOptions();

  // read options from command line
  po::options_description all_options;
  all_options.add(command_line_options).add(sim_options).
    add(ic_options).add(graph_options).add(statistics_options).
    add(model_options);
                                   
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
  // initialise random generator
  /******************************************************************/

  unsigned int seed;
  struct timeval tv;
  gettimeofday(&tv, 0);
  seed = tv.tv_sec + tv.tv_usec;
  boost::mt19937 gen(seed);

  /******************************************************************/
  // create graph variable
  /******************************************************************/
  multitype_graph graph, saved_graph;
  void (*graph_function)(const multitype_graph&, std::string,
                         const Model&, double) = 0;
            
  /******************************************************************/
  // read graph from file
  /******************************************************************/

  std::vector<std::string> fileNames(nEdgeTypes);
  
  if (vm.count("file")) {
    // read all from one file
    allFromOne = true;
    for (unsigned int i = 0; i < edgeTypes.size(); i++) {
      fileNames[i] = vm["file"].as<std::string>();
    }
  } else {
    allFromOne = false;
    for (unsigned int i = 0; i < edgeTypes.size(); i++) {
      std::stringstream s;
      s.str("");
      s << edgeTypes[i].getText() << "-file";
      if (vm.count(s.str())) {
        fileNames[i] = vm[s.str()].as<std::string>();
      } else {
        fileNames[i] = "";
      }
    }
  }
    
  for (unsigned int i = 0; i < edgeTypes.size(); i++) {
      
    onetype_graph temp_graph;
      
    // reading graph structure and initial state from file
    if (fileNames[i].size() > 0) {
      int edgesRead = read_graph(temp_graph, fileNames[i], i);
      if (edgesRead > 0) {
        
        // update number of vertices
        if (verbose) {
          std::cout << "Read " << edgesRead << " edges from graph file "
                    << fileNames[i] << std::endl;
        }
        if (num_vertices(temp_graph) > num_vertices(graph)) {
          boost::add_vertices(graph,
                              num_vertices(temp_graph) - num_vertices(graph));
	  N = num_vertices(graph);
        };
      } else if (edgesRead == 0) {
        std::cerr << "WARNING: no " << model->getEdgeTypes()[i] << "-edges to "
                  << "read from " << fileNames[i] << std::endl;
      } else {
        std::cerr << "ERROR: Could not read from graph file "
                  << fileNames[i] << std::endl;
        return 1;
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
    std::cerr << command_line_options << graph_options << std::endl;
    return 1;
  }

  // consider pairs
  if (vm.count("pairs")) {
    pairs = true;
  }
  // consider triples
  if (vm.count("triples")) {
    triples = true;
  }
  // print stats at end of run
  if (vm.count("print-stats")) {
    printStats = true;
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

  // timesteps after which to write graphviz output
  if (vm.count("graphviz")) {
    outputGraphviz = vm["graphviz"].as<double>();
  }

  // timesteps after which to write info distribution
  if (vm.count("info-dist")) {
    outputDist = vm["info-dist"].as<double>();
  }

  // paint images as lattice
  if (vm.count("lattice")) {
    graph_function = &boost::write_png;
  } else {
    graph_function = &boost::write_graph;
  }

  numSims = vm["nsims"].as<unsigned int>();
  
  unsigned int extLength = 0;
  
  if (numSims > 1) {
    std::stringstream ext;
    ext << numSims;
    extLength = ext.str().length();
    std::cout << "Running " << numSims << " simulations" << std::endl;
  }

  std::string baseString = vm["base"].as<std::string>();
  unsigned int baseState = 0;
  
  // assign baseState
  while (baseState < model->getVertexStates().size() &&
         (model->getVertexStates()[baseState].getText() != baseString)) {
    baseState++;
  }
  
  // check if baseState known
  if (baseState == model->getVertexStates().size()) {
    std::cerr << "WARNING: unknown base state: " << baseString << std::endl;
    std::cerr << "Setting base state to "
              << model->getVertexStates()[0].getText() << std::endl;
    baseState = 0;
  }              

  if (vm.count("init")) { // read initial state from file
    
    std::string icFileName = vm["init"].as<std::string>();
    int verticesRead = read_initial_graph(graph, icFileName, *model);
    if (verticesRead < 0) {
      std::cerr << "ERROR: could not read from file " << icFileName
                << std::endl;
      return 1;
    } else if (static_cast<unsigned int>(verticesRead) < num_vertices(graph)) {
      std::cerr << "ERROR: found only " << verticesRead << " vertex "
                << "states in " << icFileName << std::endl;
      return 1;
    } else if (verbose) {
      std::cout << "Read " << verticesRead << " initial states from "
                << icFileName << std::endl;
    }
    if (verbose >= 1) { 
      print_sim_status(graph, *model, pairs, triples);
    }
    generateIC = false;
  } else if (vm.count("same-ic")) {
    keepIC = true;
  }
    
  if (numSims > 0) {

    /******************************************************************/
    // initialise model
    /******************************************************************/
    if (numSims > 0) {
      model->Init(vm);
      if (verbose >=1) model->Print();
    }
    
    /******************************************************************/
    // create simulator
    /******************************************************************/
    
    if (vm.count("sim")) {
      std::string simType = vm["sim"].as<std::string>();
      if (simType == "Gillespie") {
        sim = new Simulators::GillespieSimulator<boost::mt19937, multitype_graph>
          (gen, graph, *model, verbose);
      } else if (simType == "Chris") {
        sim = new Simulators::ChrisSimulator<boost::mt19937, multitype_graph>
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
    
    stopTime = vm["tmax"].as<double>();
    stopInfections = vm["imax"].as<unsigned int>();
    stopInformations = vm["pmax"].as<unsigned int>();
    outputData = vm["data"].as<double>();
        
    // output directory
    if (vm.count("output-dir")) {
      outputDir = vm["output-dir"].as<std::string>();
      // remove existing data
      if (fs::exists(outputDir)) {
        fs::directory_iterator iter(outputDir), end_iter;
        for (; iter != end_iter; ++iter) {
          if (iter->leaf().substr(0,3) == "run") {
            if (fs::is_directory(*iter)) {
              fs::remove_all(outputDir+"/"+(iter->leaf()));
            }
          }
        }
      } else {
        mkdir(outputDir.c_str(), 0755);
      }
    }
  }
  
  for (unsigned int nSim = 1; nSim <= numSims; nSim++) {
    
    if (printStats) {
      std::cout << "----- " << "run #" << nSim << std::endl;
    } else if (numSims > 0) {
      std::cout << ".";
      std::cout.flush();
    }
    
    /******************************************************************/
    // generate new initial state
    /******************************************************************/
   
    if (generateIC) { // generate new initial state
      
      /******************************************************************/
      // set initial vertex states
      /******************************************************************/

      init.clear();
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
      
      // set random vertices of baseState to zero
      init[baseState] = 0;
      
      /******************************************************************/
      // generate vertices' state
      /******************************************************************/
      
      // set the initial state of all vertices to the base state
      boost::graph_traits<multitype_graph>::vertex_iterator vi, vi_end;      
      for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        graph[*vi].state =
          State(baseState, model->getInitDetail(baseState));
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
      
      // initialise init[i] vertices of type i
      boost::graph_traits<multitype_graph>::vertex_descriptor v;
      for (unsigned int i=0; i<model->getVertexStates().size(); i++) {
        for (unsigned int j=0; j<init[i]; j++) {
          bool inserted = false;
          while (!inserted) {
            v = boost::random_vertex(graph, gen);
            if (graph[v].state.base == baseState) {
              graph[v].state = State(i, model->getInitDetail(i));
              inserted = true;
            }
          }
          if (verbose >= 2) {
            std::cout << "Vertex #" << v << " is assigned state " 
                      << model->getVertexStates()[i] << std::endl;
          }
        }
      }
    } else {
      keepIC = true;
    }
    
    if (keepIC) {
      if (nSim == 1) {
        // save graph states
        saved_graph = graph;
        generateIC = false;
      } else {
        // recover graph states
        graph = saved_graph;
      }
    }
    
    std::stringstream runStr(std::ios::in | std::ios::out | std::ios::ate);
    runStr << "run" << std::setfill('0') << std::setw(extLength) << nSim;

    std::string runOutputDir = outputDir + "/" + runStr.str();
    
    mkdir(runOutputDir.c_str(), 0755);
    

    std::string runDataDir = runOutputDir;
    std::string runGraphDir = runOutputDir + "/images";
    std::string runDistDir = runOutputDir + "/dist";
    
    /******************************************************************/
    // open output file
    /******************************************************************/
    if (outputData > 0) {
      // create data directory
      mkdir(runDataDir.c_str(), 0755);
      outputFileName = runDataDir+"/"+runStr.str()+".sim.dat";
      
      try {
        outputFile = new std::ofstream();
        outputFile->open(outputFileName.c_str(), std::ios::out);
      }
      catch (std::exception &e) {
        std::cerr << "Unable to open output file: " << e.what() << std::endl;
      }
    }
    
    /******************************************************************/
    // initialise Simulator
    /******************************************************************/
    sim->initialise();
    
    // print time
    if (verbose)
      std::cout << "time elapsed: " << sim->getTime() << std::endl;
    
    // GraphViz output
    if (outputGraphviz >= 0) {
      // create graph directory
      mkdir(runGraphDir.c_str(), 0755);
      graph_function(graph, (runGraphDir + "/frame000000"), *model, -1);
    }
    // Information distribution output
    if (outputDist >= 0) {
      // create distribution directory
      mkdir(runDistDir.c_str(), 0755);
      write_detail_dist(graph, (runDistDir + "/dist000000"));
    }
    
    
    // prints data to outputFile
    std::string lastLine = "";
    if (outputFile) {
      lastLine = write_sim_data(graph, *model, sim->getTime(), *outputFile,
                                  pairs, triples);
    }
    if (verbose) print_sim_status(graph, *model, pairs, triples);
    
    /******************************************************************/
    // run simulation
    /******************************************************************/
    double nextDataStep = outputData;
    double nextGraphStep = outputGraphviz;
    double nextDistStep = outputDist;
    
    unsigned int steps = 0;
    unsigned int graphOutputNum = 1;
    unsigned int distOutputNum = 1;
    
    while ((stopTime >= 0 || sim->getTime()<stopTime) &&
           (stopInfections == 0 || (sim->getNumInfections() < stopInfections &&
                                    sim->getNumInfections()+1 > sim->getNumRecoveries())) &&
           (stopInformations == 0 || (sim->getNumInformations() < stopInformations && 
                                      sim->getNumInformations()+1 > sim->getNumForgettings())) &&
            sim->updateState()) {
      
      if (verbose >= 2) {
        print_sim_status(graph, *model, pairs, triples);
      }
      
      if (verbose && steps%100 == 0) {
        std::cout << "time elapsed: " << sim->getTime() << std::endl;
      }
      
      if (outputFile && sim->getTime() > nextDataStep) {
        lastLine = 
          write_sim_data(graph, *model, sim->getTime(), *outputFile, pairs, triples);
        if (outputData > 0) {
          do {
            nextDataStep += outputData;
          } while (sim->getTime() > nextDataStep);
        }
      }
      if ((outputGraphviz > 0) && (sim->getTime() > nextGraphStep)) {
        graph_function(graph,
                       generateFileName((runGraphDir +"/frame"),
                                        graphOutputNum),
                       *model,
                       sim->getTime());
        do {
	  nextGraphStep += outputGraphviz;
	} while (sim->getTime() > nextGraphStep);
        ++graphOutputNum;
      }
      if ((outputDist > 0) && (sim->getTime() > nextDistStep)) {
        write_detail_dist(graph,
                          generateFileName((runDistDir +"/dist"),
                                           distOutputNum));

        do {
	  nextDistStep += outputDist;
	} while (sim->getTime() > nextDistStep);
        ++distOutputNum;
      }
      ++steps;
    }
    
    if (verbose) std::cout << "Final status (" << sim->getTime() << "): " 
                           << std::endl;
    if (outputGraphviz >= 0) {
      graph_function(graph,
                     generateFileName((runGraphDir+"/frame"),
                                      graphOutputNum),
                     *model,
                     sim->getTime());
    }
    if (outputDist >= 0) {
      write_detail_dist(graph,
                     generateFileName((runDistDir+"/dist"),
                                      distOutputNum));
    }
    
    if (outputFile) {
      if (sim->getTime() < stopTime) {
        lastLine = 
          write_sim_data(graph, *model, sim->getTime(), *outputFile, pairs, triples);
      }
      *outputFile << stopTime << '\t' << lastLine;
      outputFile->close();
      delete outputFile;
      std::ofstream statsFile;
      std::string statsFileName = runDataDir+"/"+runStr.str()+".stats";
      statsFile.open(statsFileName.c_str(), std::ios::out);
      statsFile << "Cumulative number of infections: "
                << sim->getNumInfections() << std::endl;
      statsFile << "Cumulative number of informations: " << sim->getNumInformations() 
                << std::endl;
      boost::graph_traits<multitype_graph>::vertex_iterator vi, vi_end;
      statsFile.close();
    }
    
    if (verbose || printStats) {
      std::cout << "Cumulative number of infections: " << sim->getNumInfections() 
                << std::endl;
      std::cout << "Cumulative number of informations: " << sim->getNumInformations() 
                << std::endl;
      boost::graph_traits<multitype_graph>::vertex_iterator vi, vi_end;
    }
    
    if (verbose) print_sim_status(graph, *model, pairs, triples);
    
  }

  if (vm.count("write-file")) {

    if (stopTime == 0) stopTime = sim->getTime();
    if (stopTime > 0) {
      std::ofstream gpFile;
      std::string gpFileName = outputDir+"/"+"init.gp";
      
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
  }

  // free memory
  delete sim;
  delete model;
  
  return 0;
  
}
