/*! \file sim_statistics.hh
  \brief Functions for simulation statistics.
*/

#ifndef SIM_STATISTICS_HH
#define SIM_STATISTICS_HH

#include <fstream>
#include <sys/stat.h>
#include <math.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>

#include "network/include/nearest_infected.hh"
#include "network/include/community_structure.hh"
#include "StatRecorder.hh"

namespace fs = boost::filesystem;

//! \addtogroup sim_statistics Simulation statistics

//----------------------------------------------------------
/*! \brief Count vertices.
  
Counts the number of vertices of a given state in a graph

\param[in] g The graph containing the vertices
\param[in] nVertexStates The number of vertex states.
\return A vector of state counts.
\ingroup sim_statistics
*/
template <typename Graph>
std::vector<unsigned int>
count_vertices(Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;

  std::vector<unsigned int> counts;

  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
    if (g[*vi].state->getState()+1 > counts.size()) {
      counts.resize(g[*vi].state->getState()+1, 0);
    }
    ++counts[g[*vi].state->getState()];
  }
   
  return counts;
}

//----------------------------------------------------------
/*! \brief Count effective vertices.
  
Counts the effecitve number of vertices (weighted by the detailed state) of a
given state in a graph 

\param[in] g The graph containing the vertices
\param[in] nVertexStates The number of vertex states.
\return A vector of state counts.
\ingroup sim_statistics
*/
template <typename Graph, typename Model>
std::vector<double>
count_effective_vertices(const Graph& g, unsigned int nVertexStates)
{
  typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;
  typedef typename Model::StateType state_type;

  std::vector<double> counts;

  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
    if (g[*vi].state->getState()+1 > counts.size()) {
      counts.resize(g[*vi].state->getState()+1, 0.);
    }
    state_type* s = dynamic_cast<state_type*>(g[*vi].state);
    counts[s->getState()] += (1-s->getInfo());
  }
   
  return counts;
}

//----------------------------------------------------------
/*! \brief Calculate average trait

Calculates the average trait in each state

\param[in] g The graph containing the vertices
\param[in] nVertexStates The number of vertex states.
\return A vector of state counts.
\ingroup sim_statistics
*/
template <typename Graph, typename Model>
std::vector<double>
calc_avg_trait(const Graph& g, unsigned int nVertexStates)
{
  typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;
  typedef typename Model::StateType state_type;

  std::vector<double> avgs(nVertexStates,0.);
  std::vector<unsigned int> counts(nVertexStates,0);

  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
    state_type* s = dynamic_cast<state_type*>(g[*vi].state);
    ++counts[s->getState()];
    avgs[s->getState()] += s->getTrait(0);
  }

  for (unsigned int i = 0; i < nVertexStates; ++i) {
    if (avgs[i] > 0) avgs[i] /= counts[i];
  }
   
  return avgs;
}

//----------------------------------------------------------
/*! \brief Count state pairs.
  
Counts the number of pairs of states in a graph.

\param[in] g The graph containing the vertices and edges
\param[in] nVertexStates The number of vertex states.
\param[in] nEdgeTypes The number of edge types.
\return A 3d array of pairs, the first index of which denounces edge type and
the other two vertex states
\ingroup sim_statistics
*/
template <typename Graph>
boost::multi_array<unsigned int, 3>
count_state_pairs(Graph& g, unsigned int nVertexStates,
                  unsigned int nEdgeTypes)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
  typedef boost::multi_array<unsigned int, 3> array_type;

  array_type counts(boost::extents
                    [nEdgeTypes]
                    [nVertexStates]
                    [nVertexStates]);
    
  // count all edges of type et, including parallel
  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    std::vector<unsigned int> states(2);
    states[0] = g[source(*ei, g)].state->getState();
    states[1] = g[target(*ei, g)].state->getState();
    std::sort(states.begin(), states.end());
    ++counts[g[*ei].type][states[0]][states[1]];
  }

  return counts;
}

//----------------------------------------------------------
/*! \brief Count parallel edges.
  
Loops over all edges and counts for which edges there is another edge
connecting the same two vertices. 

\param[in] g The graph containing the vertices and edges.
\param[in] nVertexStates The number of vertex states.
\return A 2d array of parallel pairs the indices denoting the two vertex states.
\ingroup sim_statistics
*/
template <typename Graph>
boost::multi_array<unsigned int, 2>
count_parallel_edges(Graph& g, unsigned int nVertexStates)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;

  // count only PARALLEL edges
  typedef boost::multi_array<unsigned int, 2> array_type;

  array_type counts(boost::extents
                    [nVertexStates]
                    [nVertexStates]);

  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    if ((g[*ei].type == 0) && g[*ei].parallel) {
      std::vector<unsigned int> states(2);
      states[0] = g[source(*ei, g)].state->getState();
      states[1] = g[target(*ei, g)].state->getState();
      std::sort(states.begin(), states.end());
      ++counts[states[0]][states[1]];
    }
  }

  // since we looped over all edges of one type, no need to divide by 2
  return counts;
}

//----------------------------------------------------------
/*! \brief Count state triples.
  
Counts the number of triples of states in the network.

\param[in] g The graph containing the vertices.
\param[in] nVertexStates The number of vertex states.
\param[in] nEdgeTypes The number of edge types.
\ingroup sim_statistics
*/
template <typename Graph>
boost::multi_array<unsigned int, 5>
count_state_triples(Graph& g, unsigned int nVertexStates,
                    unsigned int nEdgeTypes)
{
  typedef typename boost::graph_traits<Graph>::out_edge_iterator
    out_edge_iterator;
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
  typedef typename boost::graph_traits<Graph>::vertex_iterator
    vertex_iterator;

  typedef boost::multi_array<unsigned int, 5> array_type;
    
  array_type counts(boost::extents
                    [nEdgeTypes]
                    [nEdgeTypes]
                    [nVertexStates]
                    [nVertexStates]
                    [nVertexStates]);
    
  vertex_iterator vi, vi_end;
  out_edge_iterator oi, oi_end, oi2;

  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
    if (out_degree(*vi, g) > 0) {
      for (boost::tie(oi, oi_end) = out_edges(*vi, g); oi != oi_end; oi++) {
        for(oi2 = oi + 1; oi2 != oi_end; oi2++) {
          if (target(*oi, g) != target(*oi2, g)) {

            std::vector<unsigned int> edges(2);
            std::vector<unsigned int> states(2);
            edges[0] = g[*oi].type;
            edges[1] = g[*oi2].type;
            states[0] = g[target(*oi, g)].state->getState();
            states[1] = g[target(*oi2, g)].state->getState();
            std::sort(edges.begin(), edges.end());
            std::sort(states.begin(), states.end());
            
            ++counts
              [edges[0]]
              [edges[1]]
              [g[*vi].state->getState()]
              [states[0]]
              [states[1]];

          }
        }
      }
    }
  }

  return counts;
    
}

//----------------------------------------------------------
/*! \brief Print simulation status.
  
Prints the number of vertices in all possible states and, if desired, the
number of edges of all types connecting vertices of all possible states, as
well as the number of parallel edges connecting such vertices.

\param[in] g The graph containing the vertices and edges
\param[in] m The model to be used to find all possible states
\param[in] pairs Whether pairs should be counted as well
\param[in] triples Whether triples should be counted as well
\ingroup sim_statistics
*/
template <typename Graph, typename Model>
class print_sim_status
  : public Funct<Graph>
{

public:

  print_sim_status(const Model& m, bool pairs = false, bool triples = false)
    : m(m), pairs(pairs), triples(triples)
  {;}

  void doit(const Graph& g, std::string dir, double time, unsigned int count)
  {
    unsigned int nEdgeTypes = m.getEdgeTypes().size();

    std::cout << "Time elapsed: " << time << std::endl;
    std::cout << std::endl;
    std::cout << "Vertex count: " << std::endl;
    
    // count states
    std::vector<unsigned int> vertexCount = count_vertices(g);
    unsigned int nVertexStates = vertexCount.size();
    for (unsigned int i=0; i < nVertexStates; i++) {
      if (vertexCount[i] > 0) {
        std::cout << m.getVertexState(i) << ": " << vertexCount[i]
                  << std::endl;
      }
    }
    
    std::cout << std::endl;
    if (pairs) {
      
      std::cout << "Pair count: " << std::endl;
      
      // count state pairs
      boost::multi_array<unsigned int, 3> pairCount =
        count_state_pairs(g, nVertexStates, nEdgeTypes);
      for (unsigned int i=0; i < nEdgeTypes; i++) {
        std::cout << m.getEdgeTypes()[i] << "-type: " << std::endl;
        for (unsigned int j=0; j < nVertexStates; j++) {
          for (unsigned int k=j; k < nVertexStates; k++) {
            if (pairCount[i][j][k] > 0) {
              std::cout << m.getVertexState(j) << "-" << m.getVertexState(k)
                        << ": " << pairCount[i][j][k] << std::endl;
            }
          }
        }
      }

      if (m.getEdgeTypes().size() > 1) {
        // count parallel pairs
        std::cout << "parallel:" << std::endl;
        
        boost::multi_array<unsigned int, 2> parallelCount =
          count_parallel_edges(g, nVertexStates);
        
        for (unsigned int j=0; j < nVertexStates; j++) {
          for (unsigned int k=j; k < nVertexStates; k++) {
            if (parallelCount[j][k] > 0) {
              std::cout << m.getVertexState(j) << m.getVertexState(k)
                        << ": " << parallelCount[j][k] << std::endl;
            }
          }
        }
        std::cout << std::endl;
      }
    }
    
    if (triples) {
      std::cout << "Triple count: " << std::endl;
      
      boost::multi_array<unsigned int, 5> tripleCount =
        count_state_triples(g, nVertexStates, nEdgeTypes);
      
      for (unsigned int i=0; i < nEdgeTypes; i++) {
        for (unsigned int j=i; j < nEdgeTypes; j++) {
          std::cout << m.getEdgeTypes()[i] << m.getEdgeTypes()[j] << "-triples: "
                    << std::endl;
          for (unsigned int k=0; k < nVertexStates; k++) {
            for (unsigned int l=0; l < nVertexStates; l++) {
              for (unsigned int n=l; n < nVertexStates; n++) {
                std::cout << m.getVertexState(l) << m.getVertexState(k)
                          << m.getVertexState(n) << ": "
                          << tripleCount[i][j][k][l][n] << std::endl;
              }
            }
          }
        }
      }
      std::cout << std::endl;
    }
  }

private: 

  const Model& m;
  bool pairs;
  bool triples;
  
};

template <typename Graph, typename Simulator>
class write_epi_stats
  : public Funct<Graph>
{

public:

  write_epi_stats(const Simulator& s, bool v)
    : sim(s), verbose(v)
  {;}

  void doit(const Graph& g, std::string dir, double time, unsigned int count)
  {
    if (time > 0.) {
      std::ofstream statsFile;
      std::string statsFileName = dir+"/stats.dat";
      try {
        statsFile.open(statsFileName.c_str(), std::ios::out);
      }
      catch (std::exception &e) {
        std::cerr << "... unable to open stats file "
                  << statsFileName << " for writing" << std::endl;
        std::cerr << "... Standard exception: " << e.what() << std::endl;
      }
      statsFile << "Cumulative number of infections: "
                << sim.getNumInfections() << std::endl;
      statsFile << "Cumulative number of informations: " << sim.getNumInformations()
                << std::endl;
      try {
        statsFile.close();
      }
      catch (std::exception &e) {
        std::cerr << "... unable to close stats file "
                  << statsFileName << std::endl;
        std::cerr << "... Standard exception: " << e.what() << std::endl;
      }
      
      if (verbose) {
        std::cout << "Cumulative number of infections: " << sim.getNumInfections()
                  << std::endl;
        std::cout << "Cumulative number of informations: " << sim.getNumInformations()
                  << std::endl;
      }
    }
  }

private: 

  const Simulator& sim;
  bool verbose;
  
};

template <typename Graph, typename Simulator>
class write_r0
  : public Funct<Graph>
{

public:

  write_r0(const std::vector<std::vector<unsigned int> >& gen,
           const std::vector<unsigned int>& genInf, bool v)
    : generations(gen), genInf(genInf), verbose(v)
  {;}

  void doit(const Graph& g, std::string dir, double time, unsigned int count)
  {
    if (time > 0.) {
      std::ofstream r0File;
      std::string r0FileName = dir+"/r0.dat";
      try {
        r0File.open(r0FileName.c_str(), std::ios::out);
      }
      catch (std::exception &e) {
        std::cerr << "... unable to open r0 file "
                  << r0FileName << " for writing" << std::endl;
        std::cerr << "... Standard exception: " << e.what() << std::endl;
      }
      for (unsigned int i = 0; i < genInf.size(); ++i) {
        if (generations[i].size() > 0) {
          r0File << i << " "
                 << genInf[i]/static_cast<double>(generations[i].size())
                 << std::endl;
          if (verbose) {
            std::cout << "R0 generation " << i << " "
                      << genInf[i]/static_cast<double>(generations[i].size())
                      << std::endl;
          }
        }
      }
      try {
        r0File.close();
      }
      catch (std::exception &e) {
        std::cerr << "... unable to close r0 file "
                  << r0FileName << std::endl;
        std::cerr << "... Standard exception: " << e.what() << std::endl;
      }
      
      
    }
  }

private: 

  const std::vector<std::vector<unsigned int> >& generations;
  const std::vector<unsigned int>& genInf;
  bool verbose;
  
};


//----------------------------------------------------------
/*! \brief Write simulations data statistics to a file
  
Writes a line containing the current time (passed as a parameter) to a file,
as well as the number of vertices in all possible states and, if desired, the 
number of edges of all types connecting vertices of all possible states, as
well as the number of parallel edges connecting such vertices.

\param[in] g The graph containing the vertices and edges
\param[in] m The model to be used to find all possible states
\param[in] t The current time to be written to the first column
\param[in] ofile A stream to the output file
\param[in] pairs Whether pairs should be counted as well (intensive on
computation) or not 
\param[in] triples Whether triples should be counted as well (intensive on
computation) or not 
\ingroup sim_statistics
*/
template <typename Graph, typename Model>
class write_sim_data:
  public Funct<Graph>
{
  
public:

  write_sim_data(const Model& m, bool pairs = false, bool triples = false)
    : m(m), pairs(pairs), triples(triples)
  {;}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    std::string outputFileName = dir + "/" + "singles.sim.dat";
    std::ofstream outputFile;
    
    try {
      outputFile.open(outputFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::ate);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
      std::cerr << "Will not write simulation counts to file." << std::endl;
    }
  
    unsigned int nEdgeTypes = m.getEdgeTypes().size();
    
    std::stringstream line("");

    // first in line is current time
    outputFile << time << '\t';
    
    // count vertices' states
    std::vector<unsigned int> vertexCount = count_vertices(g);
    unsigned int nVertexStates = vertexCount.size();
    for (unsigned int i = 0; i < nVertexStates; i++) {
      line << vertexCount[i] << '\t';
    }
    line << std::endl;
    outputFile << line.str();
    outputFile.close();
    line.str("");

    if (pairs) {
      outputFileName = dir + "/" + "pairs.sim.dat";
      
      try {
        outputFile.open(outputFileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::ate);
      }
      catch (std::exception &e) {
        std::cerr << "Unable to open output file: " << e.what() << std::endl;
        std::cerr << "Will not write simulation counts to file." << std::endl;
      }
      
      outputFile << time << '\t';

      boost::multi_array<unsigned int, 3> pairCount =
        count_state_pairs(g, nVertexStates, nEdgeTypes);
      // count pairs
      for (unsigned int i = 0; i < nEdgeTypes; i++) {
        for (unsigned int j = 0; j < nVertexStates; j++) {
          for (unsigned int k = j; k < nVertexStates; k++) {
            line << pairCount[i][j][k] << '\t';
          }
        }
      }

      // count parallel pairs
      if (nEdgeTypes > 1) {
        boost::multi_array<unsigned int, 2> parallelCount =
          count_parallel_edges(g, nVertexStates);
      
        for (unsigned int j = 0; j < nVertexStates; j++) {
          for (unsigned int k = j; k < nVertexStates; k++) {
            line << parallelCount[j][k] << '\t';
          }
        }
      }
      line << std::endl;
      outputFile << line.str();
      outputFile.close();
      line.str("");
    }

    if (triples) {
      outputFileName = dir + "/" + "triples.sim.dat";
      
      try {
        outputFile.open(outputFileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::ate);
      }
      catch (std::exception &e) {
        std::cerr << "Unable to open output file: " << e.what() << std::endl;
        std::cerr << "Will not write simulation counts to file." << std::endl;
      }
      
      outputFile << time << '\t';


      boost::multi_array<unsigned int, 5> tripleCount =
        count_state_triples(g, nVertexStates, nEdgeTypes);
      
      for (unsigned int i=0; i < nEdgeTypes; i++) {
        for (unsigned int j=i; j < nEdgeTypes; j++) {
          for (unsigned int k=0; k < nVertexStates; k++) {
            for (unsigned int l=0; l < nVertexStates; l++) {
              for (unsigned int n=l; n < nVertexStates; n++) {
                line << tripleCount[i][j][k][l][n] << '\t';
              }
            }
          }
        }
      }

      line << std::endl;
      outputFile << line.str();
      outputFile.close();
      line.str("");
    }
    
  }

private:
  
  const Model& m;
  bool pairs;
  bool triples;
  
};


template <typename Graph, typename Model>
class write_effective_data:
  public Funct<Graph>
{
public:
  write_effective_data(const Model& m)
    : m(m)
  {;}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    std::string outputFileName = dir + "/" + "effective.sim.dat";
    std::ofstream outputFile;
    
    try {
      outputFile.open(outputFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::ate);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
      std::cerr << "Will not write simulation counts to file." << std::endl;
    }
    
    outputFile << time << '\t';

    unsigned int nVertexStates = m.getVertexStates().size();
    std::vector<double> effVertexCount =
      count_effective_vertices<Graph, Model>(g, nVertexStates);
    for (unsigned int i = 0; i < nVertexStates; i++) {
      outputFile << effVertexCount[i] << '\t';
    }
    outputFile << std::endl;
    outputFile.close();
  }    
  
private:

  const Model& m;

};

template <typename Graph, typename Model>
class write_avg_trait:
  public Funct<Graph>
{
public:
  write_avg_trait(const Model& m)
    : m(m)
  {;}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    std::string outputFileName = dir + "/" + "avg_trait.sim.dat";
    std::ofstream outputFile;
    
    try {
      outputFile.open(outputFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::ate);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
      std::cerr << "Will not write simulation counts to file." << std::endl;
    }
    
    outputFile << time << '\t';

    unsigned int nVertexStates = m.getVertexStates().size();
    std::vector<double> avgTrait =
      calc_avg_trait<Graph, Model>(g, nVertexStates);
    for (unsigned int i = 0; i < nVertexStates; i++) {
      outputFile << avgTrait[i] << '\t';
    }
    outputFile << std::endl;
    outputFile.close();
  }    
  
private:

  const Model& m;

};

//----------------------------------------------------------
/*! \brief Detail states in pair participants
  
Calculate the average detailed state in on part of a given pair of states.

\param[in] g The graph containing the vertices and edges
\param[in] edgeType The edge type to consider
\param[in] state1 State of the first neighbour
\param[in] state2 State of the second neighbours (to be averaged)
\return A vector containg the detail values of the state2 neighbours
\ingroup sim_statistics
*/
template <typename Graph, typename Model>
std::vector<double> pair_detail_states(const Graph& g, unsigned int edgeType,
                                       unsigned int state1, unsigned int state2)
{
  typedef typename boost::graph_traits<Graph>::edge_iterator
    edge_iterator;
  typedef typename Model::StateType state_type;

  std::vector<double> detail_states;
    
  // count all edges of type et, including parallel
  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
    if (g[*ei].type == edgeType) {
      if (g[source(*ei, g)].state->getState() == state1 &&
          g[target(*ei, g)].state->getState() == state2) {
        state_type* s = dynamic_cast<state_type*>(g[target(*ei, g)].state);
        detail_states.push_back(s->getInfo());
      } else if (g[source(*ei, g)].state->getState() == state2 &&
                 g[target(*ei, g)].state->getState() == state1) {
        state_type* s = dynamic_cast<state_type*>(g[source(*ei, g)].state);
        detail_states.push_back(s->getInfo());
      }
    }
  }

  return detail_states;
}

template <typename Graph, typename Model>
class write_risk_info:
  public Funct<Graph>
{

  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {

    std::string fileName = generateFileName(dir + "/risk/risk", count);
    std::ofstream riskFile;
    if (!fs::exists(dir+"/risk")) {
      try {
        mkdir((dir+"/risk").c_str(), 0755);
      } 
      catch (std::exception &e) {
        std::cerr << "... unable to create directory "
                  << dir << "/risk" << std::endl;
        return;
      }
    }
    try {
      riskFile.open(fileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "... unable to open risk-info output file " 
                << fileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }
    
    std::vector<double> detail_states =
      pair_detail_states<Graph, Model>(g, 0, 1, 0);
    
    for (std::vector<double>::iterator vi = detail_states.begin();
         vi != detail_states.end(); ++vi) {
      riskFile << *vi << std::endl;
    }

    try {
      riskFile.close();
    }
    catch (std::exception &e) {
      std::cerr << "... unable to close risk-info output file " 
                << fileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }
    
  }
};

template <typename Graph, typename Model>
class write_detail_dist:
  public Funct<Graph>
{

  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    typedef typename Model::StateType state_type;
    
    std::string fileName = generateFileName(dir + "/dist/dist", count);
    std::ofstream distFile;
    if (!fs::exists(dir+"/dist")) {
      try {
        mkdir((dir+"/dist").c_str(), 0755);
      } 
      catch (std::exception &e) {
        std::cerr << "... unable to create directory "
                  << dir << "/dist" << std::endl;
        return;
      }
    }
    try {
      distFile.open(fileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "... unable to open output file " 
                << fileName << " for writing the information distribution"
                << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }
    
    vertex_iterator vi, vi_end;
    state_type* s = dynamic_cast<state_type*>(g[*vi].state);
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      distFile << s->getInfo() << std::endl;
    }

    try {
      distFile.close();
    }
    catch (std::exception &e) {
      std::cerr << "... unable to close output file " 
                << fileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }
  }
};

template <typename Graph, typename Model>
class write_info_dis_corr:
  public Funct<Graph>
{
public:
  write_info_dis_corr(const Model& m, unsigned int et)
    : m(m), et(et)
  {;}

  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    typedef typename boost::edge_property_type<Graph>::type::value_type
      edge_property_type;
    typedef typename Model::StateType state_type;

    std::string fileName = generateFileName(dir + "/corr/corr", count);
    std::ofstream corrFile;
    if (!fs::exists(dir+"/corr")) {
      try {
        mkdir((dir+"/corr").c_str(), 0755);
      } 
      catch (std::exception &e) {
        std::cerr << "... unable to create directory "
                  << dir << "/corr" << std::endl;
        return;
      }
    }
    try {
      corrFile.open(fileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "... unable to open correlation output file " 
                << fileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }
      
    std::vector<unsigned int> infected_distances =
      boost::nearest_infected<Graph, edge_property_type, Model>
      (g, edge_property_type(et), m);
    unsigned int vertex_counter = 0;
      
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      state_type* s = dynamic_cast<state_type*>(g[*vi].state);
      corrFile << s->getInfo() << '\t'
               << infected_distances[vertex_counter] << std::endl;
      ++vertex_counter;
    }

    try {
      corrFile.close();
    }
    catch (std::exception &e) {
      std::cerr << "... unable to close correlation output file " 
                << fileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }
      
  }
    
private:
  const Model& m;
  unsigned int et;
};

template <typename Graph, typename Model>
class write_sim_graph
  : public Funct<Graph>
{
public:
  write_sim_graph(const Model& m)
    : m(m)
  {;}

  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    if (!fs::exists(dir+"/images")) {
      try {
        mkdir((dir+"/images").c_str(), 0755);
      } 
      catch (std::exception &e) {
        std::cerr << "... unable to create directory "
                  << dir << "/images" << std::endl;
        return;
      }
    }
    std::string fileName = generateFileName(dir + "/images/frame", count, 
                                            "graph");
    write_graph(g, fileName, m,time);
  }
    
private:
  const Model& m;
};

template <typename Graph, typename Model>
class write_sim_lattice
  : public Funct<Graph>
{
public:
  write_sim_lattice(const Model& m)
    : m(m)
  {;}

  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    if (!fs::exists(dir+"/images")) {
      try {
        mkdir((dir+"/images").c_str(), 0755);
      } 
      catch (std::exception &e) {
        std::cerr << "... unable to create directory "
                  << dir << "/images" << std::endl;
        return;
      }
    }
    std::string fileName = generateFileName(dir + "/images/frame", count);
    write_png(g, fileName, m, time);
  }
    
private:
  const Model& m;
};

template <typename Graph>
class write_comp_dist
  : public Funct<Graph>
{
public:

  write_comp_dist(const Graph& g) : largestComponent(boost::num_vertices(g)) {;}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    if (!fs::exists(dir+"/CompDist")) {
      try {
        mkdir((dir+"/CompDist").c_str(), 0755);
      } 
      catch (std::exception &e) {
        std::cerr << "... unable to create directory "
                  << dir << "/CompDist" << std::endl;
        return;
      }
    }
    std::string fileName = generateFileName(dir + "/CompDist/cd", count);
    write_component_dist(g, fileName);
  }

  unsigned int getLargest() const { return largestComponent; }

private:

  unsigned int largestComponent;

};

template <typename Graph, typename Model>
class write_community_structure
  : public Funct<Graph>
{
public:

  write_community_structure(Graph& g, const Model& m,
                            bool graphs = false,
                            bool modularity = false,
                            unsigned int verbose = 0) :
    m(m), graphs(graphs), modularity(modularity), verbose(verbose)
  {}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    Graph temp_graph = g;
    double mod = boost::community_structure(g, temp_graph, 1., verbose,
                                            graphs, modularity);
    std::cout << "done" << std::endl;

    if (graphs) {
      if (!fs::exists(dir+"/comm")) {
        try {
          mkdir((dir+"/comm").c_str(), 0755);
        } 
        catch (std::exception &e) {
          std::cerr << "... unable to create directory "
                    << dir << "/comm" << std::endl;
          return;
        }
      }
      std::string fileName = generateFileName(dir + "/comm/comm", count);
      write_graph(temp_graph, fileName, m, time);
    }
    if (modularity) {
      std::string outputFileName = dir + "/" + "modularity.sim.dat";
      std::ofstream outputFile;
      
      try {
        outputFile.open(outputFileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::ate);
      }
      catch (std::exception &e) {
        std::cerr << "Unable to open output file: " << e.what() << std::endl;
        std::cerr << "Will not write modularity to file." << std::endl;
      }
      
      outputFile << time << '\t';
      outputFile << mod;
      outputFile << std::endl;
      outputFile.close();
    }
  }

private:
  const Model& m;
  bool graphs;
  bool modularity;
  unsigned int verbose;
};

template <typename Graph>
class write_same_state_components
  : public Funct<Graph>
{
public:

  write_same_state_components(Graph& g, unsigned int verbose = 0) :
    verbose(verbose)
  {}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    Graph temp_graph = g;

    std::vector<int> component(num_vertices(g));
    int num = boost::connected_components(g, &component[0]);
    if (verbose >=2 ) {
      std::cout << num << " components:" << std::endl;;
    }
    
    double fraction;
    unsigned int nonZeroVertices = 0;
    // loop over components and calculate state fractions
    for (int i = 0; i < num; ++i) {
      std::vector<unsigned int> stateDist;
      unsigned int componentSize = 0;
      for (unsigned int j = 0; j < num_vertices(g); ++j) {
        if (component[j] == i) {
          if (g[j].state->getState() > 0) {
            if (g[j].state->getState() + 1 > stateDist.size()) {
              stateDist.resize(g[j].state->getState() + 1, 0);
            }
            ++stateDist[g[j].state->getState()];
            ++nonZeroVertices;
            ++componentSize;
          }
        }
      }
      unsigned int compContr = 0;
      if (componentSize > 1) {
        for (unsigned int j = 0; j < stateDist.size(); ++j) {
          compContr +=  stateDist[j]*(stateDist[j] - 1);
        }
        fraction += compContr / static_cast<double>(componentSize - 1);
      }
    }
    if (nonZeroVertices > 0) {
      fraction /= static_cast<double>(nonZeroVertices);
    } else {
      fraction = 0;
    }
    std::string outputFileName = dir + "/" + "same_state.sim.dat";
    std::ofstream outputFile;
    
    try {
        outputFile.open(outputFileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::ate);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
      std::cerr << "Will not write same state fraction to file." << std::endl;
    }
    
    outputFile << time << '\t';
    outputFile << fraction;
    outputFile << std::endl;
    outputFile.close();
  }

private:
  unsigned int verbose;
};

template <typename Graph>
class write_components
  : public Funct<Graph>
{
public:

  write_components(Graph& g, unsigned int verbose = 0) :
    verbose(verbose)
  {}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    Graph temp_graph = g;

    std::vector<int> component(num_vertices(g));
    int num = boost::connected_components(g, &component[0]);
    if (verbose >=2 ) {
      std::cout << num << " components:" << std::endl;;
    }
    std::vector<std::vector<int> > components(num);
    for (int i = 0; i < num; ++i) {
      for (unsigned int j = 0; j < num_vertices(g); ++j) {
        if (component[j] == i) {
	  components[i].push_back(j);
        }
      }
    }
    if (!fs::exists(dir+"/comp")) {
        try {
          mkdir((dir+"/comp").c_str(), 0755);
        } 
        catch (std::exception &e) {
          std::cerr << "... unable to create directory "
                    << dir << "/comm" << std::endl;
          return;
        }
    }
    std::string fileName = generateFileName(dir + "/comp/comp", count);

    std::ofstream outputFile;
    
    try {
        outputFile.open(fileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::ate);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
      std::cerr << "Will not write same state fraction to file." << std::endl;
    }
    
    for (int i = 0; i < num; ++i) {
      for (unsigned int j = 0; j < components[i].size(); ++j) {
	outputFile << i << '\t';
	outputFile << components[i][j];
	outputFile << std::endl;
      }
    }
    outputFile.close();
  }

private:
  unsigned int verbose;
};

template <typename Graph, typename Model>
class write_group_modularity
  : public Funct<Graph>
{
public:

  write_group_modularity(Graph& g, const Model& m, unsigned int verbose = 0,
			 unsigned int edgeType = 0) :
    m(m), verbose(verbose), edgeType(edgeType)
  {}
  
  virtual void doit(const Graph& g, std::string dir, double time,
                    unsigned int count)
  {
    std::string outputFileName = dir + "/" + "groupmod.sim.dat";
    std::ofstream outputFile;
    
    try {
      outputFile.open(outputFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::ate);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open output file: " << e.what() << std::endl;
      std::cerr << "Will not write simulation counts to file." << std::endl;
    }
  
    unsigned int nEdgeTypes = m.getEdgeTypes().size();
    unsigned int nVertexStates = m.getVertexStates().size();
    
    std::stringstream line("");

    // first in line is current time
    outputFile << time << '\t';

    double modularity = 0.;

    boost::multi_array<unsigned int, 3> pairCount =
      count_state_pairs(g, nVertexStates, nEdgeTypes);
    // count pairs
    for (unsigned int j = 0; j < nVertexStates; j++) {
      unsigned int intraCount = 0;
      unsigned int interCount = 0;
      for (unsigned int k = 0; k < nVertexStates; k++) {
	if (j == k) {
	  intraCount += pairCount[edgeType][j][k];
	} else {
	  interCount += pairCount[edgeType][j][k];
	}
      }
      modularity += (intraCount / static_cast<double>(num_edges(g)))
	- pow((2*intraCount + interCount) / (2*num_edges(g)), 2);
    }
    
    outputFile << modularity << std::endl;
    outputFile.close();
  }

private:
  const Model& m;
  unsigned int verbose;
  unsigned int edgeType;
};


//----------------------------------------------------------

#endif
