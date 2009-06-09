/*! \file graph_structure.hh
  \brief Functions for the creation and modification of graphs.
*/

#ifndef GRAPH_STRUCTURE_HH
#define GRAPH_STRUCTURE_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <vector>
#include <iostream>

#include "graph_statistics.hh"
#include "Edge.hh"

//! \addtogroup graph_structure Graph structure
//! \addtogroup helper_functions Helper functions
//! \addtogroup graph_generators Graph generators

//----------------------------------------------------------
namespace boost {

  /*! \brief Add vertices to a graph.
    
  Adds a given number of (defaultly constructed) vertices to the graph.
  
  \param[in, out] g The graph to add the vertices to.
  \param[in] n The number of vertices to add.
  \ingroup graph_structure
  */
  template <typename Graph>
  void add_vertices(Graph& g, 
                    typename graph_traits<Graph>::vertices_size_type n)
  {
    for (unsigned int i = 0; i < n; i++) {
      add_vertex(g);
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Add vertices to a graph.
    
  Adds a given number of vertices with given properties to the graph.
  
  \param[in, out] g The graph to add the vertices to.
  \param[in] n The number of vertices to add.
  \param[in] v Template for the properties of the vertices to add.
  \ingroup graph_structure
  */
  template <typename Graph>
  void add_vertices(Graph& g, 
                    typename graph_traits<Graph>::vertices_size_type n,
                    const typename vertex_property_type<Graph>::type v)
  {
    for (unsigned int i = 0; i < n; i++) {
      add_vertex(v, g);
    }
  }
   
  //----------------------------------------------------------
  /*! \brief Add edges to a graph.
    
  Loops over a given iterator and add corresponding (defaultly constructed)
  edges to the graph.
  
  \param[in, out] g The graph to add the vertices to.
  \param[in] ei The begin edge iterator.
  \param[in] ei_end The end edge iterator.
  \ingroup graph_structure
  */
  template <typename Graph, typename EdgeIterator>
  void add_edge_structure(Graph& g, EdgeIterator ei, EdgeIterator ei_end)
  {
    while (ei != ei_end) {
      add_edge((*ei).first, (*ei).second, g);
      ++ei;
    }
  }

  //----------------------------------------------------------
  /*! \brief Add edges to a graph.
    
  Loops over a given iterator and add corresponding edges with given properties
  to the graph.
  
  \param[in, out] g The graph to add the vertices to.
  \param[in] ei The begin edge iterator.
  \param[in] ei_end The end edge iterator.
  \param[in] e Template for the properties of the edges to add.
  \ingroup graph_structure
  */
  template <typename Graph, typename EdgeIterator>
  void add_edge_structure(Graph& g, EdgeIterator ei, EdgeIterator ei_end,
                          const typename edge_property_type<Graph>::type e)
  {
    while (ei != ei_end) {
      add_edge((*ei).first, (*ei).second, e, g);
      ++ei;
    }
  }

  //----------------------------------------------------------
  /*! \brief Check if there is an unseen edge among stubs.
    
  Loops over all pairs of stubs and checks if they are contained in
  seen_edges. Needed for random_regular_graph and rewireEdges.
  
  \param[in] g The graph to check for unseen edges.
  \param[in] stubs The vector of stubs in the graph (basically a vector of
  vertices with each vertices occuring multiple times according to its degree).
  \param[in] seen_edges The edge pairs which are not usable anymore.
  \return true if there are still unseen vertices which can be connected with
  an edge, false if ther are not.
  \ingroup helper_functions
  */    

  template <typename Graph>
  bool suitable(Graph& g,
                std::vector<typename graph_traits<Graph>::vertex_descriptor>& 
		stubs,
                std::vector<std::vector<bool> >& seen_edges)
  {
    typedef typename std::vector
      <typename graph_traits<Graph>::vertex_descriptor>::iterator iter;
    
    for (iter s = stubs.begin(); s != stubs.end(); s++) 
      for (iter t = s + 1; t != stubs.end(); t++) 
        if (!seen_edges[*s][*t] && *s != *t) return 1; // success
    
    return 0; //failure
  }
  
  //----------------------------------------------------------
  /*! \brief Copy all edges from one graph to another
    
  Loops over all edges in a graph and copy them to antoher graph.
  
  \param[in] source_graph The graph to copy the edges from.
  \param[out] target_graph The graph to copy the egges to.
  \ingroup graph_structure
  */
  template <typename Graph1, typename Graph2>
  void copy_edges(Graph1& source_graph, Graph2& target_graph)
  {
    typedef typename graph_traits<Graph1>::edge_iterator
      edge_iterator;
    
    edge_iterator ei, ei_end;
    // loop over all edges in source graph and add to target_graph with
    // property et
    for (tie(ei, ei_end) = edges(source_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, source_graph), target(*ei, source_graph),
               source_graph[*ei], target_graph);
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Copy all edges from one graph to another, assigning them a given
    edge property
    
  Loops over all edges in a graph and copy them to another graph, assigning them
  a given edge property
  
  \param[in] source_graph The graph to copy the edges from.
  \param[out] target_graph The graph to copy the egges to.
  \param[in] et Property to assign to the newly created edges in
  target_graph.
  \ingroup graph_structure
  */
  template <typename Graph1, typename Graph2>
  void copy_edges(Graph1& source_graph, Graph2& target_graph,
                 typename edge_property_type<Graph2>::type et)

  {
    typedef typename graph_traits<Graph1>::edge_iterator
      edge_iterator;
    
    edge_iterator ei, ei_end;
    // loop over all edges in source graph and add to target_graph with
    // property et
    for (tie(ei, ei_end) = edges(source_graph); ei != ei_end; ei++) {
      add_edge(source(*ei, source_graph), target(*ei, source_graph),
               et, target_graph);
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Copy all edges with a given property from one graph to another.
    
  Loops over all edges in a graph and copy them to another graph if they possess
  a given edge property
  
  \param[in] source_graph The graph to copy the edges from.
  \param[out] target_graph The graph to copy the egges to.
  \param[in] et Property to assign to the newly created edges in
  target_graph.
  \ingroup graph_structure
  */
  template <typename Graph1, typename Graph2>
  void copy_edges_type(Graph1& source_graph, Graph2& target_graph,
                       typename edge_property_type<Graph1>::type::value_type et)

  {
    typedef typename graph_traits<Graph1>::edge_iterator
      edge_iterator;
    
    edge_iterator ei, ei_end;
    // loop over all edges in source graph and add to target_graph if 
    // they have property et
    for (tie(ei, ei_end) = edges(source_graph); ei != ei_end; ei++) {
      if (source_graph[*ei] == et) {
        add_edge(source(*ei, source_graph), target(*ei, source_graph),
                 et, target_graph);
      }
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Create a random regular graph.
    
  Creates a random regular graph by creating a vertex of stubs (a list of
  vertices with each vertex multiply according to degree) and connecting them
  randomly, as described in Kim and Vu's 2003 paper.
  
  \param[in] g The graph to fill
  \param[out] rrg_edges Vector of vertex pairs representing edges.
  \param[in] d The (constant) degree.
  \param[in] et The edge type.
  \param[in] r The random generator to be used.
  \param[in] countThreshold How many iterations to try
  \return true if generation of random graph successful, false if one is left
  without possible pairs of nodes before convergence to a full regular graph.
  \ingroup graph_generators
  */
  template <typename Graph, typename RandomGenerator>
  int random_regular_graph(Graph& g,
			   typename edge_property_type<Graph>::type et,
			   unsigned int d,
			   RandomGenerator& r,
			   unsigned int countThreshold = 100)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;

    typedef std::vector<std::pair<unsigned int, unsigned int> > GraphEdges;
    typedef typename GraphEdges::iterator GraphEdgeIter;

    GraphEdges rrg_edges;
    
    uniform_real<> uni_dist(0, 1);
    variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);

    bool success = false;
    unsigned int count = 0;

    while (!success && count < countThreshold) {

      ++count;
      // define seen_edges NxN zero matrix
      std::vector<std::vector<bool> >
	seen_edges(num_vertices(g), std::vector<bool>(num_vertices(g), false));
      
      // loop over graph and fill seen_edges
      edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
	seen_edges[source(*ei, g)][target(*ei, g)] = true;
	seen_edges[target(*ei, g)][source(*ei, g)] = true;
      }
      
      // define stubs
      std::vector<vertex_descriptor> stubs;
      
      // init stubs
      for (unsigned int i = 0; i < num_vertices(g); ++i)
	for (unsigned int j = 0; j < d; ++j)
	  stubs.push_back(i);
      
      // construct random regular graph
      bool stillSuitable = true;

      while (stubs.size() && stillSuitable) {
	
	// select source and target from stubs
	unsigned int src = 
	  static_cast<vertex_descriptor>(uni_gen() * stubs.size());
	unsigned int trg = 
	  static_cast<vertex_descriptor>(uni_gen() * stubs.size());
	unsigned int source = stubs[src];
	unsigned int target = stubs[trg];
	
	// check if pair is suitable
	if ((source != target) && !(seen_edges[source][target])) {
	  
	  // removing source and target from stubs
	  stubs.erase(stubs.begin() + src);
	  if (src > trg) {
	    stubs.erase(stubs.begin() + trg);
	  } else {
	    stubs.erase(stubs.begin() + trg - 1);
	  }
	  
	  // update seen_edges
	  seen_edges[source][target] = true;
	  seen_edges[target][source] = true;         
	  
	  // add edge to rrg_edges
	  rrg_edges.push_back(std::make_pair(source, target));
	  
	} else { // check if suitable stubs left
	  if (!suitable(g, stubs, seen_edges)) {
	  // failure - no more suitable pairs
	    stillSuitable = false; 
	  }
	}
      }   
      if (stillSuitable) {
	success = true;
      }
    }

    for (GraphEdgeIter it = rrg_edges.begin(); it != rrg_edges.end(); it++) {
      boost::add_edge((*it).first, (*it).second, et, g);
    }
    
    if (success) {
      return count; // success
    } else {
      return -1; // failure
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Create a graph following a given degree distribution
    
  Creates a graph with a degree distribution given from a file by
  creating a vertex of stubs (a list of vertices with each vertex
  multiply according to degree) and connecting them 
  randomly, as described in Kim and Vu's 2003 paper.
  
  \param[in] fileName The file to read the degree distribution from.
  \param[in] g The graph to mark the parallel edges in.
  \param[in] et The edge type.
  \param[in] d The (constant) degree.
  \param[in] r The random generator to be used.
  \param[in] countThreshold How many iterations to try
  \return true if generation of random graph successful, false if one is left
  without possible pairs of nodes before convergence to a full regular graph.
  \ingroup graph_generators
  */
  
  template <typename Graph, typename RandomGenerator>
  int configuration_graph(std::string fileName, Graph& g,
			  typename edge_property_type<Graph>::type et,
			  RandomGenerator& r,
			  unsigned int countThreshold = 100)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;
    
    typedef std::vector<std::pair<unsigned int, unsigned int> > GraphEdges;
    typedef typename GraphEdges::iterator GraphEdgeIter;
    
    GraphEdges graph_edges;
    
    uniform_real<> uni_dist(0, 1);
    variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);

    std::vector<double> degreeDist;

    std::ifstream degreeFile;
    try {
      degreeFile.open(fileName.c_str(), std::ios::in);
    } catch (std::exception &e) {
      std::cerr << e.what() << std::endl;
      return -1;
    }
    std::string line;
    double norm = 0.;

    if (degreeFile.is_open()) {
      getline(degreeFile, line);
      while(!degreeFile.eof()) {
	//read line
	if (line.size() > 0) {
	  std::istringstream linestream(line);
	  unsigned int degree = 0;
	  double prob = 0.;
	  linestream >> degree >> prob;
	  if (degree + 1 > degreeDist.size()) {
	    degreeDist.resize(degree + 1, 0.);
	    degreeDist[degree] = prob;
	    norm += prob;
	  }
	}
	getline(degreeFile, line);
      }
    } else {
      std::cerr << "ERROR: problem opening " << fileName << std::endl;
      return -1;
    }

    for (unsigned int i = 0; i < degreeDist.size(); ++i) {
      degreeDist[i] /= norm;
    }
    
    bool success = false;
    unsigned int count = 0;

    while (!success && count < countThreshold) {

      ++count;
      // define seen_edges NxN zero matrix
      std::vector<std::vector<bool> >
	seen_edges(num_vertices(g), std::vector<bool>(num_vertices(g), false));
      
      // loop over graph and fill seen_edges
      edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
	seen_edges[source(*ei, g)][target(*ei, g)] = true;
	seen_edges[target(*ei, g)][source(*ei, g)] = true;
      }
      
      // define stubs
      std::vector<vertex_descriptor> stubs;
      
      // init stubs
      unsigned int vertexCount = 0;
      for (unsigned int i = 0; i < degreeDist.size(); ++i) {
	for (unsigned int j = 0; 
	     j < round(degreeDist[i] * num_vertices(g)); ++j) {
	  for (unsigned int k = 0; k < i; ++k) {
	    stubs.push_back(vertexCount);
	  }
	  ++vertexCount;
	}
      }
      std::cout << "TEST: " << num_vertices(g) << " " << (vertexCount-1) << 
	std::endl;
      
      // construct random regular graph
      bool stillSuitable = true;

      while (stubs.size() && stillSuitable) {
	
	// select source and target from stubs
	unsigned int src = 
	  static_cast<vertex_descriptor>(uni_gen() * stubs.size());
	unsigned int trg = 
	  static_cast<vertex_descriptor>(uni_gen() * stubs.size());
	unsigned int source = stubs[src];
	unsigned int target = stubs[trg];
	
	// check if pair is suitable
	if ((source != target) && !(seen_edges[source][target])) {
	  
	  // removing source and target from stubs
	  stubs.erase(stubs.begin() + src);
	  if (src > trg) {
	    stubs.erase(stubs.begin() + trg);
	  } else {
	    stubs.erase(stubs.begin() + trg - 1);
	  }
	  
	  // update seen_edges
	  seen_edges[source][target] = true;
	  seen_edges[target][source] = true;         
	  
	  // add edge to rrg_edges
	  graph_edges.push_back(std::make_pair(source, target));
	  
	} else { // check if suitable stubs left
	  if (!suitable(g, stubs, seen_edges)) {
	  // failure - no more suitable pairs
	    stillSuitable = false; 
	  }
	}
      }   
      if (stillSuitable) {
	success = true;
      }
    }

    for (GraphEdgeIter it = graph_edges.begin(); 
	 it != graph_edges.end(); it++) {
      boost::add_edge((*it).first, (*it).second, et, g);
    }
    
    if (success) {
      return count; // success
    } else {
      return -1; // failure
    }
  }
  
  //----------------------------------------------------------
  /*! \brief Randomly rewire edges in a graph.
    
  Removes a given fraction of edges from a graph and randomly rewires them,
  preserving the degree distribution.
  
  \param[in, out] g The graph to rewiwre the edges in.
  \param[in] r The random generator to use for selecting the edges to rewire.
  \param[in] rewireFraction The fraction of edges to rewire.
  \param[in] verbose Whether to be verbose.
  \return 0 if successful, -1 if rewiring fails because one is left without
  connectable edges according to the given degree distribution before rewiring
  all the nodes.
  \ingroup graph_structure
  */
  template <typename Graph, typename RandomGenerator>
  int rewireEdges(Graph& g, RandomGenerator& r, double rewireFraction, 
                  unsigned int verbose)
  {
    typedef typename graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;
    typedef typename edge_property_type<Graph>::type
      edge_property_type;
    
    std::vector< std::vector<bool> >
      seen_edges(num_vertices(g), std::vector<bool>(num_vertices(g), false));
    std::vector<vertex_descriptor> stubs;

    edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(g);
    // initialize edge type from the first edge in the graph
    unsigned int et = g[*ei].type;
    // loop over all edges and initialise seen_edges to prevent self-loops
    for (; ei != ei_end; ei++) {
      seen_edges[source(*ei, g)][target(*ei, g)] = true;
      seen_edges[target(*ei, g)][source(*ei, g)] = true;
    }
    
    unsigned int N = num_edges(g);
    
    if (rewireFraction > 0. && rewireFraction <= 1.) {
      uniform_real<> uni_dist(0, 1);
      variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);
      // calculate number of links to rewire
      unsigned int num_rewire =
        static_cast<unsigned int>(rewireFraction * N);
      if (verbose >=1) {
        std::cout << "Randomly rewiring " << num_rewire << " edges" << std::endl;
      }
      
      // remove num_rewire random edges
      for (unsigned int i = 0; i < num_rewire; i++) {
        edge_descriptor e = random_edge(g, r);
        stubs.push_back(source(e, g));
        stubs.push_back(target(e, g));
        if (verbose >=2) {
          std::cout << "Removing edge " << source(e, g) << "--" 
	            << target(e, g) << std::endl;
        }
        boost::remove_edge(e, g);
      }
      // rewire the removed edges
      while (num_rewire > 0) {
        // randomly select two stubs
        unsigned int src =
          static_cast<unsigned int>(uni_gen() * stubs.size());
        unsigned int trg = 
          static_cast<unsigned int>(uni_gen() * stubs.size());
        // check if suitable pair
        if ((stubs[src] != stubs[trg]) && (!seen_edges[stubs[src]][stubs[trg]])) {
          if (verbose >=2) {
            std::cout << "Adding edge " << stubs[src] << "--" 
  	              << stubs[trg] << std::endl;
          }
          // if suitable, add the edge to the graph
          boost::add_edge(stubs[src], stubs[trg], edge_property_type(et), g);
          // remove stubs from vector
          stubs.erase(stubs.begin() + src);
          if (src > trg) {
            stubs.erase(stubs.begin() + trg);
          } else {
            stubs.erase(stubs.begin() + trg - 1);
          }
          --num_rewire;
          if (verbose >=2) {
            std::cout << num_rewire << " to go" << std::endl;
          }
        } else {
          if (!suitable(g, stubs, seen_edges)) {
            return -1; // failure - no more suitable pairs
          }
        }
      }
    } else {
      std::cerr << "ERROR: rewire fraction must be between 0 and 1" << std::endl;
      std::cerr << "no rewiring performed" << std::endl;
    }
    return 0;
  }
  
  //----------------------------------------------------------
  /*! \brief Randomly rewire edges in a graph, clustering them in another graph.
    
  Removes a given fraction of edges from a graph and randomly rewires them,
  preserving the degree distribution and clustering them in another graph.
  
  \param[out] cg The graph to cluster.
  \param[in, out] g The graph to rewire the edges in.
  \param[in] r The random generator to use for selecting the edges to rewire.
  \param[in] rewireFraction The fraction of edges to rewire.
  \param[in] verbose Whether to be verbose.
  \return 0 if successful, -1 if rewiring fails because one is left without
  connectable edges according to the given degree distribution before rewiring
  all the nodes.
  \ingroup graph_structure
  */
  template <typename ClusterGraph, typename Graph, typename RandomGenerator>
  int rewireClustered(const ClusterGraph& cg, Graph& g, RandomGenerator& r,
                      double rewireFraction, unsigned int verbose)
  {
    typedef typename graph_traits<Graph>::edge_iterator
      g_edge_iterator;
    typedef typename graph_traits<Graph>::edge_descriptor
      g_edge_descriptor;
    typedef typename graph_traits<Graph>::out_edge_iterator
      g_out_edge_iterator;
    typedef typename edge_property_type<Graph>::type
      g_edge_property_type;
    typedef typename graph_traits<ClusterGraph>::out_edge_iterator
      cg_out_edge_iterator;
    typedef typename graph_traits<ClusterGraph>::vertex_descriptor
      cg_vertex_descriptor;

    // put all edges into a vector for speed reasons...

    g_edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(g);

    // initialize edge type from the first edge in graph g (which is to be
    // rewired)
    unsigned int et = g[*ei].type;
    std::vector<std::pair<unsigned int, unsigned int> > temp_edges;

    // loop over all edges in graph g and put them into temp_edges if they do
    // not exist
    
    for(; ei != ei_end; ei++) {
      cg_out_edge_iterator oi, oi_end;
      bool found = false;
      for (tie(oi, oi_end) = out_edges(source(*ei, g), cg);
           (oi != oi_end) && !found; oi++) {
        if (edge(target(*oi, cg), target(*ei, g), cg).second) found = true;
      }
      if (!found) {
        temp_edges.push_back(std::make_pair(source(*ei, g), target(*ei, g)));
        std::cout << *ei << std::endl;
      }
    }

    if (verbose >=2) {
      std::cout << temp_edges.size() << " edges can be rewired" << std::endl;
    }

    if (rewireFraction > 0. && rewireFraction <= 1.) {
      uniform_real<> uni_dist(0, 1);
      variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);
      // calculate number of links to rewire
      unsigned int num_rewire =
        static_cast<unsigned int>(rewireFraction * num_edges(g));
      if (verbose >=1) {
        std::cout << "Randomly rewiring " << num_rewire << " edges, clustering "
                  << "them in the remaining graph." << std::endl;
      }
      
      // rewire edges
      while (num_rewire > 0 && temp_edges.size() > 1) {
        std::cout << "Got " << temp_edges.size() << " left." << std::endl;
        // randomly select edge
        bool reverse = false;
        bool validEdge = false;
        unsigned int rand_edge;
        while (!validEdge) {
          rand_edge =
            static_cast<unsigned int>(uni_gen() * temp_edges.size() * 2);
          if (rand_edge > temp_edges.size()-1) {
            reverse = true;
            rand_edge -= temp_edges.size();
          }
          if (edge(temp_edges[rand_edge].first,
                   temp_edges[rand_edge].second, g).second) {
            validEdge = true;
          } else {
            temp_edges.erase(temp_edges.begin() + rand_edge);
          }
        }
        
        g_edge_descriptor randomEdge = 
	  edge(temp_edges[rand_edge].first,
               temp_edges[rand_edge].second, g).first;
        cg_vertex_descriptor source_vertex;
        cg_vertex_descriptor target_vertex;
        if (reverse) {
          source_vertex = target(randomEdge, g);
          target_vertex = source(randomEdge, g);
        } else {
          source_vertex = source(randomEdge, g);
          target_vertex = target(randomEdge, g);
        }
        // select random neighbour of source vertex in cluster graph
        if (out_degree(source_vertex, cg) > 0) {
          unsigned int rand_no =
            static_cast<unsigned int>
            (uni_gen() * out_degree(source_vertex, cg));
          cg_out_edge_iterator cg_oi, cg_oi_end;
          unsigned int rand_idx = 0;
          tie(cg_oi, cg_oi_end) = out_edges(source_vertex, cg);
          while (rand_idx < rand_no) {
            cg_oi++;
            ++rand_idx;
          }
          cg_vertex_descriptor rand_temp = target(*cg_oi, cg);
          if (rand_temp != target_vertex) {
            // select random neighbour of that vertex in cluster graph
            rand_no =
              static_cast<unsigned int>(uni_gen() * out_degree(rand_temp, cg));
            rand_idx = 0;
            tie(cg_oi, cg_oi_end) = out_edges(rand_temp, cg);
            while (rand_idx < rand_no) {
              cg_oi++;
              ++rand_idx;
            }
            cg_vertex_descriptor rand_nb = target(*cg_oi, cg);
            
            // randomly select second edge
            if (out_degree(rand_nb, g) > 0 && rand_nb != source_vertex) {
              rand_no =
                static_cast<unsigned int>(uni_gen() * out_degree(rand_nb, g));
              g_out_edge_iterator g_oi, g_oi_end;
              rand_idx = 0;
              tie(g_oi, g_oi_end) = out_edges(rand_nb, g);
              while (rand_idx < rand_no) {
                g_oi++;
                ++rand_idx;
              }
              cg_vertex_descriptor rand_new = target(*g_oi, g);
              
              // rewire if both edges to be created do not exist yet
              if (rand_new != source_vertex && rand_new != target_vertex &&
                  edge(source_vertex, rand_nb, g).second == false &&
                  edge(target_vertex, rand_new, g).second == false) {
                if (verbose >= 2) {
                  std::cout << "Rewiring " << randomEdge << " and " << *g_oi;
                }
                boost::remove_edge(randomEdge, g);
                boost::remove_edge(*g_oi, g);
                temp_edges.erase(temp_edges.begin()+rand_edge);
                bool found = false;
                for(std::vector<std::pair<unsigned int, unsigned int> >::iterator it =
                      temp_edges.begin(); (it != temp_edges.end()) && (!found); it++) {
                  if ((it->first == rand_nb && it->second == rand_new) ||
                      (it->first == rand_new && it->second == rand_nb)) {
                    found = true;
                    *it = std::make_pair(target_vertex, rand_new);
                  }
                }
                
                g_edge_descriptor new_edge1 =
                  boost::add_edge(source_vertex, rand_nb,
                                  g_edge_property_type(et), g).first;
                g_edge_descriptor new_edge2 =
                  boost::add_edge(target_vertex, rand_new,
                                  g_edge_property_type(et), g).first;
                temp_edges.push_back(std::make_pair(target_vertex, rand_new));
                --num_rewire;
                if (verbose >= 2) {
                  std::cout << " to " << new_edge1 << " and "
                            << new_edge2 << std::endl;
                  std::cout << num_rewire << " to go" << std::endl;
                }
              }
            }
          } else if (out_degree(rand_temp, cg) == 1) {
            // cannot rewire this link, so remove it from the list
            std::cout << " cannot rewire this link, so remove it from the list" << std::endl;
            temp_edges.erase(temp_edges.begin() + rand_edge);
            --num_rewire;
            if (verbose >=2) {
              std::cout << num_rewire << " to go" << std::endl;
            }
          }
        }
      }
    } else {
      std::cerr << "ERROR: rewire fraction must be between 0 and 1" << std::endl;
      std::cerr << "no rewiring performed" << std::endl;
    }
    return 0;
  }
  
  //----------------------------------------------------------
  /*! \brief Randomly remove edges from a graph.
    
  \param[in, out] g The graph to remove the edges from.
  \param[in] r The random generator to use for selecting the edges to remove.
  \param[in] removeFraction The fraction of edges to remove.
  \ingroup graph_structure
  */
  template <typename Graph, typename RandomGenerator>
  void removeEdges(Graph& g, RandomGenerator& r, double removeFraction)
  {
    typedef typename graph_traits<Graph>::edge_descriptor
      edge_descriptor;

    unsigned int N = num_edges(g);
    
    if (removeFraction > 0. && removeFraction < 1.) {
      unsigned int num_remove =
        static_cast<unsigned int>(removeFraction * N);
      
      while (num_remove > 0) {
        // select random edge for removal
        edge_descriptor e = random_edge(g, r);
        boost::remove_edge(e, g);
        --num_remove;
      }
    } else {
      std::cerr << "ERROR: remove fraction must be between 0 and 1" << std::endl;
      std::cerr << "no removing performed" << std::endl;
    }
  } 
  
  //----------------------------------------------------------
  /*! \brief Add random edges to a graph
    
  \param[in, out] g The graph to remove the edges from.
  \param[in] r The random generator to use for selecting the edges to remove.
  \param[in] addFraction The number of edges to add as a fraction of existing
  edges.
  \ingroup graph_structure
  */
  template <typename Graph, typename RandomGenerator>
  void addEdges(Graph& g, RandomGenerator& r, double addFraction)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename edge_property_type<Graph>::type
      edge_property_type;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;
    
    edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(g);

    // initialize edge type from the first edge in the graph
    unsigned int et = g[*ei].type;
    unsigned int N = num_edges(g);
    
    if (addFraction > 0.) {
      unsigned int num_add =
        static_cast<unsigned int>(addFraction * N);
      
      while (num_add > 0) {
        // select random edge for adding
        bool found = false;
        vertex_descriptor v1;
        vertex_descriptor v2;
        // find a random pair of yet unconnected vertices
        while (!found) {
          v1 = random_vertex(g, r);
          do {
            v2 = random_vertex(g, r);
          } while (v2 == v1);
          found = !(edge(v1,v2,g).second);
        }
        // add edge between the chosen vertices
        add_edge(v1,v2,edge_property_type(et),g);
        --num_add;
      }
    }
  }

  //----------------------------------------------------------
  /*! \brief Reverses an edge.
  
  \param[in] e The edge to reverse.
  \param[in] g The graph containing the edge.
  \return The reversed edge.
  \ingroup graph_structure
  */
  template <class Graph>
  typename graph_traits<Graph>::edge_descriptor
  reverse_edge(typename graph_traits<Graph>::edge_descriptor e, Graph& g)
  {
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    return edge_descriptor(target(e, g), source(e, g), e.get_property());
  }
  
  //----------------------------------------------------------
  /*! \brief Choose a random undirected edge from a graph.

  This differs from boost:random_edge in that it additionally randomizes source
  and target, whereas boost:edge returns an edge with source and target as they
  have been created.
  
  \param[in, out] g The graph to remove the edges from.
  \param[in] gen The random generator to use for selecting the edges to remove.
  \return The edge_descriptor of the randomly chosen edge.
  \ingroup graph_structure
  */
  template <class Graph, class RandomNumGen>
  typename graph_traits<Graph>::edge_descriptor
  random_edge_undirected(Graph& g, RandomNumGen& gen) {
    if (num_edges(g) > 1) {
      bool reverse = false;
      
      uniform_int<> distrib(0, (num_edges(g)*2-1));
      variate_generator<RandomNumGen&, uniform_int<> > rand_gen(gen, distrib);
      typename graph_traits<Graph>::edges_size_type
        n = rand_gen();

      if (n > num_edges(g)-1) {
        reverse = true;
        n -= num_edges(g);
      }

      typename graph_traits<Graph>::edge_iterator
        i = edges(g).first;
      while (n-- > 0) ++i; // std::advance not VC++ portable
      if (reverse) {
        return reverse_edge(*i, g);
      } else {
        return *i;
      }
    } else
      return *edges(g).first;
  }

  //----------------------------------------------------------
  /*! \brief Randomise the vertices in a graph
    
  \param[in, out] g The graph to randomise the vertices in.
  \param[in] r The random generator to use.
  \param[in] fraction The fraction of vertices to randomise.
  \ingroup graph_structure
  */
  template <typename Graph, typename RandomGenerator>
  void randomise_vertices(Graph& g, RandomGenerator& r, double fraction = 1.)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    typedef typename graph_traits<Graph>::edge_iterator
      edge_iterator;

    uniform_real<> uni_dist(0, 1);
    variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);

    Graph temp_graph;
    // add vertices to temp_graph
    boost::add_vertices(temp_graph, num_vertices(g));

    // decide which vertices to randomize
    std::vector<vertex_descriptor> original_vertices;
    std::vector<vertex_descriptor> vertices_to_randomise;
    std::vector<vertex_descriptor> randomised_vertices;
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      original_vertices.push_back(*vi);
    }
    std::vector<bool> randomise_flag(num_vertices(g), false);
    for (unsigned int i = 0; i < fraction*num_vertices(g); ++i) {
      vertex_descriptor rand_no =
        static_cast<unsigned int>(uni_gen() * original_vertices.size());
      randomise_flag[original_vertices[rand_no]] = true;
      vertices_to_randomise.push_back(original_vertices[rand_no]);
      original_vertices.erase(original_vertices.begin()+rand_no);
    }
    
    // randomise vertices and store in vector
    for (unsigned int i = 0; i < num_vertices(g); ++i) {
      if (randomise_flag[i]) {
        vertex_descriptor rand_no =
          static_cast<unsigned int>(uni_gen() * vertices_to_randomise.size());
        randomised_vertices.push_back(vertices_to_randomise[rand_no]);
        vertices_to_randomise.erase(vertices_to_randomise.begin()+rand_no);
      } else {
        randomised_vertices.push_back(original_vertices[0]);
        original_vertices.erase(original_vertices.begin());
      }
    }

    // add edges according to randomisation
    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      add_edge(randomised_vertices[source(*ei, g)],
               randomised_vertices[target(*ei, g)],
               g[*ei], temp_graph);
    }

    g = temp_graph;
    
    return;
  }

  //----------------------------------------------------------
  /*! \brief Mark all parallel edges in a graph
    
  Loops over all edges in a graph and sets the parallel_edges flag for all
  edges for which there exists an additional edge between the two vertices
  connected by it. Neede by count_parallel_edges.
  
  \param[in] g The graph to mark the parallel edges in.
  \return The number double edges in the graph.
  \ingroup helper_functions
  */
  template <typename Graph>
  unsigned int mark_parallel_edges(Graph& g)
  {
    typename boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    typename boost::graph_traits<Graph>::out_edge_iterator oi, oi_end;
    typename Graph::vertex_descriptor s, t;
    
    // counter
    unsigned int parallel_edges = 0;      
    
    // loop over all edges
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      // setting source and target vertices
      s = source(*ei, g);
      t = target(*ei, g);
      
      // check if more than one edge is going to the same target
      // from the edge's source
      
      unsigned int parallel_outgoing = 0;
      
      for (boost::tie(oi, oi_end) = out_edges(s, g); oi != oi_end; oi++) {
        if (target(*oi, g) == t) parallel_outgoing++;            
      }
      
      if (parallel_outgoing > 1) {
        
        // count parallel edges
        parallel_edges++;
        
        // mark both parallel edges
        g[*ei].parallel = true;
      }      
    }
    
    // since we looped over all edges, we have counted parallels twice
    parallel_edges = parallel_edges / 2;
    
    return parallel_edges;
    
  }

  // Utility function to fill in edge index property correctly on a graph
  // without
  // it.
  template <typename Graph>
  void setup_edge_index_map( Graph& g ) {
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    int i = 0;
    BOOST_FOREACH( edge_descriptor e, edges(g) ) {
      put(&Edge::index, g, e, i++ );
    }
  }

} // namespace boost

//----------------------------------------------------------
  

#endif
