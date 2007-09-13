/*! \file graph_io.hh
  \brief Functions for reading and writing graphs.
*/

#ifndef GRAPH_IO_HH
#define GRAPH_IO_HH

#include <boost/graph/graphviz.hpp>

#include <iostream>
#include <string>
#include <sstream>

#include "Edge.hh"

//! \addtogroup graph_visualisation Graph visualisation
//! \addtogroup helper_functions Helper functions

//----------------------------------------------------------
/*! \brief The vertex writer class for graphviz output
  
\ingroup graph_visualisation
*/
class vertex_writer {
      
public:

  /*! \brief Constructor
    
  \param[in] newState The state of the vertex which is to be written to a
  graphviz file
  \param[in] newModel The model to be used for getting the graphviz option
  corresponding to a vertex state
  */
  vertex_writer() {}

  //! Operator for writing vertex entries to a graphviz files
  template<typename VertexProperty>
  void operator()(std::ostream& out, VertexProperty v) const
  {
    out << "[label=\"\"]";
  }
  
};

//----------------------------------------------------------
/*! \brief Wrapper function for the generation of a vertex_writer
  
\param[in] v The state to be passed to the vertex_writer
\param[in] m The model to be passed to the vertex_writer
\ingroup graph_visualisation
*/
inline vertex_writer
make_vertex_writer() {
  return vertex_writer();
}

//----------------------------------------------------------
/*! \brief The edge writer class for graphviz output
  
\ingroup graph_visualisation
*/
template <class EdgeProperty>
class edge_writer {

public:

  /*! \brief Constructor
    
  \param[in] newEdge The type of the edge which is to be written to a
  graphviz file
  \param[in] newModel The model to be used for getting the graphviz option
  corresponding to an edge type
  */
  edge_writer(EdgeProperty newEdge, bool s)
    : edge(newEdge), saveType(s)
  {}

  //! Operator for edge entries in graphviz files
  template <class EdgeType>
  void operator()(std::ostream& out, const EdgeType& e) const
  {
    out << "[";
    if (saveType) out << "type=" << edge[e] << " ";
    out << "len=0.1 color=white" << "]";
  }

private:

  EdgeProperty edge; //!< The type of the edge
  bool saveType; //!< Whether to save the edge type into the file or not

};

//----------------------------------------------------------
/*! \brief Wrapper function for the generation of an edge_writer
  
\param[in] e The type to be passed to the edge_writer
\param[in] m The model to be passed to the edge_writer
\ingroup graph_visualisation
*/
template <class EdgeProperty>
inline edge_writer<EdgeProperty>
make_edge_writer(EdgeProperty e, bool s) {
  return edge_writer<EdgeProperty>(e, s);
}

//----------------------------------------------------------
/*! \brief The graph writer class for graphviz output
  
\ingroup graph_visualisation
*/
struct graph_writer {
  /*! Constructor
    \param[in] t The title of the graph which is to be written to a graphviz
    file
  */
  graph_writer() {;}
  
  //! Operator for the graph entry in graphviz files
  void operator()(std::ostream& out) const {
    // graph options
    out << "graph [overlap=scalexy bgcolor=black outputorder=edgesfirst";
    // more graph options
    out << " fontcolor=white normalize=true]" << std::endl;
    // default node options
    out << "node [shape=circle style=filled width=0.2 height=0.2 "
        << "fixedsize=true]" << std::endl;
  }

};

//----------------------------------------------------------
/*! \brief Write graph to a graphviz file
  
Uses the boost function write_graphviz to save a graph to a file

\param[in] g The graph to write to the file
\param[in] m The model to be used to determine the graphviz options
corresponding to vertex and edge properties
\param[in] fileName The name of the file to be written
\param[in] timeLabel Optional number to be written to the graph as title,
corresponding to a time
\ingroup graph_visualisation
*/
template <typename Graph>
void write_graph(const Graph& g, std::string fileName, bool saveType = true)
{
  typedef typename boost::edge_property_type<Graph>::type::value_type
    edge_property_type;
  
  std::ofstream out(fileName.c_str());
  std::stringstream s;

  write_graphviz(out, g,
                 make_vertex_writer(),
                 make_edge_writer(get(&edge_property_type::type, g), saveType),
                 graph_writer());
}

//----------------------------------------------------------
/*! \brief Extract the value of a draw option from a graphviz line
  
\param[in] option The option to extract
\param[in] line The line to extract the option from
\return The value of the option
\ingroup helper_functions
*/
std::string extractDrawOption(std::string option, const std::string& line) 
{
  // find option position
  std::string::size_type optionPos = line.find(option+"=");

  // put whatever comes after option= in optionString
  std::string optionString = line.substr(optionPos+option.size()+1);

  // find first blank after option value
  std::string::size_type endPos = optionString.find(" ");

  // return option value
  if (endPos == std::string::npos) {
    return optionString;
  } else {
    return optionString.substr(0, endPos);
  }
}

//----------------------------------------------------------
/*! \brief Cast a variable to another type
  
\param[in] t The variable to be cast
\return The converted variable
\ingroup helper_functions
*/
template <class OutType, class InValue>
OutType cast_stream(const InValue& t)
{
  std::stringstream ss;
  ss << t; // first insert value to stream
  OutType result; // value will be converted to OutType
  ss >> result; // write value to result
  return result;
}

//----------------------------------------------------------
/*! \brief Read graph from a graphviz file
  
Reads a graph from a graphviz file by going through the file line-by-line.
If there are already vertices in the graph, it only connects them according to
the edges given in the file. If not, it creates the vertices and, if desired,
reads in their states. 

\param[out] g The graph to read into
\param[in] graphFileName The name of the file to be read
\param[in] edgeType The edge type to be assigned to the read file
\param[in] verbose Add verbosity
\return 0 if successful, >0 if failed
\ingroup graph_visualisation
*/
template <typename Graph>
int read_graph(Graph& g, const std::string graphFileName,
               unsigned int edgeType, unsigned int verbose = 0)
{

  std::vector<std::string> edgeStyles;
  edgeStyles.push_back("solid");
  edgeStyles.push_back("dashed");
  edgeStyles.push_back("dotted");
  
  // open file
  std::ifstream file;
  try {      
    file.open(graphFileName.c_str(), std::ios::in);
  }
  catch (std::exception &e) {
    std::cerr << "Unable to open graph file: " << e.what()
              << std::endl;            
  }
  
  // map
  typedef std::map<std::string, unsigned int> StyleToType;
  StyleToType style2type;
  
  // initialize style2type edge map
  for (unsigned int i = 0; i < edgeStyles.size(); i++) {
    std::string style =
      extractDrawOption("style",edgeStyles[i]);
    style2type.insert(std::make_pair(style,i));
  }

  // initialize addVertices option
  bool addVertices = false;
  if (num_vertices(g) == 0) {
    addVertices = true;
  }

  // read file
  std::string line = "";
  std::stringstream sstr; // for typecasting
  
  std::map<unsigned int, unsigned int> edge_count;
  std::map<unsigned int, unsigned int> state_count;
  
  if (file.is_open()) {
    while(!file.eof()) {
      //read line
      getline(file, line);
      
      // if edge or vertex
      if (line.find("];") != std::string::npos) { // found ];
        
        std::string::size_type lpos = line.find("--");
        std::string::size_type bpos = line.find("[");
        std::string s, t;
        unsigned int type, src, trg;
        
        // if edge
        if (lpos != std::string::npos) {
          
          // extract source
          s = line.substr(0, lpos);
          
          // extract target
          t = line.substr(lpos+2, bpos-lpos-3);
          
          // extract type
          if (line.find("style") != std::string::npos) {
            type = style2type[extractDrawOption("style", line)];
          } else if (line.find("type") != std::string::npos) {
            std::istringstream typeStream(extractDrawOption("type", line));
            typeStream >> type;
          } else {
            // no type is given, so we automatically add
            type = edgeType;
          }
          
          // typecasting string to int
          src = cast_stream<unsigned int>(s);
          trg = cast_stream<unsigned int>(t);
          
          // check if desired type and valid edge
          if (type == edgeType &&
              src < num_vertices(g) &&
              trg < num_vertices(g)) {
                
            // add edge using Edge constructor
            add_edge(src, trg, Edge(type), g);
            
            // count
            edge_count[type]++;
          }
          
        } else { // if vertex

          if (addVertices) {
            add_vertex(g);
          }
        }
      }
    }
      
    if (verbose) {
      // print message
      std::map<unsigned int, unsigned int>::const_iterator it;
      
      std::cout << std::endl
                <<"data from read_graph file = " << graphFileName
                << std::endl;
      
      unsigned int tmpSum = 0;
        
      for (it = edge_count.begin(); it != edge_count.end(); ++it) {         
        std::cout << "# of edges of type " << it->first << " read: "
                  << it->second << std::endl;
        tmpSum += it->second;
      }
      std::cout << "Total # of edges read: " << tmpSum << std::endl
                << std::endl;
      
    }
  } else { // is_open = false
    return 1;
  }
  
  file.close();
  
  return 0;
  
}

//----------------------------------------------------------
/*! \brief Write the degree distribution to a file
  
Loops over all vertices and records both total degree distribution and degree
distribution for each edge type. 

\param[in] g The graph to calculate the degree distribution for
\param[in] m The model to be used to determine edge types
\param[in] degreeFileName The name of the file to write the degree
distribution to
\return 0 if successful, <0 if failed
\ingroup graph_visualisation
*/
template <typename Graph>
unsigned int write_degree(const Graph& g, unsigned int nEdgeTypes,
                          const std::string degreeFileName)
{
  typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
  typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;   

  // counters
  typedef std::map<unsigned int, unsigned int> degree_map;
  std::vector<degree_map*> degree(nEdgeTypes+1); // including parallel
   
  for (unsigned int i = 0; i < degree.size(); i++)
    degree[i] = new degree_map;
   
  // open file
  std::ofstream file;
  try {      
    file.open(degreeFileName.c_str(), std::ios::out);
  }
  catch (std::exception &e) {
    std::cerr << "Unable to open degree file: " << e.what()
              << std::endl;
    return -1;
  }
   
  // loop over all vertices
  for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      
    // tmp sums
    unsigned int vertex_deg[nEdgeTypes+1];
    for (unsigned int i = 0; i < nEdgeTypes+1; i++) 
      vertex_deg[i] = 0;
      
    // loop over the vertex out edges
    for (tie(ei, ei_end) = out_edges(*vi, g); ei != ei_end; ei++) {
         
      // count all edges including parallel
      vertex_deg[nEdgeTypes] += 1;
         
      // count all other edges (whether they are parallel or not)
      vertex_deg[g[*ei].type] += 1;
    }

    // update global degree
    for (unsigned int i = 0; i < degree.size(); i++)
      (*degree[i])[vertex_deg[i]] += 1;
  }

  // printing
  degree_map::const_iterator it;      
   
  // setting max_degree
  unsigned int max_degree = 0;
  for (unsigned int i =0; i < degree.size(); ++i) {
    it = (*degree[i]).end();
    --it;
    if ((*it).first > max_degree)
      max_degree = (*it).first;
  }

  file << "# max_degree = " << max_degree << std::endl;
  file << "# degree | d-edges | i-edges | total-edges\n";
  file << "# ------ | ------- | ------- | -----------\n";
   
  // enumerate over all degrees and all maps
  for (unsigned int i = 0; i <= max_degree; i++) {
    file << i;
    for (unsigned int j = 0; j < degree.size(); j++)
      file << " " << (*degree[j])[i];
    file << std::endl;
  }    
   
  // close file
  file.close();
   
  // delete
  for (unsigned int i = 0; i <= nEdgeTypes; i++) 
    delete degree[i];
   
  return 0;
   
}


//----------------------------------------------------------

#endif
