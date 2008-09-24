/*! \file graph_io.hh
  \brief Functions for reading and writing graphs.
*/

#ifndef GRAPH_IO_HH
#define GRAPH_IO_HH

#include <boost/graph/graphviz.hpp>

#include <iostream>
#include <string>
#include <sstream>

#define NO_FREETYPE // don't include freetype fonts
#include <pngwriter.h>
  
#include "Edge.hh"

//! \addtogroup graph_io Graph I/O
//! \addtogroup graph_statistics Graph statistics
//! \addtogroup helper_functions Helper functions

namespace boost {
  
  //----------------------------------------------------------
  /*! \brief Vertex writer class which does not write vertex properties.
    \ingroup graph_io
  */
  class vertex_writer {
      
  public:

    /*! \brief Constructor
     */
    vertex_writer() {}

    //! Operator for writing vertex entries to a graphviz files
    template<typename VertexID>
    void operator()(std::ostream& out, VertexID v) const
    {
      out << "[label=\"\"]";
    }
  
  };


  //----------------------------------------------------------
  /*! \brief Vertex writer class for graphviz output which colours vertices
    according to a model providing getVertexStates().
    \ingroup graph_io
  */
  template <typename VertexProperty, typename Model>
  class vertex_state_writer {
  
  public:
  
    /*! \brief Constructor
    
    \param[in] newState The state of the vertex which is to be written to a
    graphviz file
    \param[in] newModel The model to be used for getting the graphviz option
    corresponding to a vertex state
    */
    vertex_state_writer(VertexProperty newState, const Model& newModel)
      : state(newState), model(newModel) {}

    //! Operator for writing vertex entries to a graphviz files
    template <class VertexID>
    void operator()(std::ostream& out, const VertexID& v) const
    {
      std::stringstream rgbString;
      double darkening = 1 - (1 - state[v].detail)*4/5;
      int red = static_cast<unsigned int>
        (model.getVertexStates()[state[v].base].getRGB(0) * darkening);
      int green = static_cast<unsigned int>
        (model.getVertexStates()[state[v].base].getRGB(1) * darkening);
      int blue = static_cast<unsigned int>
        (model.getVertexStates()[state[v].base].getRGB(2) * darkening);
      
      rgbString << "#" 
                << std::hex << std::setw(2) << std::setfill('0') 
                << (red > 0 ? red : 0)
                << std::hex << std::setw(2) << std::setfill('0') 
                << (green > 0 ? green : 0)
                << std::hex << std::setw(2) << std::setfill('0') 
                << (blue > 0 ? blue : 0);

      if (model.getEdgeTypes()[state[v].base].getDrawOption().find("fillcolor") != std::string::npos) {
        out << "[" << model.getEdgeTypes()[state[v].base].getDrawOption();
      } else {
        out << "[fillcolor=\"" << rgbString.str() << "\"";
      }
      out << " label=\"\"]";
    }
    
  private:
  
    VertexProperty state; //!< The state of the vertex 
    //!< The model used for determining graphviz options
    const Model& model;
  
  };

  //----------------------------------------------------------
  /*! \brief Wrapper function for the generation of a vertex_writer
    \ingroup graph_io
  */
  inline vertex_writer
  make_vertex_writer() {
    return vertex_writer();
  }

  //----------------------------------------------------------
  /*! \brief Wrapper function for the generation of a vertex_state_writer.
  
  \param[in] v The state to be passed to the vertex_writer
  \param[in] m The model to be passed to the vertex_writer
  \ingroup graph_io
  */
  template<class VertexProperty, class Model>
  inline vertex_state_writer<VertexProperty, Model>
  make_vertex_writer(VertexProperty v, const Model& m) {
    return vertex_state_writer<VertexProperty, Model>(v, m);
  }

  //----------------------------------------------------------
  /*! \brief Edge writer class for graphviz output.
  
  \ingroup graph_io
  */
  template <class EdgeProperty>
  class edge_writer {

  public:

    /*! \brief Constructor
    
    \param[in] newEdge The type of the edge which is to be written to a
    graphviz file
    \param[in] s Whether the edge type should be written or not.
    */
    edge_writer(EdgeProperty newEdge, bool s)
      : edge(newEdge), saveType(s)
    {}

    //! Operator for edge entries in graphviz files
    template <class EdgeID>
    void operator()(std::ostream& out, const EdgeID& e) const
    {
      out << "[";
      if (saveType) out << "type=" << edge[e] << " ";
//      out << "len=0.1 color=white" << "]";
      out << "color=white" << "]";
    }

  private:

    EdgeProperty edge; //!< The type of the edge
    bool saveType; //!< Whether to save the edge type into the file or not

  };

  //----------------------------------------------------------
  /*! \brief Edge writer class for graphviz output which writes edges according to
    the types given by a model providing getEdgeTypes.
  
    \ingroup graph_io
  */
  template <class EdgeProperty, class Model>
  class edge_type_writer {

  public:
    /*! \brief Constructor
    
    \param[in] newEdge The type of the edge which is to be written to a
    graphviz file
    \param[in] newModel The model to be used for getting the graphviz option
    corresponding to an edge type
    */
    edge_type_writer(EdgeProperty newEdge, const Model& newModel)
      : edge(newEdge), model(newModel) {}
  
    //! Operator for edge entries in graphviz files
    template <class EdgeType>
    void operator()(std::ostream& out, const EdgeType& e) const
    {
      out << "[" << model.getEdgeTypes()[edge[e]].getDrawOption()
          << " len=0.1 color=white" << "]";
    }
  
  private:
  
    EdgeProperty edge; //!< The type of the edge
    const Model& model; //!< The model used for determining graphviz options

  };

  //----------------------------------------------------------
  /*! \brief Wrapper function for the generation of an edge_writer
  
  \param[in] e The edge type to be passed to the edge_writer
  \param[in] s Whether to put the type into the file or not
  \ingroup graph_io
  */
  template <class EdgeProperty>
  inline edge_writer<EdgeProperty>
  make_edge_writer(EdgeProperty e, bool s) {
    return edge_writer<EdgeProperty>(e, s);
  }

  //----------------------------------------------------------
  /*! \brief Wrapper function for the generation of an edge_type_writer
  
  \param[in] e The type to be passed to the edge_type_writer
  \param[in] m The model to be passed to the edge_type_writer
  \ingroup graph_io
  */
  template <class EdgeProperty, class Model>
  inline edge_type_writer<EdgeProperty, Model>
  make_edge_writer(EdgeProperty e, const Model& m) {
    return edge_type_writer<EdgeProperty, Model>(e, m);
  }

  //----------------------------------------------------------
  /*! \brief The graph writer class for graphviz output
  
  \ingroup graph_io
  */
  struct graph_writer {
    /*! Constructor
      \param[in] t The title of the graph which is to be written to a graphviz
      file
    */
    graph_writer(std::string t = "") : title(t) {;}
  
    //! Operator for the graph entry in graphviz files
    void operator()(std::ostream& out) const {
      // graph options
//      out << "graph [overlap=scalexy bgcolor=black outputorder=edgesfirst";
      out << "graph [overlap=false bgcolor=black outputorder=edgesfirst";
      // graph title
      if (title.size() > 0) out << " label=\"" << title << "\" labelloc=\"t\"";
      // more graph options
      out << " fontcolor=white normalize=true]" << std::endl;
      // default node options
//      out << "node [shape=circle style=filled width=0.2 height=0.2 "
//          << "fixedsize=true]" << std::endl;
      out << "node [shape=circle style=filled]" << std::endl;
    }

    std::string title; //!< The title of the graph
  };

  //----------------------------------------------------------
  /*! \brief Write graph to a graphviz file
  
  Uses the boost function write_graphviz to save a graph to a file

  \param[in] g The graph to write to the file
  \param[in] fileName The name of the file to be written
  \param[in] saveType Whether the edge type should be saved in the file.
  \ingroup graph_io
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
  /*! \brief Write graph to a graphviz file according to a given model.
  
  Uses the boost function write_graphviz to save a graph to a file

  \param[in] g The graph to write to the file
  \param[in] m The model to be used to determine the graphviz options
  corresponding to vertex and edge properties
  \param[in] fileName The name of the file to be written
  \param[in] timeLabel Optional number to be written to the graph as title,
  corresponding to a time
  \ingroup graph_io
  */
  template <typename Graph, typename Model>
  void write_graph(const Graph& g, std::string fileName,
                   const Model& m, double timeLabel=0.)
  {
    typedef typename boost::vertex_property_type<Graph>::type::value_type
      vertex_property_type;
    typedef typename boost::edge_property_type<Graph>::type::value_type
      edge_property_type;
  
    std::ofstream out(fileName.c_str());
    std::stringstream s;
    if (timeLabel > 0) s << "Time: " << timeLabel;
    if (timeLabel < 0) s << "Time: start";

    write_graphviz(out, g,
                   make_vertex_writer(get(&vertex_property_type::state, g), m),
                   make_edge_writer(get(&edge_property_type::type, g), m),
                   graph_writer(s.str()));
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

    if (optionPos == std::string::npos) {
      return "";
    } else {
      // put whatever comes after option= in optionString
      std::string optionString = line.substr(optionPos+option.size()+1);
      
      // find first blank after option value
      std::string::size_type endPos = optionString.find_first_of(" ]");
      
      // return option value
      if (endPos == std::string::npos) {
        return optionString;
      } else {
        return optionString.substr(0, endPos);
      }
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
  \return The number of edges which have been read.
  \ingroup graph_io
  */
  template <typename Graph>
  int read_graph(Graph& g, const std::string graphFileName,
                          unsigned int edgeType)
  {

    std::vector<std::string> edgeStyles;
    edgeStyles.push_back("style=\"solid\"");
    edgeStyles.push_back("style=\"dashed\"");
    edgeStyles.push_back("style=\"dotted\"");
  
    // open file
    std::ifstream file;
    try {      
      file.open(graphFileName.c_str(), std::ios::in);
    }
    catch (std::exception &e) {
      return -1;
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
  
    int edgeCount = 0;
  
    if (file.is_open()) {
      while(!file.eof()) {
        //read line
        getline(file, line);
      
        // if edge or vertex
        if (line.find("];") != std::string::npos) { // found ];
        
          std::string::size_type lpos = line.find("--");
          std::string::size_type bpos = line.find("[");
        
          // if edge
          if (lpos != std::string::npos) {
          
            // extract source
            std::string s = line.substr(0, lpos);
          
            // extract target
            std::string t = line.substr(lpos+2, bpos-lpos-3);

            // set default if no edge type is given 
            unsigned int type = edgeType;
            if (line.find("style") != std::string::npos) {
              type = style2type[extractDrawOption("style", line)];
            } else if (line.find("type") != std::string::npos) {
              type = cast_stream<unsigned int>(extractDrawOption("type", line));
            } 
          
            // typecasting string to int
            unsigned int src = cast_stream<unsigned int>(s);
            unsigned int trg = cast_stream<unsigned int>(t);
          
            // check if desired type and valid edge
            if (type == edgeType &&
                src < num_vertices(g) &&
                trg < num_vertices(g)) {
                
              // add edge using Edge constructor
              add_edge(src, trg, Edge(type), g);
            
              // count
              ++edgeCount;
            }
          
          } else { // if vertex

            if (addVertices) {
              std::string s = line.substr(0, bpos);
              unsigned int src = cast_stream<unsigned int>(s);
            
	      // add vertices until we have enough to accomodate what is in graph file
              while (src >= num_vertices(g)) add_vertex(g);
            }
          }
        }
      }
    } else { // is_open = false
      return -1;
    }
  
    file.close();
  
    return edgeCount;
  
  }

  //----------------------------------------------------------
  /*! \brief Read initial states from a graphviz file
  
  Reads a graph from a graphviz file by going through the file line-by-line. The
  states are associated with states according to a given model.
  
  \param[out] g The graph to read into
  \param[in] graphFileName The name of the file to be read
  \param[in] m The model to be used to assign vertex colours to states.
  \param[in] defaultState The default state of vertices in the graph
  \return The number of edges which have been read.
  \ingroup graph_io
  */
  template <typename Graph, typename Model>
  int read_initial_graph(Graph& g, const std::string graphFileName,
                         const Model& m, unsigned int defaultState = 0)
  {

    // open file
    std::ifstream file;
    try {      
      file.open(graphFileName.c_str(), std::ios::in);
    }
    catch (std::exception &e) {
      return -1;
    }

    int vertexCount = 0;
    std::string line = "";

    // initialize color2state vertex map
    typedef std::map<std::string, unsigned int> ColorToState;
    ColorToState color2state;
    typedef std::map<unsigned int, unsigned int> RgbToState;
    RgbToState rgb2state;
    
    for (unsigned int i = 0; i < m.getVertexStates().size(); i++) {
      std::string color = extractDrawOption("fillcolor",
                                            m.getVertexStates()[i].getDrawOption());
      if (color == "") {
        // color information is rgb type
        unsigned int colorCode =
          (m.getVertexStates()[i].getRGB(0) > 0)*1 +
          (m.getVertexStates()[i].getRGB(1) > 0)*2 +
          (m.getVertexStates()[i].getRGB(2) > 0)*4;
        rgb2state.insert(std::make_pair(colorCode, i));
      } else {
        // color information is name type
        color2state.insert(std::make_pair(color,i));
      }
    }
    
    if (file.is_open()) {
      while(!file.eof()) {
        //read line
        getline(file, line);
        
        // if edge or vertex
        if (line.find("];") != std::string::npos) { // found edge or vertex
          std::string::size_type lpos = line.find("--");
          std::string::size_type bpos = line.find("[");
          if (lpos == std::string::npos) { // found vertex
            // extract vertex index
            std::string s = line.substr(0, bpos);
            
            // set default if no drawoption is given
            unsigned int state = defaultState;
            double detail = 1.;
            if (line.find("fillcolor") != std::string::npos) {
              std::string color = extractDrawOption("fillcolor", line);
              if (color.find("#") != std::string::npos) {
                // color information is rgb type
                unsigned int red, green, blue;
                std::stringstream(color.substr(2,2)) >> std::hex >> red;
                std::stringstream(color.substr(4,2)) >> std::hex >> green;
                std::stringstream(color.substr(6,2)) >> std::hex >> blue;
                state = rgb2state[(red > 0)*1 + (green > 0)*2 + (blue > 0) * 4];
                detail = 1 - 5/4. * (1 - (red+green+blue)/ static_cast<double>
                                    (m.getVertexStates()[state].getRGB(0) +
                                     m.getVertexStates()[state].getRGB(1) +
                                     m.getVertexStates()[state].getRGB(2)));
                // brush up rounding errors
                detail = detail > 0 ? detail : 0;
              } else {
                // color information is name type
                state = color2state[color];
              }
            } else if (line.find("state") != std::string::npos) {
              state =
                cast_stream<unsigned int>(extractDrawOption("state", line));
            }
              
            // typecasting string to int
            unsigned int src = cast_stream<unsigned int>(s);
            
	    // add vertices until we have enough to accomodate what is in ic file
            while (src >= num_vertices(g)) add_vertex(g);
            g[src].state.base = state;
            g[src].state.detail = detail;
            ++vertexCount;
          }
        }
      }
    } else {
      return -1;
    }
          
    file.close();

    return vertexCount;
  
  }

  //----------------------------------------------------------
  /*! \brief Get the RGB colour code corresponding to an ASCII colour
    
  \param[in] colour The ASCII colour code
  \return A vector of size three, corresponding to red, green and blue colour
  intensity
  \ingroup helper_functions
  */
  std::vector<double> getColourCode(std::string colour)
  {
    std::vector<double> code(3, 0.0);
    if (colour == "00;32") {
      // dark blue
      code[2]=0.3;
    } else if (colour == "00;31") {
      // dark red
      code[0]=0.3;
    } else if (colour == "00;34") {
      // dark green
      code[1]=0.3;
    } else if (colour == "01;32") {
      // light blue
      code[2]=1.0;
    } else if (colour == "01;31") {
      // light red
      code[0]=1.0;
    } else if (colour == "01;34") {
      // light green
      code[1]=1.0;
    }
    return code;
  }

  //----------------------------------------------------------
  /*! \brief Paint lattice to a png file
    
  \param[in] g The graph to write to the file
  \param[in] m The model to be used to determine the vertex colours
  corresponding to their states
  \param[in] fileName The name of the file to be written (without extension .png)
  \param[in] timeLabel Optional number to be written to the graph as title,
  corresponding to a time
  \ingroup graph_visualisation
  */
  template <typename Graph, typename Model>
  void write_png(const Graph& g, std::string fileName, const Model& m,
                    double timeLabel = 0.)
  {   
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    
    unsigned int sideLength = static_cast<unsigned int>(sqrt(num_vertices(g)));
    
    // create canvas
    pngwriter lattice_image(sideLength, sideLength, 0.0,
                            (fileName+".png").c_str());
    unsigned x = 1;
    unsigned y = 1;
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      // get colour code corresponding to the state of the current vertex
      std::vector<double> colourCode =
        getColourCode(m.getVertexStates()[g[*vi].state.base].getColour());
      for (unsigned int i = 0; i < colourCode.size(); ++i) {
        // if we have a light colour, shade by detailed state
        if (colourCode[i] == 1.) {
          colourCode[i] -= (1 - g[*vi].state.detail)*0.8;
        }
      }
      lattice_image.plot(x,y,colourCode[0], colourCode[1], colourCode[2]);
      ++x;
      if (x > sideLength) {
        // go to next row
        x = 1;
        ++y;
      }
    }
    lattice_image.close();
  }

  template <typename Graph>
  void write_detail_dist(const Graph& g, std::string fileName)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::ofstream distFile;
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
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      distFile << g[*vi].state.detail << std::endl;
    }
  }


  template <typename Graph, typename EdgeType, typename Model>
  void write_info_dis_corr(const Graph& g, const Model& m,
                           EdgeType et, std::string fileName)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    std::ofstream corrFile;
    try {
      corrFile.open(fileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "... unable to open correlation output file " 
                << fileName << std::endl;
      std::cerr << "... Standard exception: " << e.what() << std::endl;      
      return;
    }

    std::vector<unsigned int> infected_distances = nearest_infected(g, et, m);
    unsigned int vertex_counter = 0;
    
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      corrFile << g[*vi].state.detail << '\t'
               << infected_distances[vertex_counter] << std::endl;
      ++vertex_counter;
    }
  }
}
  
//----------------------------------------------------------

#endif
