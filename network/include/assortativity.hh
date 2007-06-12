/******************************************************************/
// assortativity.hh
// contains routines to control assortativity
/******************************************************************/
#ifndef ASSORTATIVITY_HH
#define ASSORTATIVITY_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <math.h>

//----------------------------------------------------------

namespace boost {

  // calculate out_degree of a given type (as opposed to total out_degree)
  template <typename Graph, typename VertexDescriptor, typename EdgeType>
  unsigned int out_degree_type(Graph& g, VertexDescriptor v, EdgeType et)
  {
    unsigned int deg = 0;
    
    typename graph_traits<Graph>::out_edge_iterator oi, oi_end;
    for (tie(oi, oi_end) = out_edges(v, g); oi != oi_end; oi++) {
      if (g[*oi].type == et.type) ++deg;
    }
    return deg;
  }
  

  // calculate Hamiltonian for given assortativity factor and two
  // degree types whose (dis-)assortativity is to be considered
  template <typename Graph, typename EdgeType>
  double hamiltonian(Graph& g, double J,
                     EdgeType deg_type1, EdgeType deg_type2)
  {
    double ham = 0;
    typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;

    // calculate correlations between degrees
    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      // we want the correlation of deg_type1 and deg_type2 degrees, but at
      // the two ends of a deg_type1 link, so we check here if type is deg_type1
      if (g[*ei].type == deg_type1.type) {
        ham += out_degree_type(g, source(*ei, g), deg_type1)*
          out_degree_type(g, target(*ei, g), deg_type2);
      }
    }
    // multiply with desired assortivity factor
    ham *= -J/2;
    return ham;
  }

  // rewire Edges of deg_type2 to create assortativity between the deg_type1
  // and deg_type2 degrees at the two ends of a deg_type1 edge --
  // rewiring is controlled by parameter J: positive J leads to
  // correlation, negative J to anticorrelation, J=0 to no correlation
  template <typename RandomGenerator, typename Graph, typename EdgeType>
  void rewire_assortatively(Graph& g, RandomGenerator& r, double J,
                            EdgeType deg_type1, EdgeType deg_type2,
                            unsigned int verbose = 0)
  {
    typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    
    const unsigned int max_steps = 100000;
    unsigned int steps = 0;

    boost::uniform_01<boost::mt19937, double> uni_gen(r);

    double current_assortativity =
      assortativity(g, deg_type1, deg_type1, deg_type2);
    if (verbose >= 1) {
      std::cout << "Rewiring graph with assortativity "
                << current_assortativity << std::endl;
    }
    
    double H = hamiltonian(g, J, deg_type1, deg_type2);

    bool converged = false;
    
    for (unsigned int i = 0; (i < max_steps) && (!converged); i++) {
      // choose two random edges of type deg_type2 which do not share a node for
      // rewiring
      edge_descriptor edge1, edge2, rewired1, rewired2;
      vertex_descriptor source1, target1, source2, target2;

      // choose first random edge for swapping
      do {
        edge1 = random_edge(g, r);
      } while (g[edge1].type != deg_type2.type);
      source1 = source(edge1, g);
      target1 = target(edge1, g);

      // choose second random edge for swapping
      bool valid_swap;
      do {
        valid_swap = false;
        edge2 = random_edge(g, r);
        source2 = source(edge2, g);
        target2 = target(edge2, g);
        // make sure that the edge has the correct type
        valid_swap = (g[edge2].type == deg_type2.type);
        
        // make sure we have 4 distinct vertices
        valid_swap &= !(source1 == source2 || source1 == target2 ||
                        source2 == target1 || target1 == target2);
        
        // make sure that the none of the possibly new edges is already there
        out_edge_iterator oi, oi_end;
        for (tie(oi, oi_end) = out_edges(source1, g);
             (oi != oi_end && valid_swap); oi++) {
          if (target(*oi, g) == target2 && g[*oi].type == deg_type2.type) {
            valid_swap = false;
          }
        }
        for (tie(oi, oi_end) = out_edges(source2, g);
             (oi != oi_end && valid_swap); oi++) {
          if (target(*oi, g) == target1 && g[*oi].type == deg_type2.type) {
            valid_swap = false;
          }
        }
      } while (!valid_swap);

      // rewire edges
      remove_edge(edge1, g);
      remove_edge(edge2, g);
      rewired1 = (add_edge(source1, target2, deg_type2, g)).first;
      rewired2 = (add_edge(source2, target1, deg_type2, g)).first;

      // calculate rewired hamiltonian 
      double temp_H = hamiltonian(g, J, deg_type1, deg_type2);
      
      // accept change with probability min(1, exp(-(H(G')-H(G))))
      bool accept = false;
      double prob = exp(H - temp_H);
      if (prob > 1 || uni_gen() < prob) {
        accept = true;
      }
      if (verbose >= 2) {
        std::cout << i << " (" << source1 << " " << target1 << ")<->("
                  << source2 << " " << target2 << ") ";
        if (accept) std::cout << "accepted";
        else std::cout << "rejected";
        std::cout << " (rewiring probability " << (prob > 1. ? 1. : prob)
                  << ", H=" << H << ", H'=" << temp_H << ")" << std::endl;
      }
      
      // if we do not accept the change, revert back to previous graph
      if (accept) {
        H = temp_H;
      } else {
        remove_edge(rewired1, g);
        remove_edge(rewired2, g);
        add_edge(source1, target1, deg_type2, g);
        add_edge(source2, target2, deg_type2, g);
      }

      // check if assortativity is converging
      ++steps;
      if (steps%2000 == 0) {
        double new_assortativity =
          assortativity(g, deg_type1, deg_type1, deg_type2);
        // less than 5% change in 1000 events
        if (new_assortativity < (current_assortativity * 1.05)) {
          converged = true;
        }
        current_assortativity = new_assortativity;
      }
    }

    if (verbose >= 1) {
      std::cout << "Rewiring done: assortativity "
                << current_assortativity
                << " obtained after " << steps << " steps" << std::endl;
    }
  }

  // calculates assortativity as defined by Newman
  // (DOI: 10.1103/PhysRevLett.89.208701I)
  template <typename Graph, typename EdgeType>
  double assortativity(Graph& g, EdgeType et,
                       EdgeType deg_type1, EdgeType deg_type2)
  {

    typedef typename boost::graph_traits<Graph>::edge_iterator
      edge_iterator;
    
    double avg_corr = 0.;
    double avg_deg = 0.;
    double avg_sq_deg = 0.;

    unsigned int count = 0;

    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
      if (g[*ei].type == et.type) {
        ++count;
        avg_corr += out_degree_type(g, source(*ei, g), deg_type1)*
          out_degree_type(g, target(*ei, g), deg_type2) +
          out_degree_type(g, target(*ei, g), deg_type1) *
          out_degree_type(g, source(*ei, g), deg_type2);
        avg_deg += out_degree_type(g, source(*ei, g), deg_type1) +
          out_degree_type(g, target(*ei, g), deg_type2) +
          out_degree_type(g, source(*ei, g), deg_type2) +
          out_degree_type(g, target(*ei, g), deg_type1);
        avg_sq_deg += pow(out_degree_type(g, source(*ei, g), deg_type1), 2) +
          pow(out_degree_type(g, target(*ei, g), deg_type2), 2) +
          pow(out_degree_type(g, source(*ei, g), deg_type2), 2) +
          pow(out_degree_type(g, target(*ei, g), deg_type1), 2);
      }
    }

    avg_corr /= 2*count;
    avg_deg /= 4*count;
    avg_sq_deg /= 4*count;
    
    return (avg_corr - pow(avg_deg, 2))/(avg_sq_deg - pow(avg_deg, 2));
  }
  
} // namespace boost

//----------------------------------------------------------
  
#endif
