/*! \file assortativity.hh
  \brief contains routines to control assortativity
*/
#ifndef ASSORTATIVITY_HH
#define ASSORTATIVITY_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>

#include "graph_statistics.hh"

//! \addtogroup graph_structure Graph structure
//! \addtogroup helper_functions Helper functions
//! \addtogroup graph_statistics Graph statistics

namespace boost {

  //----------------------------------------------------------
  /*! \brief Calculate assortativity Hamiltonian

  Calculates the hamiltonian for given assortativity factor and two 
  degree types whose (dis-)assortativity is to be considered

  \param[in] g The graph to consider.
  \param[in] J Assortativity factor (\sa boost::rewire_assortatively)
  \param[in] et The edge type at the end of which vertices whill be considered.
  \param[in] deg_type1 The first adjacent edge type to consider.
  \param[in] deg_type2 The second adjacent edge type to consider.
  \ingroup helper_functions
  */
  template <typename Graph, typename EdgeType>
  double assortativity_hamiltonian(Graph& g, double J, EdgeType et,
                                   EdgeType deg_type1, EdgeType deg_type2)
  {
    double ham = 0;
    typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;

    // calculate correlations between degrees
    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      // we want the correlation of deg_type1 and deg_type2 degrees, but at
      // the two ends of an et link, so we check here if type is et
      if (g[*ei].type == et) {
        ham += out_degree_type(g, source(*ei, g), deg_type1)*
          out_degree_type(g, target(*ei, g), deg_type2) +
          out_degree_type(g, target(*ei, g), deg_type1) *
          out_degree_type(g, source(*ei, g), deg_type2);
      }
    }
    // multiply with desired assortivity factor
    ham *= -J;
    return ham;
  }

  //----------------------------------------------------------
  /*! \brief Rewire edges assortatively

  Randomly rewires edges of type et to create assortativity between the
  deg_type1 and deg_type2 degrees at the two ends of an et edge --
  rewiring is controlled by parameter J: positive J leads to
  correlation, negative J to anticorrelation, J=0 to no correlation

  \param[in] g The graph to consider.
  \param[in] r The random generator to use.
  \param[in] J Assortativity factor.
  \param[in] et The edge type at the end of which vertices whill be considered.
  \param[in] deg_type1 The first adjacent edge type to consider.
  \param[in] deg_type2 The second adjacent edge type to consider.
  \param[in] verbose Whether to be verbose.
  \ingroup graph_structure
  */
  template <typename RandomGenerator, typename Graph, typename EdgeType>
  void rewire_assortatively(Graph& g, RandomGenerator& r, double J,
                            EdgeType et, EdgeType deg_type1, EdgeType deg_type2,
                            unsigned int verbose = 0)
  {
    typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_iterator
      edge_iterator;
    typedef typename boost::graph_traits<Graph>::edges_size_type
      edges_size_type;
    typedef typename edge_property_type<Graph>::type
      edge_property_type;
    
    const unsigned int max_steps = 100000;
    unsigned int steps = 0;

    boost::uniform_01<boost::mt19937, double> uni_gen(r);

    double current_assortativity =
      assortativity(g, et, deg_type1, deg_type2);
    if (verbose >= 1) {
      std::cout << "Rewiring graph with assortativity "
                << current_assortativity << std::endl;
    }
    
    double H = assortativity_hamiltonian(g, J, et, deg_type1, deg_type2);

    bool converged = (fabs(current_assortativity - J)/J < 0.1);

    // we will need to choose a lot of random edges, so we do not want to use
    // the slow boost::random_edge. Instead, we create a vector of edges from
    // which we can pick elements quickly

    std::vector<edge_descriptor> et_edges;

    edge_iterator ei, ei_end;
    
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      if (g[*ei].type == et) et_edges.push_back(*ei);
    }

    unsigned int nEdges = et_edges.size();
    
    uniform_int<> dist(0, 2*nEdges-1);
    variate_generator<RandomGenerator&, uniform_int<> > rand_gen(r, dist);

//     for (unsigned int i = 0; (i < max_steps) && (!converged); i++) {
    for (unsigned int i = 0; !converged; i++) {
      // choose two random edges of type et which do not share a node for
      // rewiring

      edges_size_type n1, n2;
    
      edge_descriptor edge1, edge2, rewired1, rewired2;
      vertex_descriptor source1, target1, source2, target2;

      
      // choose first random edge for swapping
      n1 = rand_gen();
      if (n1 > nEdges-1) {
        n1 -= nEdges;
        source1 = target(et_edges[n1], g);
        target1 = source(et_edges[n1], g);
      } else {
        source1 = source(et_edges[n1], g);
        target1 = target(et_edges[n1], g);
      }
      
      // choose second random edge for swapping
      bool valid_swap;
      do {
        n2 = rand_gen();
        if (n2 > nEdges-1) {
          n2 -= nEdges;
          source2 = target(et_edges[n2], g);
          target2 = source(et_edges[n2], g);
        } else {
          source2 = source(et_edges[n2], g);
          target2 = target(et_edges[n2], g);
        }
      
        valid_swap = (n1 != n2);
        
        // make sure we have 4 distinct vertices
        valid_swap = !(source1 == source2 || source1 == target2 ||
                       source2 == target1 || target1 == target2);
        
        if (valid_swap) {
          // make sure that the none of the possibly new edges is already
          // there 
          if (edge(g, source1, target2, et).second ||
              edge(g, source2, target1, et).second) {
            valid_swap = false;
          }
        }
      } while (!valid_swap);

      // calculated hamiltonian as it would be after rewiring
      double s1_d1 = out_degree_type(g, source1, deg_type1);
      double s1_d2 = out_degree_type(g, source1, deg_type2);
      double t1_d1 = out_degree_type(g, target1, deg_type1);
      double t1_d2 = out_degree_type(g, target1, deg_type2);
      double s2_d1 = out_degree_type(g, source2, deg_type1);
      double s2_d2 = out_degree_type(g, source2, deg_type2);
      double t2_d1 = out_degree_type(g, target2, deg_type1);
      double t2_d2 = out_degree_type(g, target2, deg_type2);

      double sum = s1_d1*t2_d2 + t2_d1*s1_d2 + s2_d1*t1_d2 + t1_d1*s2_d2 -
        s1_d1*t1_d2 - t1_d1*s1_d2 - s2_d1*t2_d2 - t2_d1*s2_d2;

      double temp_H = H - J/4 * sum;
      
      // accept change with probability min(1, exp(-(H(G')-H(G))))
      bool accept = false;
      double prob = exp(H - temp_H);
      if (prob > 1 || uni_gen() < prob) {
        accept = true;
      }
      // if we do accept the change, do the swap
      if (verbose >= 2) {
        std::cout << i << " (" << source1 << " " << target1 << ")<->("
                  << source2 << " " << target2 << ") ";
        if (accept) std::cout << "accepted";
        else std::cout << "rejected";
        std::cout << " (rewiring probability " << (prob > 1. ? 1. : prob)
                  << ", H=" << H << ", H'=" << temp_H << ")" << std::endl;
      }
      
      if (accept) {
        H = temp_H;
        remove_edge(et_edges[n1], g);
        remove_edge(et_edges[n2], g);
        et_edges[n1] =
          (add_edge(source1, target2, edge_property_type(et), g)).first;
        et_edges[n2] =
          (add_edge(source2, target1, edge_property_type(et), g)).first;
        et_edges[n1+nEdges] = et_edges[n1];
        et_edges[n2+nEdges] = et_edges[n2];
      }
      
      // check if assortativity is converging
      ++steps;
      if (steps%(num_vertices(g)*10) == 0) {
        double new_assortativity =
          assortativity(g, deg_type1, deg_type1, deg_type2);

        // less than 5% change in 10000 events
//         converged =
//           (fabs(new_assortativity - current_assortativity)/
//            current_assortativity < 0.05);
//         current_assortativity = new_assortativity;

        if (verbose >=1) {
          std::cout << "Within " << fabs(new_assortativity - J)*100/J
                    << "% of convergence (assortativity " << new_assortativity
                    << ")" << std::endl;
        }
      }
    }

    if (verbose >= 1) {
      std::cout << "Rewiring done: assortativity "
                << assortativity(g, deg_type1, deg_type1, deg_type2)
                << " obtained after " << steps << " steps" << std::endl;
    }
  }

  //----------------------------------------------------------
  /*! \brief Calculate assortativity.

  Calculates assortativity in a graph as defined by Newman
  (DOI: 10.1103/PhysRevLett.89.208701I)

  \param[in] g The graph to consider.
  \param[in] et The edge type vertices at the end of which are to be considered.
  \param[in] deg_type1 The first edge type to consider when looking at the
  degree of vertices at the end of the original edge.
  \param[in] deg_type2 The second edge type to consider when looking at the
  degree of vertices at the end of the original edge.
  \return The assortativity as defined by Newman.
  */
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
      if (g[*ei].type == et) {
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
