/*! \file degree_overlap.hh
  \brief contains routines to control degree overlap
*/
#ifndef DEGREE_OVERLAP_HH
#define DEGREE_OVERLAP_HH

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
  /*! \brief Calculate Hamiltonian for degree overlap

  Calculates the hamiltonian for given degree overlap and two 
  degree types whose degree overlap is to be considered

  \param[in] g The graph to consider.
  \param[in] J overlap factor (\sa boost::rewire_degree_overlap)
  \param[in] deg_type1 The first edge type to consider.
  \param[in] deg_type2 The second edge type to consider.
  \ingroup helper_functions
  */
  template <typename Graph, typename EdgeType>
  double overlap_hamiltonian(Graph& g, double J,
                     EdgeType deg_type1, EdgeType deg_type2)
  {
    double ham = 0;
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    // calculate correlations between degrees
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      // we want the correlation of deg_type1 and deg_type2 degrees
      ham += out_degree_type(*vi, g, deg_type1)*
        out_degree_type(*vi, g, deg_type2);
    }
    // multiply with desired assortivity factor
    ham *= -J;
    return ham;
  }

  //----------------------------------------------------------
  /*! \brief Rewire edges according to degree overlap.

  Randomly swaps links of type deg_type2 to create degree overlap between the
  deg_type1 and deg_type2 degrees of vertices --
  rewiring is controlled by parameter J: positive J leads to
  correlation between degrees, negative J to anticorrelation,
  J=0 to no correlation

  \param[in] g The graph to consider.
  \param[in] r The random generator to use.
  \param[in] ol The desired degree overlap.
  \param[in] deg_type1 The first edge type to consider.
  \param[in] deg_type2 The second edge type to consider.
  \param[in] verbose Whether to be verbose.
  \ingroup graph_structure
  */
  template <typename RandomGenerator, typename Graph, typename EdgeType>
  void rewire_degree_overlap(Graph& g, RandomGenerator& r, double ol,
                             EdgeType deg_type1, EdgeType deg_type2,
                             unsigned int verbose = 0)
  {
    typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge_descriptor;
    typedef typename std::vector<edge_descriptor>
      edge_vector;
    typedef typename std::vector<edge_descriptor>::iterator
      edge_vector_iterator;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator
      out_edge_iterator;
    typedef typename edge_property_type<Graph>::type
      edge_property_type;
    
    unsigned int steps = 0;

    uniform_real<> uni_dist(0, 1);
    variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);

    double current_degree_overlap =
      degree_overlap(g, deg_type1, deg_type2);
      
    if (verbose >= 1) {
      std::cout << "Rewiring graph with degree overlap "
                << current_degree_overlap << std::endl;
    }
    
    double J = (ol > current_degree_overlap) ? 1 : 0;
    
    double H = overlap_hamiltonian(g, J, deg_type1, deg_type2);

    bool converged = (fabs(current_degree_overlap - ol) < 0.01);
    bool converging = true;
    
    for (unsigned int i = 0; !converged; i++) {
      // choose two vertices for considering swapping links of deg_type2
      vertex_descriptor vertex1, vertex2;

      // choose first random vertex for swapping
      vertex1 = random_vertex(g, r);

      // choose second random vertex for swapping
      do {
        vertex2 = random_vertex(g, r);
      } while (vertex1 == vertex2);
      
      // calculated hamiltonian as it would be after rewiring
      double s1 = out_degree_type(vertex1, g, deg_type1) *
        out_degree_type(vertex2, g, deg_type2);
      double s2 = out_degree_type(vertex2, g, deg_type1) *
        out_degree_type(vertex1, g, deg_type2);
      double s3 = out_degree_type(vertex1, g, deg_type1) *
        out_degree_type(vertex1, g, deg_type2);
      double s4 = out_degree_type(vertex2, g, deg_type1) *
        out_degree_type(vertex2, g, deg_type2);

      double sum = s1 + s2 - s3 - s4;

      double temp_H = H - J * sum;
      
      // accept change with probability min(1, exp(-(H(G')-H(G))))
      bool accept = false;
      double prob = exp(H - temp_H);
      if (prob > 1 || uni_gen() < prob) {
        accept = true;
      }
      // if we do accept the change, do the swap
      if (verbose >= 2) {
        std::cout << i << " (" << vertex1 << ")<->("
                  << vertex2 << ") ";
        if (accept) std::cout << "accepted";
        else std::cout << "rejected";
        std::cout << " (swapping probability " << (prob > 1. ? 1. : prob)
                  << ", H=" << H << ", H'=" << temp_H << ")" << std::endl;
      }
      
      if (accept) {
        H = temp_H;
        // swap all edges of type deg_type2 between vertex1 and vertex2
        edge_vector edges1, edges2;
        out_edge_iterator oi, oi_end;

        // store existing edges
        for (tie(oi, oi_end) = out_edges(vertex1, g); oi != oi_end; oi++) {
          if (g[*oi].type == deg_type2 && target(*oi, g) != vertex2) edges1.push_back(*oi);
        }
        for (tie(oi, oi_end) = out_edges(vertex2, g); oi != oi_end; oi++) {
          if (g[*oi].type == deg_type2 && target(*oi, g) != vertex1) edges2.push_back(*oi);
        }

        // swap edges
        for (edge_vector_iterator it = edges1.begin();
             it != edges1.end(); it++) {
          remove_edge(*it, g);
          add_edge(vertex2, target(*it, g), edge_property_type(deg_type2), g);
        }
        for (edge_vector_iterator it = edges2.begin();
             it != edges2.end(); it++) {
          remove_edge(*it, g);
          add_edge(vertex1, target(*it, g), edge_property_type(deg_type2), g);
        }
      }
      
      // check if degree_overlap is converging
      ++steps;
      if (steps%num_vertices(g) == 0) {
        double new_degree_overlap =
          degree_overlap(g, deg_type1, deg_type2);

        // converge if degree overlap is going into the wrong direction for two
        // times in a row
        bool now_converging =
          (fabs(ol-new_degree_overlap) < fabs(ol-current_degree_overlap));
        converged = !now_converging && !converging;
        converging = now_converging;

        current_degree_overlap = new_degree_overlap;

        if (verbose >= 1) {
          std::cout << "degree overlap " << new_degree_overlap << " step " << steps
                    << std::endl;
        }
        
      }
    }

    if (verbose >= 2) {
      std::cout << "Rewiring done: degree overlap "
                << current_degree_overlap
                << " obtained after " << steps << " steps" << std::endl;
    }
  }

  //----------------------------------------------------------
  /*! \brief Calculate degree overlap.

  Calculates the degree overlap as Pearson correlation coefficient of the two
  degrees of vertices.
  
  \param[in] g The graph to consider.
  \param[in] deg_type1 The first edge type to consider. 
  \param[in] deg_type2 The second edge type to consider.
  \return The correlation between the two degrees in vertices.
  */
  template <typename Graph, typename EdgeType>
  double degree_overlap(Graph& g, EdgeType deg_type1, EdgeType deg_type2)
  {

    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;
    
    double corr_sum(0.), deg_sum_1(0.), deg_sum_2(0.);
    double deg_sq_sum_1(0.), deg_sq_sum_2(0.);

    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
      corr_sum += out_degree_type(*vi, g, deg_type1)*
        out_degree_type(*vi, g, deg_type2);
      deg_sum_1 += out_degree_type(*vi, g, deg_type1);
      deg_sum_2 += out_degree_type(*vi, g, deg_type2);
      deg_sq_sum_1 += pow(out_degree_type(*vi, g, deg_type1), 2);
      deg_sq_sum_2 += pow(out_degree_type(*vi, g, deg_type2), 2);
    }

    double correlation =
      (corr_sum - (deg_sum_1*deg_sum_2)/(double)num_vertices(g))/
      sqrt((deg_sq_sum_1 - pow(deg_sum_1, 2)/(double)num_vertices(g))*
           (deg_sq_sum_2 - pow(deg_sum_2, 2)/(double)num_vertices(g)));
    
    return correlation;
  }
  
} // namespace boost

//----------------------------------------------------------
  
#endif
