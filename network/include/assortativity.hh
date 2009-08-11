/*! \file assortativity.hh
  \brief contains routines to control assortativity
*/
#ifndef ASSORTATIVITY_HH
#define ASSORTATIVITY_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
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
        ham += out_degree_type(source(*ei, g), g, deg_type1)*
          out_degree_type(target(*ei, g), g, deg_type2) +
          out_degree_type(target(*ei, g), g, deg_type1) *
          out_degree_type(source(*ei, g), g, deg_type2);
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
  \param[in] ass The desired assortativity.
  \param[in] et The edge type at the end of which vertices whill be considered.
  \param[in] deg_type1 The first adjacent edge type to consider.
  \param[in] deg_type2 The second adjacent edge type to consider.
  \param[in] verbose Whether to be verbose.
  \ingroup graph_structure
  */
  template <typename RandomGenerator, typename Graph, typename EdgeType>
  void rewire_assortatively(Graph& g, RandomGenerator& r, double ass,
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

    unsigned int steps = 0;

    double current_assortativity =
      assortativity(g, et, deg_type1, deg_type2);
    if (verbose >= 1) {
      std::cout << "Rewiring graph with assortativity "
                << current_assortativity << std::endl;
    }

    double J = (ass > current_assortativity) ? 1 : -1;
    
    double H = assortativity_hamiltonian(g, J, et, deg_type1, deg_type2);

    bool converged = (fabs(current_assortativity - ass) < 0.01);
    bool converging = true;

    // we will need to choose a lot of random edges, so we do not want to use
    // the slow boost::random_edge. Instead, we create a vector of edges from
    // which we can pick elements quickly

    std::vector<edge_descriptor> et_edges;

    edge_iterator ei, ei_end;
    
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
      if (g[*ei].type == et) et_edges.push_back(*ei);
    }

    unsigned int nEdges = et_edges.size();
    
    uniform_int<> int_dist(0, 2*nEdges-1);
    variate_generator<RandomGenerator&, uniform_int<> > int_gen(r, int_dist);
    uniform_real<> uni_dist(0, 1);
    variate_generator<RandomGenerator&, uniform_real<> > uni_gen(r, uni_dist);

    for (unsigned int i = 0; !converged; i++) {
      // choose two random edges of type et which do not share a node for
      // rewiring

      edges_size_type n1, n2;
    
      edge_descriptor edge1, edge2, rewired1, rewired2;
      vertex_descriptor source1, target1, source2, target2;

      
      // choose first random edge for swapping
      n1 = int_gen();
      if (n1 > nEdges-1) {
        n1 -= nEdges;
        edge1 = reverse_edge(et_edges[n1], g);
      } else {
        edge1 = et_edges[n1];
      }
      
      source1 = source(et_edges[n1], g);
      target1 = target(et_edges[n1], g);

      // choose second random edge for swapping
      bool valid_swap;
      do {
        n2 = int_gen();
        if (n2 > nEdges-1) {
          n2 -= nEdges;
          edge2 = reverse_edge(et_edges[n2], g);
        } else {
          edge2 = et_edges[n2];
        }
      
        source2 = source(et_edges[n2], g);
        target2 = target(et_edges[n2], g);

        // make sure we have 4 distinct vertices
        valid_swap = !(source1 == source2 || source1 == target2 ||
                       target1 == source2 || target1 == target2);
        
        if (valid_swap) {
          // make sure that the none of the possibly new edges is already
          // there 
          if (edge(source1, target2, g, et).second ||
              edge(source2, target1, g, et).second) {
              valid_swap = false;
          }
        }
      } while (!valid_swap);
      
      // calculated hamiltonian as it would be after rewiring
      double s1_d1 = out_degree_type(source1, g, deg_type1);
      double s1_d2 = out_degree_type(source1, g, deg_type2);
      double t1_d1 = out_degree_type(target1, g, deg_type1);
      double t1_d2 = out_degree_type(target1, g, deg_type2);
      double s2_d1 = out_degree_type(source2, g, deg_type1);
      double s2_d2 = out_degree_type(source2, g, deg_type2);
      double t2_d1 = out_degree_type(target2, g, deg_type1);
      double t2_d2 = out_degree_type(target2, g, deg_type2);

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
      }
      
      // check if assortativity is converging
      ++steps;    
      if (steps%num_vertices(g) == 0) {
        double new_assortativity =
          assortativity(g, et, deg_type1, deg_type2);

        // converge if assortativity is going into the wrong direction for two
        // times in a row
        bool now_converging =
          (fabs(ass-new_assortativity) < fabs(ass-current_assortativity));
        converged = !now_converging && !converging;
        converging = now_converging;

        current_assortativity = new_assortativity;

        if (verbose >= 1) {
          std::cout << "assortativity " << new_assortativity << " step " << steps
                    << std::endl;
        }
      }
    }

    if (verbose >= 2) {
      std::cout << "Rewiring done: assortativity "
                << current_assortativity
                << " obtained after " << steps << " steps" << std::endl;
    }
  }

  //----------------------------------------------------------
  /*! \brief Calculate assortativity.

  Calculates assortativity in a graph as Pearson correlation coefficient of the
  two degrees. In the case of a single network, this reduces to assortativity as
  defined by Newman (DOI: 10.1103/PhysRevLett.89.208701I)

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
    
    double corr_sum(0.), deg_sum_1(0.), deg_sum_2(0.);
    double deg_sq_sum_1(0.), deg_sq_sum_2(0.);

    unsigned int count(0);

    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
      if (g[*ei].type == et) {
        // we actually add two contributions, according to both directions of
        // the edge
        count += 2;
        corr_sum += out_degree_type(source(*ei, g), g, deg_type1)*
          out_degree_type(target(*ei, g), g, deg_type2) +
          out_degree_type(target(*ei, g), g, deg_type1) *
          out_degree_type(source(*ei, g), g, deg_type2);
        deg_sum_1 += out_degree_type(source(*ei, g), g, deg_type1) +
          out_degree_type(target(*ei, g), g, deg_type1);
        deg_sum_2 += out_degree_type(source(*ei, g), g, deg_type2) +
          out_degree_type(target(*ei, g), g, deg_type2);
        deg_sq_sum_1 += pow(out_degree_type(source(*ei, g), g, deg_type1), 2) +
          pow(out_degree_type(target(*ei, g), g, deg_type1), 2);
        deg_sq_sum_2 += pow(out_degree_type(source(*ei, g), g, deg_type2), 2) +
          pow(out_degree_type(target(*ei, g), g, deg_type2), 2);
      }
    }

    double correlation =
      (corr_sum - (deg_sum_1*deg_sum_2)/(double)count)/
      sqrt((deg_sq_sum_1 - pow(deg_sum_1, 2)/(double)count)*
           (deg_sq_sum_2 - pow(deg_sum_2, 2)/(double)count));
    
    return correlation;
  }
  
  //----------------------------------------------------------
  /*! \brief Calculate e_jk's.

  */
  template <typename Graph>
  std::string print_corr_matrices(Graph& g, unsigned int nEdgeTypes = 2)
  {

    typedef typename boost::graph_traits<Graph>::edge_iterator
      edge_iterator;
    typedef boost::multi_array<boost::numeric::ublas::matrix<unsigned int>, 3>
      array_type;
    array_type e(boost::extents[nEdgeTypes][nEdgeTypes][nEdgeTypes]);

    typedef boost::multi_array<unsigned int, 3>
      degree_type;
    degree_type d(boost::extents[nEdgeTypes][nEdgeTypes][0]);
    std::stringstream s;

    edge_iterator ei, ei_end;

    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {

      std::vector<unsigned int> source_degree(nEdgeTypes);
      std::vector<unsigned int> target_degree(nEdgeTypes);
      
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        source_degree[i] = out_degree_type(source(*ei, g), g, i);
        target_degree[i] = out_degree_type(target(*ei, g), g, i);
	unsigned int maxdeg = std::max(source_degree[i], target_degree[i]);
	if (maxdeg + 1 >  d.shape()[2]) {
	  d.resize(boost::extents[nEdgeTypes][nEdgeTypes][maxdeg + 1]);
	}
	++d[g[*ei].type][i][source_degree[i]];
	++d[g[*ei].type][i][target_degree[i]];
      }

      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        for (unsigned int j = i; j < nEdgeTypes; ++j) {
	  
	  unsigned int size1 = e[g[*ei].type][i][j].size1();
	  unsigned int size2 = e[g[*ei].type][i][j].size2();
	  
	  unsigned int maxi = std::max(source_degree[i], target_degree[i]);
	  
	  if (maxi+1 > size1) {
	    e[g[*ei].type][i][j].resize(maxi+1, size2, true);
	    for (unsigned k = size1; k < maxi+1; ++k) {
	      for (unsigned l = 0; l < size2; ++l) {
		e[g[*ei].type][i][j](k, l) = 0;
	      }
	    }
	    size1 = maxi+1;
	  }
	  
	  unsigned int maxj = std::max(source_degree[j], target_degree[j]);
	  if (maxj+1 > size2) {
	    e[g[*ei].type][i][j].resize(size1, maxj+1, true);
	    for (unsigned k = size2; k < maxj+1; ++k) {
	      for (unsigned l = 0; l < size1; ++l) {
		e[g[*ei].type][i][j](l, k) = 0;
	      }
	    }
	  }

          ++e[g[*ei].type][i][j](source_degree[i],target_degree[j]);
          ++e[g[*ei].type][i][j](target_degree[i],source_degree[j]);
        }
      }
    }
    
    for (unsigned int et = 0; et < nEdgeTypes; ++et) {
      for (unsigned int i = 0; i < nEdgeTypes; ++i) {
        for (unsigned int j = i; j < nEdgeTypes; ++j) {
	  s << et << " " << i << " " << j << ":" << std::endl;
	  for (unsigned int k = 0; k < e[et][i][j].size1(); ++k) {
	    for (unsigned int l = 0; l < e[et][i][j].size2(); ++l) {
	      if (e[et][i][j](k, l) > 0) {
		s << k << "\t" << l << "\t" << (e[et][i][j](k, l) / static_cast<double>(d[et][i][k])) << std::endl;
	      }
	    }
	  }
        }
      }
    }

    return s.str();
  }
  
} // namespace boost

//----------------------------------------------------------
  
#endif
