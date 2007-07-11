/*! \file cluster_coeffs.hh
  \brief Routines for calculation of clustering coefficients.
*/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <vector>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graph_traits.hpp>

#include "Model.hh"
#include "Vertex.hh"

/*! \brief Sparse matrix routines.
  \ingroup helper_functions
*/
namespace SparseMatrix {
  
  //--------------------------------------------------------------------
  //! brief Sparse matrix entry.
  template <class T, class C>
  struct Entry {
    
    template <class P, class Q>
    friend std::ostream& operator<< (std::ostream& os, const Entry<P, Q>& entry);
  
    template <class P, class Q>
    friend bool operator< (const Entry<P, Q>& lhs, const Entry<P, Q>& rhs);
  
    /*! \brief Constructor.
      \param[in] r row_i initialiser.
      \param[in] c col_i initialiser.
      \param[in] v value initialiser.
    */
    Entry(unsigned int r, unsigned int c, T v) : row_i(r), col_i(c), value(v) {}
    //! Destructor.
    ~Entry() {}
  
    unsigned int row_i; //!< The row index of the matrix entry.
    unsigned int col_i; //!< The column index of the matrix entry.
 
    T value; //!< The value contained in the matrix entry.
  
  };

  //--------------------------------------------------------------------
  template <class T, class C>
  std::ostream& operator<< (std::ostream& os, const Entry<T, C>& entry)
  {
    os << entry.row_i << "  " << entry.col_i << "  " << entry.value << std::endl;
    return os;
  }

  //--------------------------------------------------------------------
  template <class T, class C> 
  bool operator< (const Entry<T, C>& lhs, const Entry<T, C>& rhs)
  {
    return C::compare(lhs.row_i, lhs.col_i, rhs.row_i, rhs.col_i);
  }

  //--------------------------------------------------------------------
  //! \brief The comparison operator for row major
  struct RowMajor {
    static bool compare(unsigned int lhs_row_i, unsigned int lhs_col_i,
                        unsigned int rhs_row_i, unsigned int rhs_col_i)
    {
      if (lhs_row_i < rhs_row_i) return true;
      if (lhs_row_i > rhs_row_i) return false;
      if (lhs_col_i < rhs_col_i) return true;
      return false;  
    }  
  };

  //--------------------------------------------------------------------
  //! \brief The comparison operator for column major
  struct ColMajor {
    static bool compare(unsigned int lhs_row_i, unsigned int lhs_col_i,
                        unsigned int rhs_row_i, unsigned int rhs_col_i)
    {
      if (lhs_col_i < rhs_col_i) return true;
      if (lhs_col_i > rhs_col_i) return false;
      if (lhs_row_i < rhs_row_i) return true;
      return false;  
    }  
  };

  //--------------------------------------------------------------------
  //! \brief Class holding the sparse matrix.
  template <class T>
  struct SparseMatrix {
  
    typedef Entry<T, RowMajor> Mij;
    typedef Entry<T, ColMajor> Mji;
    typedef std::list<Mij> Row;
    typedef std::list<Mji> Col;
    typedef std::vector<Row> RowMajorMatrix;
    typedef std::vector<Col> ColMajorMatrix;
  
    /*! \brief Constructor.
      \param[in] size Matrix size.
    */
    SparseMatrix(unsigned int size) {
      mr.resize(size);
      mc.resize(size);
    }

    //! Matrix destructor
    ~SparseMatrix() {}

    //! Insert a value into the matrix.
    void insert(unsigned int i, unsigned int j, T val) {
      mr[i].push_back(Mij(i, j, val));
      mc[j].push_back(Mji(i, j, val));
    }
  
    void write_J(std::string fileName);

    //! Sort the matrix.
    void sort() {
      for (unsigned int i = 0; i < mr.size(); ++i) mr[i].sort();    
      for (unsigned int j = 0; j < mc.size(); ++j) mc[j].sort();
    }

    //! Empty the matrix
    void clear() {
      for (unsigned int i = 0; i < mr.size(); ++i) mr[i].clear();    
      for (unsigned int j = 0; j < mc.size(); ++j) mc[j].clear();
    }
    
    RowMajorMatrix mr; //!< The row major matrix.
    ColMajorMatrix mc; //!< The column major matrix.
    
  };
  
  //--------------------------------------------------------------------
  /*! \brief Write the matrix to a file.
    \param[in] fileName The name of the file to write.
  */
  template <class T>
  void SparseMatrix<T>::write_J(std::string fileName) 
  {
    std::ofstream file;          
  
    // open file
    try {      
      file.open(fileName.c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open Jd file: " << e.what() << std::endl;
    }
    
    // print row major
    typename Row::const_iterator aij;
    for (unsigned int i = 0; i < mr.size(); ++i)
      for (aij = mr[i].begin(); aij != mr[i].end(); aij++)
        file << *aij;    
  
    // close file
    file.close();
  
  
    //     // column major
    //     typename Col::const_iterator aji;
    //     for (unsigned int j = 0; j < mc.size(); ++j)
    //       for (aji = mc[j].begin(); aji != mc[j].end(); aji++)
    //         std::cout << *aji;    
  
  }

  //--------------------------------------------------------------------
  /*! \brief Multiply two sparse matrices.
    \param[in] A The first matrix to multiply.
    \param[in] B The second matrix to multiply.
    \param[out] C The result of the multiplication.
  */
  template <class T>
  void prod(SparseMatrix<T>& A, SparseMatrix<T>& B, SparseMatrix<T>& C)
  {  
    typename SparseMatrix<T>::Row::iterator row_it, row_end;
    typename SparseMatrix<T>::Col::iterator col_it, col_end;
    T sum;
  
    for (unsigned int i = 0; i < A.mr.size(); i++)
      for (unsigned int j = 0; j < B.mc.size(); j++) {
        row_it = A.mr[i].begin();
        row_end = A.mr[i].end();
        col_it = B.mc[j].begin();
        col_end = B.mc[j].end();
        sum = 0;
        bool assign = false;
      
        while ((row_it != row_end) && (col_it != col_end)) {
          if (row_it->col_i == col_it->row_i) {
            sum += row_it->value * col_it->value;
            row_it++;
            col_it++;
            assign = true;
          }
          else if(row_it->col_i < col_it->row_i)
            row_it++;
          else
            col_it++;
        }
        if (assign) {
          C.insert(i, j, sum);
          assign = false;
        }            
      }

    // sort
    C.sort();
  
  }

  //--------------------------------------------------------------------
  /*! \brief Trace of a sparse matrix.
    \param[in] A The matrix to calculate the trace of.
    \result The trace of the matrix.
  */
  template <class T>
  T trace(SparseMatrix<T>& A)
  { 
    typename SparseMatrix<T>::Row::iterator row_it, row_end;
    typename SparseMatrix<T>::Col::iterator col_it, col_end;
    T sum = 0;
  
    for (unsigned int k = 0; k < A.mr.size(); k++) {
      row_it = A.mr[k].begin();
      row_end = A.mr[k].end();
      col_it = A.mc[k].begin();
      col_end = A.mc[k].end();
      
      while ((row_it != row_end) && (col_it != col_end) &&
             (row_it->col_i <= k) && (col_it->row_i <= k)) {
        if ((row_it->col_i == k) && (col_it->row_i == k)) {
          sum += row_it->value;
          break;
        }
        else if(row_it->col_i < col_it->row_i)
          row_it++;
        else
          col_it++;
      }      
      
    }
  
    return sum;  
  }

  //--------------------------------------------------------------------
  /*! \brief Norm of a sparse matrix.
    \param[in] A The matrix to calculate the norm of.
    \result The norm of the matrix.
  */
  template <class T>
  T norm(SparseMatrix<T>& A)
  {  
    typename SparseMatrix<T>::Row::iterator row_it, row_end;
    T sum = 0;
  
    for (unsigned int i = 0; i < A.mr.size(); i++) {
      row_it = A.mr[i].begin();
      row_end = A.mr[i].end();
    
      while (row_it != row_end) sum += row_it++->value;
    }
  
    return sum;  
  }

  //--------------------------------------------------------------------
  /*! \brief Clustering coefficient of three adjacency matrices.

  Calculates the clustering coefficient from three adjacency matrices using
  sparse matrix multiplication.
  \param[in] J1 First adjacency matrix to consider.
  \param[in] J2 Second adjacency matrix to consider.
  \param[in] J3 Third adjacency matrix to consider.
  \result The clustering coefficient
  */
  template <class T>
  double clustering(SparseMatrix<T>& J1, SparseMatrix<T>& J2,
                    SparseMatrix<T>& J3)
  {
    double nom, denom;
    unsigned int n = J1.mr.size();
      
    SparseMatrix<T> JJ(n), JJJ(n);
      
    prod(J1, J2, JJ);
    prod(JJ, J3, JJJ);
    nom = static_cast<double>(trace(JJJ));
      
    JJ.clear();
      
    prod(J1, J2, JJ);
    denom = static_cast<double>(norm(JJ) - trace(JJ));
      
    return .5 * nom / denom ;
  }
  
}

//--------------------------------------------------------------------

  //--------------------------------------------------------------------
  /*! \brief Class for conversion adjacency list to matrices.

  Converts a boost-type adjacency list with two different edge types to two
  separate matrices.
  
  \ingroup helper_functions
  */
  template <typename Graph, typename Matrix>
  struct ListToMatrix {
    
    /*! \brief Constructor.
      \param[in] g_ g initialiser.
      \param[in] Jd_ Jd initialiser.
      \param[in] Ji_ Ji initialiser.
    */
    ListToMatrix(Graph& g_, Matrix& Jd_, Matrix& Ji_) :
      g(g_), Jd(Jd_), Ji(Ji_) {}
  
    /*! \brief Operator() to perform the actual conversion for one edge.
      \param[in] e The edge under consideration
    */
    typedef typename boost::graph_traits<Graph>::edge_descriptor ed;
    void operator() (const ed& e) {
    
      if (g[e].type == 0) {
        Jd.insert(source(e, g), target(e, g), 1);
        Jd.insert(target(e, g), source(e, g), 1);
      } else if (g[e].type == 1) {            
        Ji.insert(source(e, g), target(e, g), 1);
        Ji.insert(target(e, g), source(e, g), 1);
      } else { // should never get here
        std::cerr << "ERROR: corrupted edge\n";
      }
    
    }
  
    Graph& g; //!< The graph containing the adjacency list.
    Matrix& Jd; //!< The matrix storing the edges of type 0.
    Matrix& Ji; //!< The matrix storing the edges of type 1.
  };

//--------------------------------------------------------------------

namespace boost {

  //--------------------------------------------------------------------
  /*! \brief Calculate clustering coefficients.

  Calculates the eight possible clustering coefficients for a network of two
  edge types
  \param[in] g The graph containing the edges.
  \param[in] baseFileName The base name of the files to write the clustering
  coefficients and (if desired) adjacency matrices to.
  \param[in] writeJs Whether to save the adjacency matrices to disk.
  \ingroup graph_statistics
  */
  template <typename Graph>
  bool cluster_coeff(Graph& g, const std::string baseFileName,
                     const bool writeJs)    
  {
    // Adjacency matrix
    typedef SparseMatrix::SparseMatrix<unsigned int> Matrix;
    
    // size
    unsigned int n = boost::num_vertices(g);
  
    // creating the adjacecny matrices of size N
    Matrix Jd(n), Ji(n);
  
    // filling adjacecny matrices from g
    std::for_each(boost::edges(g).first, boost::edges(g).second,
                  ListToMatrix<Graph, Matrix>(g, Jd, Ji)); 
  
    // sort Matrices
    Jd.sort();
    Ji.sort();
  
    // print J's
    if (writeJs) {
      Jd.write_J(baseFileName + ".Jd");
      Ji.write_J(baseFileName + ".Ji");
    }
  
    // open cluster coeff file
    std::ofstream file;
    try {      
      file.open((baseFileName + ".cluster").c_str(), std::ios::out);
    }
    catch (std::exception &e) {
      std::cerr << "Unable to open cluster coeff file: " << e.what()
                << std::endl;
    }
  
    // write clustering coefficients
  
    file << "Cidi = " << SparseMatrix::clustering(Ji, Jd, Ji) << std::endl;
    file << "Cidd = " << SparseMatrix::clustering(Ji, Jd, Jd) << std::endl;
    file << "Cdii = " << SparseMatrix::clustering(Jd, Ji, Ji) << std::endl;
    file << "Cdid = " << SparseMatrix::clustering(Jd, Ji, Jd) << std::endl;
    file << "Cddi = " << SparseMatrix::clustering(Jd, Jd, Ji) << std::endl;
    file << "Cddd = " << SparseMatrix::clustering(Jd, Jd, Jd) << std::endl;
    file << "Ciii = " << SparseMatrix::clustering(Ji, Ji, Ji) << std::endl;
    file << "Ciid = " << SparseMatrix::clustering(Ji, Ji, Jd) << std::endl;
  
    // close file
    file.close();
  
    return 0;
  }

  //--------------------------------------------------------------------
  
} // namespace boost

  
