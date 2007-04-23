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

namespace boost {
  
  
  //--------------------------------------------------------------------
  
  template <class T, class C>
  struct Entry {
    
    template <class P, class Q>
    friend std::ostream& operator<< (std::ostream& os, const Entry<P, Q>& entry);
  
    template <class P, class Q>
    friend bool operator< (const Entry<P, Q>& lhs, const Entry<P, Q>& rhs);
  
    Entry(unsigned int r, unsigned int c, T v) : row_i(r), col_i(c), value(v) {}
    ~Entry() {}
  
    unsigned int row_i, col_i;
 
    T value;
  
  };

  //--------------------------------------------------------------------
  // operator<<

  template <class T, class C>
  std::ostream& operator<< (std::ostream& os, const Entry<T, C>& entry)
  {
    os << entry.row_i << "  " << entry.col_i << "  " << entry.value << std::endl;
    return os;
  }

  //--------------------------------------------------------------------
  // operator<

  template <class T, class C> 
  bool operator< (const Entry<T, C>& lhs, const Entry<T, C>& rhs)
  {
    return C::compare(lhs.row_i, lhs.col_i, rhs.row_i, rhs.col_i);
  }

  //--------------------------------------------------------------------
  // operator< for row major

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
  // operator< for column major
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

  template <typename Graph, typename Matrix>
  struct ListToMatrix {
  
    ListToMatrix(Graph& g_, Matrix& Jd_, Matrix& Ji_) : g(g_), Jd(Jd_), Ji(Ji_) {}
  
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
  
    Graph& g;
    Matrix& Jd;
    Matrix& Ji;
  };

  //--------------------------------------------------------------------

  template <class T>
  struct SparseMatrix {
  
    typedef Entry<T, RowMajor> Mij;
    typedef Entry<T, ColMajor> Mji;
    typedef std::list<Mij> Row;
    typedef std::list<Mji> Col;
    typedef std::vector<Row> RowMajorMatrix;
    typedef std::vector<Col> ColMajorMatrix;
  
    SparseMatrix(unsigned int size) {
      mr.resize(size);
      mc.resize(size);
    }
  
    ~SparseMatrix() {}
  
    void insert(unsigned int i, unsigned int j, T val) {
      mr[i].push_back(Mij(i, j, val));
      mc[j].push_back(Mji(i, j, val));
    }
  
    void write_J(std::string fileName);
  
    void sort() {
      for (unsigned int i = 0; i < mr.size(); ++i) mr[i].sort();    
      for (unsigned int j = 0; j < mc.size(); ++j) mc[j].sort();
    }

    void clear() {
      for (unsigned int i = 0; i < mr.size(); ++i) mr[i].clear();    
      for (unsigned int j = 0; j < mc.size(); ++j) mc[j].clear();
    }
    
    RowMajorMatrix mr;
    ColMajorMatrix mc;
    
  };
  
  //--------------------------------------------------------------------
  // write_J
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
  // general sparse matrix multiplication

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
  // general sparse matrix trace

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
  // general sparse matrix norm

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
  // clustering coefficient
    
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
  

  //--------------------------------------------------------------------

  template <typename Graph>
  bool cluster_coeff(Graph& g, const std::string baseFileName,
                     const bool writeJs)    
  {
    // Adjacency matrix
    typedef boost::SparseMatrix<unsigned int> Matrix;
    
    // size
    unsigned int n = boost::num_vertices(g);
  
    // creating the adjacecny matrices of size N
    Matrix Jd(n), Ji(n);
  
    // filling adjacecny matrices from g
    std::for_each(boost::edges(g).first, boost::edges(g).second,
                  boost::ListToMatrix<Graph, Matrix>(g, Jd, Ji)); 
  
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
  
    file << "Cidi = " << boost::clustering(Ji, Jd, Ji) << std::endl;
    file << "Cidd = " << boost::clustering(Ji, Jd, Jd) << std::endl;
    file << "Cdii = " << boost::clustering(Jd, Ji, Ji) << std::endl;
    file << "Cdid = " << boost::clustering(Jd, Ji, Jd) << std::endl;
    file << "Cddi = " << boost::clustering(Jd, Jd, Ji) << std::endl;
    file << "Cddd = " << boost::clustering(Jd, Jd, Jd) << std::endl;
    file << "Ciii = " << boost::clustering(Ji, Ji, Ji) << std::endl;
    file << "Ciid = " << boost::clustering(Ji, Ji, Jd) << std::endl;
  
    // close file
    file.close();
  
    return 0;
  }

  //--------------------------------------------------------------------
  
} // namespace boost

  
