/*! \file Tree.hh
  \brief The Tree class.
*/
#ifndef TREE_HH
#define TREE_HH

#include <vector>
#include "Bin.hh"

#include <boost/graph/adjacency_list.hpp>

/*! \brief Classes related to tree structure.
  \ingroup gillespie_simulator
*/
namespace Tree {

  /*! \brief Structure for storing rates in a tree for quick access.
    
  This is a tree containing of a vector of leaves (lowest level items) holding a
  variable of arbitrary type and a number (here: rate), connected in two to
  higher level bins which have the sum of the rates of their children as their
  own sum. This allows to choose events quickly in the Gillespie algorithm by
  starting from the top bin and going down to the lowest level Leaf
  corresponding to a given random number.
  
  */
  template <class T>
  class Tree
  {
    
  public:
    //! Constructor.
    Tree(): lastLeaf(0), topBin(0) {;}
    /*! \brief Destructor.
      
    Deletes the top bin, which will in turn cause the whole tree to be deleted
    in a cascade of deletions.
    */
    ~Tree() { clear(); }
    void clear()
    { if (topBin) delete topBin; leaves.clear(); lastLeaf = 0; topBin = 0; }
    
    T* pickRandomElement(unsigned int& randNo);

    //! Accessor for the topBin variable
    Bin* getTopBin() { return topBin; }
    
    void addLeaf(double rateSum, T leafItem);

    //! Accessor for the leaves variable
    std::vector<Leaf<T>*>& getLeaves() { return leaves; }
    
  private:
    Leaf<T>* lastLeaf; //!< The latest leaf added of the tree
    Bin* topBin; //!< top Bin of the tree
    std::vector<Leaf<T>*> leaves; //!< A map from numbers (vertex id) to leaves
    
  };
  
  //----------------------------------------------------------
  /*! \brief Pick a random leaf from the tree.

  Picks a random leaf from the tree, based on a random number between 0 and
  1. It starts with the top Bin and propagates down the tree to find the
  corresponding Leaf, each Leaf weighed by its rateSum.
  \param[in] randNo A random between 0 and 1 to base the choice of Leaf on.
  \return A pointer to the item contained in the chosen Leaf.
  */
  template <class T>  
  T* Tree<T>::pickRandomElement(unsigned int& randNo)
  {
    // check if we have a topBin (if not we either have no vertex
    // or something else is seriously wrong)
    if (topBin) {
      // get the sum of the two slots for the topBin, corresponding
      // to the total sum of rates of the whole tree
      unsigned int rateSum = topBin->getRateSum();
      if (rateSum > 0) {
        // call pickBin of the topBin, which will traverse down the tree
        // using randNo, store the vertex eventually picked in eventBin
        // have leafBin us randNo (a random number in [0, eventBin->rateSum))
        // to pick its own event
        Leaf<T>* pickedLeaf =
          dynamic_cast<Leaf<T>*>(topBin->pickChild(randNo));
        if (pickedLeaf) {
          return pickedLeaf->getItem();
        } else {
          std::cout << "picked bin is not a leaf, very strange..."
                    << std::endl;
        }
      } else {
        std::cout << "nothing can happen" << std::endl;
      }
    } else {
      std::cout << "Could not find top bin." << std::endl;
    }
    // we should never end up here
    return 0;
  }
  
  //----------------------------------------------------------
  /*! \brief Add an item to the Tree.

  Creates the first item in the Tree or, if there are already items, expands the
  Tree to accomodate the new item.
  
  \param[in] rateSum The rate corresponding to the new item.
  \param[in] leafItem The item to add to the Tree.
   */
  template <class T>
  void Tree<T>::addLeaf(double rateSum, T leafItem)
  {
    // create the new leaf with given elementNumber (position in the
    // vector) 
    Leaf<T>* newLeaf =
      new Leaf<T>(leafItem);
    if (lastLeaf) {
      // we are not creating the first vertex, so we try to an empty slot
      // in the tree
    
      // call generateBinEntry for the Bin the lowest Bin is sitting in.
      // This traverses the tree and accomodates gVertex in an empty slot
      // and returns 0 or, if it cannot find an empty slot, creates a new
      // top bin and a branch stretching down to gVertex, returning a
      // pointer to the new top bin
      Bin* newTop =
        lastLeaf->getParent()->generateBinEntry(*newLeaf);
      // if generateBinEntry returned non-zero we have a new top bin
      if (newTop) {
        topBin = newTop;
      }
    } else {
      // no vertex exists yet, so create the first one
      Bin* newBin = new Bin();
      // add gVertex to the new bin and
      // set the new bin as top-level bin for gVertex
      newBin->addChild(*newLeaf);
      newLeaf->setParent(newBin);
      // the first bin is also the first top bin
      topBin = newBin;
    }
    // set the new leaf to be the new last leaf
    lastLeaf = newLeaf;
    getLeaves().push_back(newLeaf);
    return;
  }

  //----------------------------------------------------------
  /*! \brief Generate a Tree from a graph.

  Loops over all vertices in the graph and creates a Tree assigning each vertex
  a place in it. This allows to later choose events quickly based on the random
  number generated in the GillespieSimulator.

  \param[out] t An empty Tree variable to create the leaves in.
  \param[in] g The graph containing the vertices to be accomodated in the Tree.
  \param[in] rates A property map containing the rates corresponding to the
  vertices. 
  \param[in] indices A property map containing the vertex indices.
   */
  template <class T, class Graph, class ratesPropertyMap,
            class indexPropertyMap>
  void generateTree(Tree<T>& t, const Graph& g,
                    ratesPropertyMap rates, indexPropertyMap indices)
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;

    t.clear();
   
    vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g);
         vi != vi_end; ++vi) {
      // add item to tree
      t.addLeaf(get(rates, *vi), get(indices, *vi));
      // update rate sum
      t.getLeaves().back()->updateRateSum(get(rates, *vi));
    }
  }

}

#endif
