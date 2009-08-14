/*! \file Bin.hh
  \brief The Tree::Bin class.
*/

#ifndef BIN_HH
#define BIN_HH

#include <vector>
#include <iostream>

namespace Tree{

  template<typename T>
  class Leaf;

  /*! \brief A bin which can be stored in a (binary) tree.
    
  This is used for the Bin containers needed for the tree
  structure in the GillespieAlgorithm. Each Bin has a pointer to the higher
  level Bin (parent) (except for the top Bin), and a vector of pointers to the
  lower level Bins occupying its own slots (children). At the lowest level,
  these two pointers are empty, and the Bins are actually Leaves (which inherit
  from Bin). Each Bin has a member rateSum, which is the sum of the two rateSum
  variables of the children.

  */
  class Bin
  {
  
  public:

    //! Constructor
    Bin(): parent(0), rateSum(0), containerSize(2) {;}

    /*! \brief Destructor.
      
    Deletes the two children, which will in turn delete their children.
    */
    virtual ~Bin()
    { for (unsigned int i = 0; i < children.size(); i++) delete children[i]; }
  
    //! Mutator of the parent variable.
    void setParent(Bin* const newParent) { parent = newParent; }
  
    Bin* generateBinEntry(Bin& linkToBin);

    //! Accessor for the parent variable.
    Bin* getParent() const { return parent; } 
    //! Accessor for the rateSum variable.
    unsigned int getRateSum() const { return rateSum; } 
    //! Accessor for the children variable.
    short getNChildren() const { return children.size(); } 
  
    short addChild(Bin& newChild);
    
    bool updateRateSum(int rate);
  
    Bin* pickChild(unsigned int& randNo);
  
  private:
  
    Bin* parent; //!< Pointer to the parent Bin
    unsigned rateSum; //!< The sum of all rates contained in the Bin
  
    std::vector<Bin*> children; //!< Vector of pointers to the children
    const unsigned short containerSize; //!< Size of the Bin (usually 2)
  };


  //----------------------------------------------------------
  /*! \brief A leaf, i.e. a bin containing an item of arbitrary type
    
  A member of this class is the lowest end in a tree, i.e. it has no
  children. Instead, it contains an item of arbitrary type (in the case of the
  GillespieSimulator, an unsigned int corresponding to the vertex index), and
  its rateSum corresponds to the sum of the rates of the events that can happen
  to the vertex.
  */
  template<typename T>
  class Leaf :
    public Bin
  {

  public:
    /*! \brief Constructor.
      
    \param[in] newItem The item to be held by the Leaf.
    */
    Leaf(T newItem): item(newItem) {;}

    //! Get a pointer to the item held by the leaf
    T* getItem() { return &item; }
  private:
    T item; //!< The item held by the leaf
  
  };
  
}

#endif
