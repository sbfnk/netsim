/*! \file Bin.cc
  \brief Implementation of the Tree::Bin class
*/

#include <iostream>

#include "Bin.hh"

namespace Tree {

  //----------------------------------------------------------
  /*!
    \brief Generate a bin entry to accomodate a new Leaf.

    This is called by the Tree for the parent of the most recently added Leaf if
    a new Leaf is to be added. The Bin checks whether it can still accomodate a
    child. If yes, it adds the Leaf to its children and sets itself as a
    parent. If, on the other hand, it cannot hold any more children (because it
    already has 2), it creates a new Bin which becomes the parent of the new
    leaf, and tries to generate a Bin entry for that new (parent) Bin with its
    own parent by calling generateBinEntry recursively. This proceeds up the
    tree structure until a vacant Bin is found or the top Bin is reached, in
    which case a new top Bin is created and new branch added to tree.
  
    \param[in] linkToBin The bin to be added to the tree
    \return 0 if the new Leaf could be accomodated without generating a new
    level, a pointer to the top bin if a new top bin had to be created.
  */
  Bin* Bin::generateBinEntry(Bin& linkToBin)
  {
    if (getNChildren() < containerSize) {
      // we have an empty place, so we add linkToBin to our children
      // and we set ourselves as parent for linkToBin
      addChild(linkToBin);
      linkToBin.setParent(this);
      return 0;
    } else {
      // we have no empty place, so we have to create a new bin
      Bin* newBin(new Bin());
      // the new bin is the new parent for linkToBin
      newBin->addChild(linkToBin);
      linkToBin.setParent(newBin);
      if (parent) {
        // we have not been top bin, so we generate an entry at one
        // level higher for the new bin
        return parent->generateBinEntry(*newBin);
      } else {
        // we have been top bin, so we need a new one
        Bin* topBin = new Bin();
        // we set ourselves in the first slot of the new top bin
        topBin->addChild(*this);
        setParent(topBin);
        // and the new bin in the second one
        topBin->addChild(*newBin);
        newBin->setParent(topBin);
        // we return the new topBin for the GillespieGraph to
        // set its own topBin variable
        return topBin;
      }
    }
  }

  //----------------------------------------------------------
  /*!
    \brief Add a child.

    Adds a new child to the Bin and update its ratesum.
    \param[in] newChild A reference to the new child.
    \return The index of the new child in the children vector if succeeded, -1
    if failed 
  */
  short Bin::addChild(Bin& newChild)
  {
    if (children.size() < containerSize) {
      children.push_back(&newChild);
      if (newChild.getRateSum() != 0) {
        // if lower has a rate, update the rates in our tree branch
        updateRateSum(newChild.getRateSum());
      }
      return children.size() - 1;
    } else {
      return -1;
    }
  }

  //----------------------------------------------------------
  /*!
    \brief Update the rate sum.

    This updates both the rate sum of the parent and the rate sum of the Bin by
    adding the ifference between the old and new rate sum, passed as a
    parameter. 

    \param[in] rate The difference between old and new rate sum
  */
  bool Bin::updateRateSum(int rate)
  {

    bool result = true;
    
    if (parent) {
      // we are not the top level bin, so we update the higher level bin
      result &= parent->updateRateSum(rate);
    }
    // update own rate sum
    rateSum += rate;
    if (rate > 0 && rateSum < rate) {
      result = false;
    }

    return result;
  }
      
  //----------------------------------------------------------
  /*!
    \brief Pick a child according to a given random number.

    Selects a bin from its two children using a number between 0 and its ratesum
    and recusively calls pickChild for that bin, stripping the random number of
    the rates corresponding to the branches not traversed.

    \param[in] randNo A number between 0 and the ratesum of the bin, determining
    which Child should be picked.
    \return A pointer the leaf ultimately selected.
  */
  Bin* Bin::pickChild(unsigned int& randNo)
  {
    // check if we are at the lowest level already
    if (children.size() > 0) {
      // we are not yet on the lowest level, so we further go down the tree
      std::vector<Bin*>::iterator it = children.begin();
      unsigned int randSum(0);
      // get the rate represented by the first slot
      unsigned int tempSum((*it)->getRateSum());
      // find correct bin. This could be a simple comparison for the binary
      // case, but like this it can account for any container size
      while (it != children.end() && randNo > tempSum) {
        randSum = tempSum;
        it++;
        tempSum += (*it)->getRateSum();
      }
      // pick a bin on the lower level, subtracting the rate representing
      // the branches we are not traversing
      randNo -= randSum;
      return (*it)->pickChild(randNo);
    } else {
      // we are at the lowest level and should be a leaf
      return this;
    }
  }

}
