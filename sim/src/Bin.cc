/*******************************************************************/
//
// Class Bin
// --------------------
//
// This is used for the bin containers needed for the tree
// structure in the GillespieAlgorithm. Each bin has a pointer parent
// to the higher level bin it is sitting in, and a vector of two pointers
// lowerBin, pointing to the lower level bins occupying its own slots. At the
// lowest level, these two pointers are empty, and the Bins are actually
// GillespieVertices (which inherit from Bin). Each bin also has a variable
// rateSum containing the rate it represents 
//
/******************************************************************/

#include <iostream>

#include "Bin.hh"

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

void Bin::updateRateSum(double rate)
{
   if (parent) {
      // we are not the top level bin, so we update the higher level bin
      parent->updateRateSum(rate);
   }
   // update own rate sum
   rateSum += rate;
}
      
/******************************************************************/
// Bin::pickBin
// selects a bin from its two children using randNo and returns
// the leaf ultimately selected, stripping randNo of the rates representing
// the branches not treversed.
/******************************************************************/
Bin* Bin::pickBin(double& randNo)
{
   // check if we are at the lowest level already
   if (children.size() > 0) {
      // we are not yet on the lowest level, so we further go down the tree
      std::vector<Bin*>::iterator it = children.begin();
      double randSum(.0);
      // get the rate represented by the first slot
      double tempSum((*it)->getRateSum());
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
      return (*it)->pickBin(randNo);
   } else {
      // we are at the lowest level and should be a leaf
      return this;
   }
}

/******************************************************************/
// Bin::setRateSum
// update our rate and the rates in the branch we are sitting in
// rate represents the new rate, so it should be >= 0
/******************************************************************/
void Bin::setRateSum(double rate)
{
   if (parent) {
      // we are not the top level bin, so we update the higher level bin
      parent->updateRateSum(rate - rateSum);
   }
   // set own rate sum
   rateSum = rate;
}
