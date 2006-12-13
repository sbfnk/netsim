/******************************************************************/
// Tree.hh
// contains the Tree class, used by the Gillespie class to store
// event rates
/******************************************************************/
#ifndef TREE_HH
#define TREE_HH

#include <vector>
#include "Bin.hh"

#include <boost/graph/adjacency_list.hpp>

template <class T>
class Tree
{
   public:
      Tree(): lastLeaf(0), topBin(0) {;}

      T& pickRandomElement(double& randNo)
      {
         // check if we have a topBin (if not we either have no vertex
         // or something else is seriously wrong)
         if (topBin) {
            // get the sum of the two slots for the topBin, corresponding
            // to the total sum of rates of the whole tree
            double rateSum = topBin->getRateSum();
            if (rateSum > 0.) {
               // call pickBin of the topBin, which will traverse down the tree
               // using randNo, store the vertex eventually picked in eventBin
               // and make randNo a random number in [0, eventBin->rateSum),
               // which the leafBin can use to pick its own event
               randNo *= rateSum;
               // **** debugging printout ****
               std::cout << "rateSum: " << rateSum << " randNo: " << randNo
                         << std::endl;

               Leaf<T>* pickedLeaf =
                  dynamic_cast<Leaf<T>*>(topBin->pickBin(randNo));
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
         return *(new T);
      }

      Bin* getTopBin() { return topBin; }
      
      Leaf<T>* lastLeaf; // latest leaf added of the tree

      // reverse map from elements to leaves
      std::vector<Leaf<T>*> leaves;

      void addLeaf(double rateSum, T leafItem)
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
         leaves.push_back(newLeaf);
         return;
      }
   private:
      Bin* topBin; // top Bin of the tree
};

template <class T, class Graph, class ratesPropertyMap, class indexPropertyMap>
void generateTree(Tree<T>& t, const Graph& g,
                  ratesPropertyMap rates, indexPropertyMap indices)
{
   typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iterator;
   vertex_iterator vi, vi_end;
   for (tie(vi, vi_end) = vertices(g);
        vi != vi_end; ++vi) {
      t.addLeaf(get(rates, *vi), get(indices, *vi));
      t.leaves.back()->updateRateSum(get(rates, *vi));
   }
}

#endif
