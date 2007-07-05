#ifndef BIN_HH
#define BIN_HH

#include <vector>
#include <iostream>

const unsigned short containerSize = 2;
template<typename T>
class Leaf;

class Bin
{
  
public:
  Bin(): parent(0), rateSum(0.) {;}
  virtual ~Bin()
  { for (unsigned int i = 0; i < children.size(); i++)
      if (children[i]) delete children[i];
  }
  
  void setParent(Bin* const newParent) { parent = newParent; }
  
  Bin* generateBinEntry(Bin& linkToBin);
  Bin* getParent() const { return parent; }
  double getRateSum() const { return rateSum; }
  short getNChildren() const { return children.size(); }
  
  short addChild(Bin& newChild);
  
  void updateRateSum(double rate);
  
  Bin* pickBin(double& randNo);
  void setRateSum(double rate);
  
private:
  
  Bin* parent;
  double rateSum;
  
  std::vector<Bin*> children;
};

template<typename T>
class Leaf :
   public Bin
{
   public:
      Leaf(T newItem): item(newItem) {;}
      T* getItem() { return &item; }
   private:
      T item;
};

#endif
