/*! \file Model.hh
  \brief The Model class.
*/
#ifndef STATRECORDER_HH
#define STATRECORDER_HH

#include <list>
#include <vector>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>
#include <Vertex.hh>

//! \addtogroup models Models

namespace po = boost::program_options;

std::string generateFileName(std::string nameBase, unsigned int id,
                             std::string ext = "");

// functoid base class
template <class Graph>
class Funct {
public:
  virtual void doit(const Graph& g, std::string s,
                    double time, unsigned int count) = 0;
  virtual ~Funct() = 0;
};

template <class Graph>
inline Funct<Graph>::~Funct() { }

template <typename Graph>
class StatRecorder 
{
public:
  
  StatRecorder(Funct<Graph>* s, double f) :
    statFunc(s), outputFreq(f), nextStep(f),
    seqNum(0), outputDir("")
  {;}

  ~StatRecorder()
  { delete statFunc; }
  
  void reset(std::string dir)
  {
    nextStep = outputFreq;
    outputDir = dir;
    seqNum = 0;  
  }

  void update(const Graph& g, double time, bool force = false) 
  {
    if (force || (nextStep > 0. && time > nextStep) || (nextStep < 0.)) {
      statFunc->doit(g, outputDir, time, seqNum);
      if (nextStep > 0. && time > nextStep) { nextStep += outputFreq; }
      ++seqNum;
    }
  }

  const Funct<Graph>* getStatFunc() const {return statFunc;}
  
private:
  Funct<Graph>* statFunc;
  double outputFreq;

  double nextStep;
  unsigned int seqNum;

  std::string outputDir;
};

#endif
