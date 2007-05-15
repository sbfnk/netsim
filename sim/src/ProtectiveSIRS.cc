/*******************************************************************/
//
// Class ProtectiveSIRS
// --------------------
//
// This is the class implementing the SIRS model, i.e. the parameters
// and the dependance of the rates of various processes on node states
// and edge properties
//
/******************************************************************/

#include <iostream>
#include <fstream>

#include "ProtectiveSIRS.hh"

/******************************************************************/
// ProtectiveSIRS constructor
/******************************************************************/
ProtectiveSIRS::ProtectiveSIRS()
{
  // susceptible uninformed
  vertexStates.push_back(Label("S-","00;32", 0, "fillcolor=\"royalblue4\""));
  // susceptible informed
  vertexStates.push_back(Label("S+","01;32", 1, "fillcolor=\"royalblue\""));
  // infected 
  vertexStates.push_back(Label("I","00;31", 2, "fillcolor=\"red\""));
  // recovered
  vertexStates.push_back(Label("R","00;34", 3, "fillcolor=\"green\""));

  edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

  model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate")
    ("sigma", po::value<double>(),
     "reduction of transmission rate due to information")
    ("gamma", po::value<double>(),
     "recovery rate")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ;

  params.insert(std::make_pair("beta", &beta));
  params.insert(std::make_pair("sigma", &sigma));
  params.insert(std::make_pair("gamma", &gamma));
  params.insert(std::make_pair("delta", &delta));
  params.insert(std::make_pair("nu", &nu));
  params.insert(std::make_pair("lambda", &lambda));
}

/******************************************************************/
// Vertex destructor
/******************************************************************/
ProtectiveSIRS::~ProtectiveSIRS()
{}

/******************************************************************/
// ProtectiveSIRS::getNodeEvents
// get the events that can happen for a given state of a node. Stores the
// events in the events list and returns the sum of their rates
/******************************************************************/
double ProtectiveSIRS::getNodeEvents(eventList& events,unsigned int state) const
{
   double rateSum(.0);

   // loss of immunity
   if (getDisease(state) == Recovered) {
      event immunityLoss;
      immunityLoss.rate = delta;
      immunityLoss.newState = getState(Susceptible, getInfo(state));
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
      }
   } else
   if (getDisease(state) == Infected) {
      // recovery
      event recovery;
      recovery.rate = gamma;
      recovery.newState = getState(Recovered, getInfo(state));
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
      }
   }
   // information loss
   if (getInfo(state) == Informed) {
      event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState = getState(getDisease(state), Uninformed);
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
      }
   }

   return rateSum;
}

/******************************************************************/
// ProtectiveSIRS::getNodeEvents
// get the events that can happen for a given edge type and the states of the
// neighbouring nodes. Stores the events in the events list and returns the sum
// of their rates 
/******************************************************************/
double ProtectiveSIRS::getEdgeEvents(eventList& events,
                            unsigned int state, unsigned int edge,
                            unsigned int nbState) const
{
   double rateSum(.0);
   if (edge == Disease) {
      // infection
      if (getDisease(state) == Susceptible &&
          getDisease(nbState) == Infected) {
         event infection;
         infection.rate = beta;
         if (getInfo(state) == 1) infection.rate *= sigma;
         infection.newState = getState(Infected, getInfo(state));
         if (infection.rate > 0) {
           events.push_back(infection);
           rateSum += infection.rate;
         }
      }
   } else if (edge == Information) {
     if (getDisease(state) == Susceptible &&
         getInfo(state) == Uninformed && getDisease(nbState) == Infected) {
      // information generation
         event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState = getState(getDisease(state), Informed);
         if (infoGeneration.rate > 0) {
           events.push_back(infoGeneration);
           rateSum += infoGeneration.rate;
         }
      }
   }

   return rateSum;
}
