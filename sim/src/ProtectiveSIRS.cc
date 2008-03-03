/*! \file ProtectiveSIRS.cc
  \brief Implementation of the Models:: class
*/
#include <iostream>
#include <fstream>

#include "ProtectiveSIRS.hh"

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
Models::ProtectiveSIRS::ProtectiveSIRS(unsigned int v)
  : Model(v)
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

  rates.insert(std::make_pair("beta", &beta));
  rates.insert(std::make_pair("gamma", &gamma));
  rates.insert(std::make_pair("delta", &delta));
  rates.insert(std::make_pair("nu", &nu));
  rates.insert(std::make_pair("lambda", &lambda));
  params.insert(std::make_pair("sigma", &sigma));
}

//----------------------------------------------------------
unsigned int Models::ProtectiveSIRS::getNodeEvents(eventList& events,
                                                   State state,
                                                   unsigned int nb) const
{
   unsigned int rateSum(0);

   // loss of immunity
   if (getDisease(state) == Recovered) {
      Event immunityLoss;
      immunityLoss.rate = delta;
      immunityLoss.newState = State(getState(Susceptible, getInfo(state)));
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
      }
   } else
   if (getDisease(state) == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma;
      recovery.newState = State(getState(Recovered, getInfo(state)));
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
      }
   }
   // information loss
   if (getInfo(state) == Informed) {
      Event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState = State(getState(getDisease(state), Uninformed));
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
      }
   }

   return rateSum;
}

//----------------------------------------------------------
unsigned int Models::ProtectiveSIRS::getEdgeEvents(eventList& events,
                                                   State state,
                                                   unsigned int edge,
                                                   State nbState,
                                                   unsigned int nb) const
{
   unsigned int rateSum(0);
   if (edge == Disease) {
      // infection
      if (getDisease(state) == Susceptible &&
          getDisease(nbState) == Infected) {
         Event infection;
         infection.rate = beta;
         if (getInfo(state) == 1) {
           infection.rate =
             static_cast<unsigned int>(infection.rate * sigma + .5);
         }
         infection.newState = State(getState(Infected, getInfo(state)));
         if (infection.rate > 0) {
           events.push_back(infection);
           rateSum += infection.rate;
         }
      }
   } else if (edge == Information) {
     if (getDisease(state) == Susceptible &&
         getInfo(state) == Uninformed && getDisease(nbState) == Infected) {
      // information generation
         Event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState = State(getState(getDisease(state), Informed));
         if (infoGeneration.rate > 0) {
           events.push_back(infoGeneration);
           rateSum += infoGeneration.rate;
         }
      }
   }

   return rateSum;
}
