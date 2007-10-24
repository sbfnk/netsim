/*! \file SingleSIRS.cc
  \brief Implementation of the Models::SingleSIRS class
*/
#include <iostream>
#include <fstream>

#include "SingleSIRS.hh"

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
Models::SingleSIRS::SingleSIRS(unsigned int v)
  : Model(v)
{
  /*************************************/
  // define vertex classes
  /************************************/
  // susceptible
  vertexStates.push_back(Label("S","00;32", 0, "fillcolor=\"royalblue4\""));
  // infected
  vertexStates.push_back(Label("I","00;31", 1, "fillcolor=\"red\""));
  // recovered
  vertexStates.push_back(Label("R","00;34", 2, "fillcolor=\"green\""));

  /*************************************/
  // define edge types
  /************************************/
  edgeTypes.push_back(Label("0", "", 0, "style=\"solid\""));
  edgeTypes.push_back(Label("1", "", 1, "style=\"dashed\""));

  /*************************************/
  // define model parameters
  /************************************/
  model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate")
    ("gamma", po::value<double>(),
     "recovery rate")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ;

  /*************************************/
  // assign model parameters to variables
  /************************************/
  params.insert(std::make_pair("beta", &beta));
  params.insert(std::make_pair("gamma", &gamma));
  params.insert(std::make_pair("delta", &delta));
}

//----------------------------------------------------------
double Models::SingleSIRS::getNodeEvents(eventList& events,
                                       State state,
                                       unsigned int nb) const
{
   double rateSum(.0);

   if (state.base == Recovered) {
     // loss of immunity
      Event immunityLoss;
      immunityLoss.rate = delta;
      immunityLoss.newState = State(Susceptible);
      immunityLoss.nb = nb;
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
        if (verbose >= 2) {
          std::cout << "Adding loss of immunity event with rate " 
	            << immunityLoss.rate << std::endl;
        }
      }
   } else if (state.base == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma;
      recovery.newState = State(Recovered);
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (verbose >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate
                    << std::endl;
        }
      }
   }
   return rateSum;
}

//----------------------------------------------------------
double Models::SingleSIRS::getEdgeEvents(eventList& events,
                                       State state,
                                       unsigned int edge,
                                       State nbState,
                                       unsigned int nb) const
{
   double rateSum(.0);
   // infection
   if (state.base == Susceptible &&
       nbState.base == Infected) {
     Event infection;
     infection.rate = beta;
     infection.newState = State(Infected);
     infection.nb = nb;
     infection.et = edge;
     if (infection.rate > 0) {
       events.push_back(infection);
       rateSum += infection.rate;
       if (verbose >= 2) {
         std::cout << "Adding infection event with rate " << infection.rate
                   << std::endl;
       }
     }
   }

   return rateSum;
}
