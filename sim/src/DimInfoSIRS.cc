/*! \file DimInfoSIRS.cc
  \brief Implementation of the Models::DimInfoSIRS class
*/
#include <iostream>
#include <fstream>

#include <math.h>

#include "DimInfoSIRS.hh"

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
Models::DimInfoSIRS::DimInfoSIRS(unsigned int v)
  : Model(v)
{
  /*************************************/
  // define vertex classes
  /************************************/
  // susceptible
  vertexStates.push_back(Label("S","01;32", 0, "", Label::rgbColour(0,0,255)));
  // infected 
  vertexStates.push_back(Label("I","01;31", 1, "", Label::rgbColour(255,0,0)));
  // recoveredn
  vertexStates.push_back(Label("R","01;34", 2, "", Label::rgbColour(0,255,0)));

  /*************************************/
  // define edge types
  /************************************/
  edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

  /*************************************/
  // define model parameters
  /************************************/
  model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate uninformed->uninformed")
    ("rho", po::value<double>(),
     "loss of quality in information")
    ("gamma", po::value<double>(),
     "recovery rate")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ("alpha", po::value<double>(),
     "information transmission rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("omega", po::value<double>(),
     "local information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ("threshold", po::value<double>(),
     "threshold below which not to consider information transmission")
    ("dim_alpha", po::value<double>()->default_value(0.),
     "whether alpha diminishes or not")
    ;

  /*************************************/
  // assign model parameters to variables
  /************************************/
  rates.insert(std::make_pair("beta", &beta));
  rates.insert(std::make_pair("gamma", &gamma));
  rates.insert(std::make_pair("delta", &delta));
  rates.insert(std::make_pair("alpha", &alpha));
  rates.insert(std::make_pair("nu", &nu));
  rates.insert(std::make_pair("omega", &omega));
  rates.insert(std::make_pair("lambda", &lambda));

  params.insert(std::make_pair("rho", &rho));
  params.insert(std::make_pair("threshold", &threshold));
  params.insert(std::make_pair("dim_alpha", &dim_alpha));
}

//----------------------------------------------------------
void Models::DimInfoSIRS::Init(po::variables_map& vm)
{
  Model::Init(vm);
}

//----------------------------------------------------------
unsigned int Models::DimInfoSIRS::getNodeEvents(eventList& events,
                                                State state,
                                                unsigned int nb) const
{
   unsigned int rateSum(0);

   if (getDisease(state) == Recovered) {
     // loss of immunity
      Event immunityLoss;
      immunityLoss.rate = delta;
      immunityLoss.newState = State(Susceptible, -1);
      immunityLoss.nb = nb;
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
        if (verbose >= 2) {
          std::cout << "Adding loss of immunity event with rate " 
	            << immunityLoss.rate/1e+4 << std::endl;
        }
      }
   } else
   if (getDisease(state) == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma;
      recovery.newState = State(Recovered, -1);
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (verbose >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate/1e+4
                    << std::endl;
        }
      }
      if (state.detail < 1) {
        // local information generation
        Event localInfo;
        localInfo.rate = omega;
        localInfo.newState = State(state.base, 1);
        localInfo.nb = nb;
        if (localInfo.rate > 0) {
          events.push_back(localInfo);
          rateSum += localInfo.rate;
          if (verbose >= 2) {
            std::cout << "Adding local information event with rate " 
	              << localInfo.rate/1e+4 << std::endl;
          }
        }
      }
   }
   // information loss
   if (state.detail > 0) {
      Event infoLoss;
      infoLoss.rate = lambda;
      double detailUpdate = state.detail * rho;
      infoLoss.newState = State(state.base, detailUpdate);
      infoLoss.nb = nb;
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
        if (verbose >= 2) {
          std::cout << "Adding information loss event with rate " 
	            << infoLoss.rate/1e+4 << std::endl;
        }
      }
   }

   return rateSum;
}

//----------------------------------------------------------
unsigned int Models::DimInfoSIRS::getEdgeEvents(eventList& events,
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
         infection.rate =
           static_cast<unsigned int>((1 - state.detail) * beta + .5);
         infection.newState = State(Infected, -1);
         infection.nb = nb;
         infection.et = edge;
         if (infection.rate > 0) {
           events.push_back(infection);
           rateSum += infection.rate;
           if (verbose >= 2) {
             std::cout << "Adding infection event with rate "
                       << infection.rate << std::endl;
           }
         }
      }
   } else if (edge == Information) {
      // information transmission
      if (state.detail < nbState.detail * rho && nbState.detail > threshold) {
         Event infoTransmission;
         if (dim_alpha) {
           infoTransmission.rate =
             static_cast<unsigned int>(nbState.detail * alpha + .5);
         } else {
           infoTransmission.rate = alpha;
         }
         infoTransmission.newState = State(state.base, nbState.detail * rho);
         infoTransmission.nb = nb;
         infoTransmission.et = edge;
         if (infoTransmission.rate > 0) {
           events.push_back(infoTransmission);
           rateSum += infoTransmission.rate;
           if (verbose >= 2) {
             std::cout << "Adding information transmission event with rate " 
                       << infoTransmission.rate/1e+4 << std::endl;
           }
         }
      }
      // information generation
      if (state.detail < 1 && getDisease(state) == Susceptible && 
          getDisease(nbState) == Infected) {
         Event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState = State(state.base, 1);
         infoGeneration.et = edge;
         infoGeneration.nb = nb;
         if (infoGeneration.rate > 0) {
           events.push_back(infoGeneration);
           rateSum += infoGeneration.rate;
           if (verbose >= 2) {
             std::cout << "Adding information generation event with rate " 
                       << infoGeneration.rate/1e+4 << std::endl;
           }
         }
      }
   }

   return rateSum;
}
