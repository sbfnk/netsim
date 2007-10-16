/*! \file InfoSIRS.cc
  \brief Implementation of the Models::InfoSIRS class
*/
#include <iostream>
#include <fstream>

#include "InfoSIRS.hh"

//----------------------------------------------------------
/*! \brief Constructor.

\sa Model::Model
*/
Models::InfoSIRS::InfoSIRS(unsigned int v)
  : Model(v)
{
  /*************************************/
  // define vertex classes
  /************************************/
  // susceptible uninformed
  vertexStates.push_back(Label("S-","00;32", 0, "fillcolor=\"royalblue4\""));
  // infected uninformed
  vertexStates.push_back(Label("I-","00;31", 1, "fillcolor=\"red4\""));
  // recovered uninformed
  vertexStates.push_back(Label("R-","00;34", 2, "fillcolor=\"green4\""));
  // susceptible informed
  vertexStates.push_back(Label("S+","01;32", 3, "fillcolor=\"royalblue\""));
  // infected informed
  vertexStates.push_back(Label("I+","01;31", 4, "fillcolor=\"red\""));
  // recovered informed
  vertexStates.push_back(Label("R+","01;34", 5, "fillcolor=\"green\""));

  /*************************************/
  // define edge types
  /************************************/
  edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

  /*************************************/
  // define model parameters
  /************************************/
  model_options.add_options()
    ("beta--", po::value<double>(),
     "disease transmission rate uninformed->uninformed")
    ("beta+-", po::value<double>(),
     "disease transmission rate informed->uninformed")
    ("beta-+", po::value<double>(),
     "disease transmission rate uninformed->informed")
    ("beta++", po::value<double>(),
     "disease transmission rate informed->informed")
    ("gamma-", po::value<double>(),
     "recovery rate of uninformed")
    ("gamma+", po::value<double>(),
     "recovery rate of informed")
    ("delta-", po::value<double>(),
     "loss of immunity rate of uninformed")
    ("delta+", po::value<double>(),
     "loss of immunity rate of informed")
    ("alpha", po::value<double>(),
     "information transmission rate")
    ("nu", po::value<double>(),
     "information generation rate")
    ("omega", po::value<double>(),
     "local information generation rate")
    ("lambda", po::value<double>(),
     "loss of information rate")
    ("sigma", po::value<double>(),
     "ratio between uninformed/uninformed susceptibilities")
    ;

  /*************************************/
  // assign model parameters to variables
  /************************************/
  params.insert(std::make_pair("beta--", &beta[0][0]));
  params.insert(std::make_pair("beta+-", &beta[1][0]));
  params.insert(std::make_pair("beta-+", &beta[0][1]));
  params.insert(std::make_pair("beta++", &beta[1][1]));
  params.insert(std::make_pair("gamma-", &gamma[0]));
  params.insert(std::make_pair("gamma+", &gamma[1]));
  params.insert(std::make_pair("delta-", &delta[0]));
  params.insert(std::make_pair("delta+", &delta[1]));
  params.insert(std::make_pair("alpha", &alpha));
  params.insert(std::make_pair("nu", &nu));
  params.insert(std::make_pair("omega", &omega));
  params.insert(std::make_pair("lambda", &lambda));
  params.insert(std::make_pair("sigma", &sigma));
}

//----------------------------------------------------------
void Models::InfoSIRS::Init(po::variables_map& vm)
{
  Model::Init(vm);
  if (vm.count("sigma")) {
    // if sigma is defined, beta+- and beta++ are overwritten
    beta[1][0]=sigma*beta[0][0];
    beta[1][1]=sigma*beta[0][1];
    if (verbose >= 1) {
      std::cout << "sigma given, setting beta+-=" << beta[1][0]
                << " and beta++=" << beta[1][1] << std::endl;
    }
  }
}

//----------------------------------------------------------
double Models::InfoSIRS::getNodeEvents(eventList& events,
                                       unsigned int state,
                                       unsigned int nb) const
{
   double rateSum(.0);

   if (getDisease(state) == Recovered) {
     // loss of immunity
      Event immunityLoss;
      immunityLoss.rate = delta[getInfo(state)];
      immunityLoss.newState = getState(Susceptible, getInfo(state));
      immunityLoss.nb = nb;
      if (immunityLoss.rate > 0) {
        events.push_back(immunityLoss);
        rateSum += immunityLoss.rate;
        if (verbose >= 2) {
          std::cout << "Adding loss of immunity event with rate " 
	            << immunityLoss.rate << std::endl;
        }
      }
   } else
   if (getDisease(state) == Infected) {
      // recovery
      Event recovery;
      recovery.rate = gamma[getInfo(state)];
      recovery.newState = getState(Recovered, getInfo(state));
      recovery.nb = nb;
      if (recovery.rate > 0) {
        events.push_back(recovery);
        rateSum += recovery.rate;
        if (verbose >= 2) {
          std::cout << "Adding recovery event with rate " << recovery.rate
                    << std::endl;
        }
      }
      if (getInfo(state) == Uninformed) {
        // local information generation
        Event localInfo;
        localInfo.rate = omega;
        localInfo.newState = getState(getDisease(state), Informed);
        localInfo.nb = nb;
        if (localInfo.rate > 0) {
          events.push_back(localInfo);
          rateSum += localInfo.rate;
          if (verbose >= 2) {
            std::cout << "Adding local information event with rate " 
	              << localInfo.rate << std::endl;
          }
        }
      }
   }
   // information loss
   if (getInfo(state) == Informed) {
      Event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState = getState(getDisease(state), Uninformed);
      infoLoss.nb = nb;
      if (infoLoss.rate > 0) {
        events.push_back(infoLoss);
        rateSum += infoLoss.rate;
        if (verbose >= 2) {
          std::cout << "Adding information loss event with rate " 
	            << infoLoss.rate << std::endl;
        }
      }
   }

   return rateSum;
}

//----------------------------------------------------------
double Models::InfoSIRS::getEdgeEvents(eventList& events,
                                       unsigned int state,
                                       double detail,
                                       unsigned int edge,
                                       unsigned int nbState,
                                       double nbDetail,
                                       unsigned int nb) const
{
   double rateSum(.0);
   if (edge == Disease) {
      // infection
      if (getDisease(state) == Susceptible &&
          getDisease(nbState) == Infected) {
         Event infection;
         infection.rate = beta[getInfo(state)][getInfo(nbState)];
         infection.newState = getState(Infected, getInfo(state));
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
   } else if (edge == Information) {
      // information transmission
      if (getInfo(state) == Uninformed && getInfo(nbState) == Informed) {
         Event infoTransmission;
         infoTransmission.rate = alpha;
         infoTransmission.newState = getState(getDisease(state), Informed);
         infoTransmission.nb = nb;
         infoTransmission.et = edge;
         if (infoTransmission.rate > 0) {
           events.push_back(infoTransmission);
           rateSum += infoTransmission.rate;
           if (verbose >= 2) {
             std::cout << "Adding information transmission event with rate " 
                       << infoTransmission.rate << std::endl;
           }
         }
      }
      // information generation
      if (getInfo(state) == Uninformed && getDisease(nbState) == Infected) {
         Event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState = getState(getDisease(state), Informed);
         infoGeneration.et = edge;
         if (infoGeneration.rate > 0) {
           events.push_back(infoGeneration);
           rateSum += infoGeneration.rate;
           if (verbose >= 2) {
             std::cout << "Adding information generation event with rate " 
                       << infoGeneration.rate << std::endl;
           }
         }
      }
   }

   return rateSum;
}
