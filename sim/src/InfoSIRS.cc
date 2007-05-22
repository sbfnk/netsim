/*******************************************************************/
//
// Class InfoSIRS
// --------------------
//
// This is the class implementing the SIRS model, i.e. the parameters
// and the dependance of the rates of various processes on node states
// and edge properties
//
/******************************************************************/

#include <iostream>
#include <fstream>

#include "InfoSIRS.hh"

/******************************************************************/
// InfoSIRS constructor
/******************************************************************/
InfoSIRS::InfoSIRS(unsigned int v)
  : Model(v)
{
  // susceptible uninformed
  vertexStates.push_back(Label("S-","00;32", 0, "fillcolor=\"green4\""));
  // infected uninformed
  vertexStates.push_back(Label("I-","00;31", 1, "fillcolor=\"red4\""));
  // recovered uninformed
  vertexStates.push_back(Label("R-","00;34", 2, "fillcolor=\"royalblue4\""));
  // susceptible informed
  vertexStates.push_back(Label("S+","01;32", 3, "fillcolor=\"green\""));
  // infected informed
  vertexStates.push_back(Label("I+","01;31", 4, "fillcolor=\"red\""));
  // recovered informed
  vertexStates.push_back(Label("R+","01;34", 5, "fillcolor=\"royalblue\""));

  edgeTypes.push_back(Label("d", "", 0, "style=\"solid\""));
  edgeTypes.push_back(Label("i", "", 1, "style=\"dashed\""));

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

/******************************************************************/
// Vertex destructor
/******************************************************************/
InfoSIRS::~InfoSIRS()
{}

/******************************************************************/
// InfoSIRS::Init
// initializes the model, uses sigma appropriately
/******************************************************************/
void InfoSIRS::Init(po::variables_map& vm)
{
  Model::Init(vm);
  if (vm.count("sigma")) {
    beta[1][0]=sigma*beta[0][0];
    beta[1][1]=sigma*beta[0][1];
    if (verbose >= 1) {
      std::cout << "sigma given, setting beta+-=" << beta[1][0]
                << " and beta++=" << beta[1][1] << std::endl;
    }
  }
}

/******************************************************************/
// InfoSIRS::getNodeEvents
// get the events that can happen for a given state of a node. Stores the
// events in the events list and returns the sum of their rates
/******************************************************************/
double InfoSIRS::getNodeEvents(eventList& events,unsigned int state) const
{
   double rateSum(.0);

   // loss of immunity
   if (getDisease(state) == Recovered) {
      event immunityLoss;
      immunityLoss.rate = delta[getInfo(state)];
      immunityLoss.newState = getState(Susceptible, getInfo(state));
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
      event recovery;
      recovery.rate = gamma[getInfo(state)];
      recovery.newState = getState(Recovered, getInfo(state));
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
        event localInfo;
        localInfo.rate = omega;
        localInfo.newState = getState(getDisease(state), Informed);
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
      event infoLoss;
      infoLoss.rate = lambda;
      infoLoss.newState = getState(getDisease(state), Uninformed);
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

/******************************************************************/
// InfoSIRS::getNodeEvents
// get the events that can happen for a given edge type and the states of the
// neighbouring nodes. Stores the events in the events list and returns the sum
// of their rates 
/******************************************************************/
double InfoSIRS::getEdgeEvents(eventList& events,
                            unsigned int state, unsigned int edge,
                            unsigned int nbState) const
{
   double rateSum(.0);
   if (edge == Disease) {
      // infection
      if (getDisease(state) == Susceptible &&
          getDisease(nbState) == Infected) {
         event infection;
         infection.rate = beta[getInfo(state)][getInfo(nbState)];
         infection.newState = getState(Infected, getInfo(state));
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
         event infoTransmission;
         infoTransmission.rate = alpha;
         infoTransmission.newState = getState(getDisease(state), Informed);
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
         event infoGeneration;
         infoGeneration.rate = nu;
         infoGeneration.newState = getState(getDisease(state), Informed);
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
