/******************************************************************/
// Model.hh
// 
// this contains all classes related to the model:
//   -- VertexState
//        for the state a vertex can assume
//   -- event
//        for the events that can happen in the simulation
//   -- Model
//        for the model parameters and the events associated with it
/******************************************************************/
#ifndef MODEL_HH
#define MODEL_HH

#include <list>
#include <vector>
#include <string>

enum diseaseStates {Susceptible,Infected,Recovered};
enum infoStates {Uninformed, Informed};
enum edgeTypes {Disease, Information};

class VertexState
{
   public:

      VertexState() { disease = Susceptible; info = Uninformed; }
      VertexState(int dState, int iState) { disease = dState; info = iState; }
      ~VertexState() {;}
      
      int getDisease() const { return disease; }
      int getInfo() const { return info; }

      void setDisease(int newDisease) { disease = newDisease; }
      void setInfo(int newInfo) { info = newInfo; }

      bool operator==(const VertexState& rhs) const
      { return disease==rhs.getDisease() && info == rhs.getInfo(); }

      VertexState& operator=(const VertexState& rhs) 
      { setDisease(rhs.getDisease()); setInfo(rhs.getInfo()); return *this; }

      bool operator<(const VertexState& rhs) const
      { return ((getInfo() < rhs.getInfo()) ||
            (getInfo() == rhs.getInfo() && getDisease() < rhs.getDisease())); }
      
      void set(std::string s);
      std::string getString() const;
      void print(std::ostream& os) const;

   private:
      int disease; // Suspected, Infected or Recovered
      int info; // Informed or Uninformed
};

// streams a letter associated with the state, needed for easy printout
std::ostream& operator<<(std::ostream& os, const VertexState& v);
      
struct event
{
      double rate; // rate at which an event occurs (depends on
                   // model parameters) 
      VertexState newState; // state an event will change a vertex to
};

typedef std::list<event> eventList;

class Model
{
   public:
      Model();
      ~Model();
      
      void InitDefaultParams();
      int InitFromFile(std::string fileName);
         
      double getNodeEvents(eventList& events,
                           VertexState state) const;
      double getEdgeEvents(eventList& events,
                           VertexState state, int edge,
                           VertexState nbState) const;

      std::vector<VertexState> getPossibleStates();

      double gamma[2], delta[2]; // model parameters
      double beta[2][2], alpha, nu, lambda; // model parameters
};

#endif
