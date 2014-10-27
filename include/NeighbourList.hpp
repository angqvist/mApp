#include <vector>
#include <string>
#include "Triplet.hpp"
#include "Quatuplet.hpp"
//class Neighbour; //fishy - should work without but doesnt

class NeighbourList
{
public:
  NeighbourList();
  NeighbourList(int);
  NeighbourList(int,class LatticeList,class ParameterList);
  void addNeighbour(class Neighbour);
  void printList();
  double getLocalEnergy(class LatticeList);
  int isThisMatchingNeighbour(class LatticeList, class Neighbour);
  double getCurrentLocalEnergy();
  void setCurrentLocalEnergy(double);
  void calcCurrentLocalEnergy(class LatticeList);
  
  
private:
  std::vector<class Atom> singletList;
  void readSingletList(class ParameterList);
  
  double getMatchingEnergy(class Neighbour, std::vector< class Neighbour>);
  double localEnergy;
  double distanceLimit;
  int thisIndex;
  int newNeighbour(class Neighbour);
  void findNeighbours(class LatticeList,class ParameterList);
  std::vector<class Neighbour> nbrList;
  std::vector<std::vector< class Neighbour> > manyNbrLists;
  std::vector<int> indexList;  
  void printManyNbr(std::vector<class Neighbour>);
  std::vector<double> singletEnergy;
  std::vector< std::string > singletType;
  double currentLocalEnergy;
  
  std::vector<int> pair_index;
  std::vector<double> pair_dist;
  std::vector<std::vector<std::string> > pair_elements;
  std::vector<std::vector<class Pair> > pair_vector;

  std::vector<std::vector<int> > trip_index;
  std::vector<std::vector<double> > trip_dists;
  std::vector<std::vector< std::string > > trip_elements;
  std::vector<std::vector<class Triplet > > trip_vector;

  std::vector<std::vector<int> > quat_index;
  
  std::vector<std::vector<class Quatuplet> > quat_vector;

  void findTripletNeighbours(class LatticeList,class ParameterList);
  void findQuatupletNeighbours(class LatticeList,class ParameterList);
  double offset_value;
  void setOffset(double );
  double getOffset();

};
   
