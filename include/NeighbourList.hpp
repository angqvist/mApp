#include <vector>
#include <string>
class Neighbour; //fishy - should work without but doesnt

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
  double getMatchingEnergy(class Neighbour, std::vector<Neighbour>);
  double localEnergy;
  double distanceLimit;
  int thisIndex;
  int newNeighbour(class Neighbour);
  void findNeighbours(class LatticeList,class ParameterList);
  std::vector<Neighbour> nbrList;
  std::vector<std::vector<Neighbour> > manyNbrLists;
  std::vector<int> indexList;  
  void printManyNbr(std::vector<Neighbour>);
  std::vector<double> singletEnergy;
  std::vector< std::string > singletType;
  double currentLocalEnergy;
};
   
