#include <vector>
#include <string>

class TripletList
{
public:
  TripletList();
  int getNbrOfTriplets();
  class Triplet& getTriplet(int);
  int updateTriplet(class Triplet,bool);
  // string with possible elements
  void initializeTriplets(class LatticeList, std::vector<std::string>,double);
  void printList();
  void resetCounts();
  std::vector<std::vector<double> > getUniqueDistances();
  
private:
  int nbrOfTriplets;
  std::vector<Triplet> tripletList;
  int isAtomInSubElements(std::string,std::vector<std::string>);
  void sortOrder(std::vector<double> &,std::vector<int> &);
  bool isTripletUnique(class Triplet, std::vector<std::vector<double > > &,bool);
  
};
