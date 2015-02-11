#include <vector>
#include <string>

class TripletList
{
public:
  TripletList();
  TripletList(class LatticeList &, std::vector<std::string> &,double &);
  int getNbrOfTriplets();
  class Triplet& getTriplet(int &);
  int updateTriplet(class Triplet &,bool &,int &);
  // string with possible elements
  void initializeTriplets(class LatticeList &, std::vector<std::string> &,double &);
  void printList();
  void resetCounts();
  std::vector<std::vector<double> > getUniqueDistances(double &);
  void sortTripletList();
  std::vector<double> getClusterVector(std::vector<std::string > &,double &,bool &);
  
  
private:
  int nbrOfTriplets;
  void reset();
  std::vector<Triplet> tripletList;
  int isAtomInSubElements(std::string ,std::vector<std::string> &);
  void sortOrder(std::vector<double> &,std::vector<int> &);
  bool isTripletUnique(class Triplet &, std::vector<std::vector<double > > &,bool &);
};
