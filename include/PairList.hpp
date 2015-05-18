#pragma once

#include <vector>
#include <string>
//#include "Pair.hpp"
//#include "LatticeList.hpp"
//hehe
class PairList
{

public:
  PairList();
  int getNbrOfPairs();
  class Pair& getPair(int &);
  int updatePair(class Pair &,bool &);
  // string with possible elements
  void initializePairs(class LatticeList &, std::vector<std::string> &,double &);
  void printList();
  void resetCounts();
  void divideCountByTwo();
  std::vector<double> getClusterVector(std::vector<std::string> &, double &,bool &);
  std::vector<double> getUniqueDistances(double &);
  
private:
  int nbrOfPairs;
  std::vector<Pair> pairList;
  int isAtomInSubElements(std::string ,std::vector<std::string> &);
  void sortPairs();
};
