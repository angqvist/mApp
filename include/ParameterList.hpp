#include <vector>
#include <string>
#include "Triplet.hpp"
#include "Atom.hpp"
//#include "Quatuplet.hpp"
class ParameterList
{
public:
  int getNbrOfParams();
  class Pair& getPair(int);
  ParameterList();
  ParameterList(std::string ,double );
  ParameterList(class PairList);
  ParameterList(std::string,double,std::vector<std::string >);
  void printList();
  void printPair(int);

private:
  double eCutOff;
  void readParams();
  void readParams_new();
  int nbrOfParams;
  std::string fileName;
  std::vector<class Atom> paramList_0;
  std::vector<Pair> paramList;
  std::vector<Triplet> paramList_3;
  void readParamsWithPL(class PairList);
  void readParamsATATStyle(std::vector<std::string >);
  void readParams_new(std::vector<std::string > );

  std::vector<Atom> atomList;
  double offset_value;
};
