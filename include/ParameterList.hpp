#include <vector>
#include <string>
#include "Triplet.hpp"
#include "Atom.hpp"
#include "Quatuplet.hpp"
class ParameterList
{
public:
  int getNbrOfPairs();
  class Pair& getPair(int);
  ParameterList();
  ParameterList(std::string ,double );
  ParameterList(class PairList);
  ParameterList(std::string,double,std::vector<std::string >);
  void printList();
  void printPair(int);
  
  class Atom& getSinglet(int);
  int getNbrOfSinglets();
  
  class Triplet& getTriplet(int);
  int getNbrOfTriplets();  
  class Quatuplet& getQuatuplet(int);
  int getNbrOfQuatuplets();
  double getOffsetValue();
  std::vector<class Atom> returnSingletVector();
private:
  void unwrapSinglets(std::vector<Atom>,std::vector<std::string>);
  void unwrapPairs(std::vector<Pair>, std::vector<std::string>);
  void unwrapTriplets(std::vector<Triplet> ,std::vector<std::string> );
  void unwrapQuatuplets(std::vector<Quatuplet> ,std::vector<std::string> );

  double eCutOff;
  void readParams();
  void readParams_new();
  int nbrOfParams;
  std::string fileName;
  std::vector<class Atom> paramList_0;
  std::vector<Pair> paramList;
  std::vector<Triplet> paramList_3;
  std::vector<Quatuplet> paramList_4;

  void readParamsWithPL(class PairList);
  void readParamsATATStyle(std::vector<std::string >);
  void readParams_new(std::vector<std::string > );

  std::vector<Atom> atomList;
  double offset_value;
};
