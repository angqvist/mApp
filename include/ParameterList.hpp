#include <vector>
#include <string>
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
  int nbrOfParams;
  std::string fileName;
  std::vector<Pair> paramList;
  void readParamsWithPL(class PairList);
  void readParamsATATStyle(std::vector<std::string >);
};
