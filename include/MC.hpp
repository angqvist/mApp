#include <vector>
#include <string>

class MC
{
public:
  MC();
  double step(class LatticeList &,std::vector<class NeighbourList>);
  double step(int,class LatticeList &,std::vector<class NeighbourList>);
  double averageStep(int,int,class LatticeList &,std::vector<std::vector<class NeighbourList> >,std::vector<double> &,std::string);
  // double averageStep(int,int,class LatticeList &,std::vector<class NeighbourList>,double &,double &,double &,double &,double &,double &,double &,double &);


  // do steps with another ensamble/method?
  
  double getSwapRatio(); //returns nbrOfSuccesfulSwaps/totalSwaps
  double getTemperature();
  std::string getEnsemble();

  void resetSwapCounts();
  void setTemperature(double);
  void setEnsemble(std::string);
  void setKBeta(double);
  double T;
  double kbeta;//8.6173324e-5;
  
  std::string ensemble;
  int nbrOfSuccesfulSwaps;
  int totalSwaps;

private:
  double stepSGC(class LatticeList &,std::vector<class NeighbourList>);
  void printInfo(std::string,class LatticeList, double, double);

};
  
