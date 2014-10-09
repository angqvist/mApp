#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

class LatticeList
{
public:
  LatticeList(); // default consructor
  LatticeList(int , int , int );
  LatticeList(int , int , int ,int , std::string );
  LatticeList(int,int,int,int,int,std::string,std::vector<std::string>);
  double getDistance(int,int);
  int getNbrOfSites();
  void printList(); //Mostly for bug checking
  void printList(std::string); //Mostly for bug checking
  void setEnergy(double);
  void setBandGap(double);
  void setVolume(double);
  void setAverageLatticeConstant(double);
  double getEnergy();
  double getBandGap();
  double getVolume();
  double getAverageLatticeConstant();
  std::string getSite(int);
  void setSite(int,std::string);
  void setRandomSites(int, std::string, std::string);
  void setRandomSites(int,int, std::string, std::string,std::string);
  void setRandomSites(int,int,int, std::string, std::string,std::string,std::string);
  double getLx();
  double getLy();
  double getLz();
  std::string getAtomInfo(int, double &, double &, double &);
  double getProperty(int);
  int getNumberOfProperties();
  std::vector<double> properties;
  double getConcentration(std::string);
private:
  void readIdealPos(); // read positions from file  
  void readIdealPos2(); // read positions from file  
  void readIdealPos3(std::vector<std::string>); // read positions from file  

  std::vector<std::string> elements;
  std::vector<int> elementCounts;
  
  int cellSizeX;
  int cellSizeY;
  int cellSizeZ;
  int nbrOfSites;
  int nbrOfAtoms; //atoms in conf of newFileName
  std::string fileName;
  double latticeConstant; //MODIFY_ME 
  std::vector<double> posList2;
  std::vector<std::string> atomTypeList;
  //box length used for getting boundary condition distances
  double Lx;
  double Ly;
  double Lz;
  double energy;
  double bandGap;
  double volume;
  double averageLatticeConstant;
};
