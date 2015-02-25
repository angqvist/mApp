#pragma once
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

class LatticeList
{
public:
  LatticeList(); // default consructor
  LatticeList(std::vector<int>,std::vector<double>,std::vector<std::string>,int,int);
  LatticeList(std::vector<std::string>,std::vector<double>,std::vector<std::string>,std::vector<double>);

  LatticeList(int,int,int);
  LatticeList(int,int,int,int,std::string);
  LatticeList(int,int,int,int,int,std::string,std::vector<std::string>);
  LatticeList(int,int,int,int,int,std::string,std::vector<std::string>,std::vector<std::vector<double> >);
  LatticeList(std::vector<std::string> , std::vector<double> , std::vector<double>);

  //LatticeList(std::vector<std::string> , std::vector<double> , std::vector<std::string>,std::vector<double> );

  std::vector<double> getPeriodicDistance(int &, int &,  double &);
  double getDistance(const int &,const int &);
  
  int getNbrOfSites() ;

  void printList(); //Mostly for bug checking
  void printList(std::string); //Mostly for bug checking
  void setEnergy(double &);
  void setBandGap(double);
  void setVolume(double);
  void setAverageLatticeConstant(double);
  void setWykPos(std::vector<std::string>);
  std::vector<double> getWyckoffOccupancy(const std::vector<std::string> &, const std::string &);
  double getEnergy();
  double getBandGap();
  double getVolume();
  double getAverageLatticeConstant();
  std::string getSite(int &);
  void set_property(std::vector<double> &);
  void setSite(int &,std::string &);
  void setRandomSites(int &, std::string &, std::string &);
  void setRandomSites(int &,int &, std::string &, std::string &,std::string &);
  void setRandomSites(int &,int &,int &, std::string &, std::string &,std::string &,std::string &);
  double getLx();
  double getLy();
  double getLz();
  std::string getAtomInfo(int, double &, double &, double &);
  double getProperty(int);
  int getNumberOfProperties();
  double getConcentration(std::string);
  void calculate_lookup_table();
  double fast_distance( int &,  int &) const;
  std::vector<std::vector<double> > getLookupTable();
  void append_atom(std::string,std::vector<double>, int);
  void clear_lookup_table();
  int get_original_atoms_count();
  int get_ghost_atoms_count();
  std::vector<double> getCellMatrix();
private:
  std::vector<double> cellMatrix;
  int number_of_original_atoms;
  int number_of_ghost_atoms;
  bool distance_table_init;
  std::vector<std::vector<double> > distance_table;
  std::vector<double> properties;
  std::vector<std::string> wyckoffSite;
  void readIdealPos(); // read positions from file  
  void readIdealPos2(); // read positions from file  
  void readIdealPos3(std::vector<std::string>); // read positions from file  

  void create_tag_list();
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
  std::vector<int> tags;
  std::vector<std::vector<int>> tag_list;
};
