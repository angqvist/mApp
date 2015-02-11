#pragma once
#include <string>
#include <vector>

class Quatuplet
{
public:
  Quatuplet();
  Quatuplet(double &,double &,double &,double &,double &,double &,std::string &,std::string &,std::string &,std::string &);
  Quatuplet(std::vector<double> &,std::vector<std::string> &);

  double getDistance1();
  double getDistance2();
  double getDistance3();
  double getDistance4();
  double getDistance5();
  double getDistance6();

  double getDistance(int &);
  


  std::string getSite1();
  std::string getSite2();
  std::string getSite3();
  std::string getSite4();

  std::string getSite(int &);
  double getEnergy();
  int getCount();
  
  std::vector<double> getDists();
  std::vector<std::string> getElements();
  void print();
  
  void setDistance1(double );
  void setDistance2(double );
  void setDistance3(double );
  void setDistance4(double );
  void setDistance5(double );
  void setDistance6(double );
  
  void setDistance(int &,double &);

  void setAll(std::vector<double> &, std::vector<std::string> &);

  void setSite1(std::string &);
  void setSite2(std::string &);
  void setSite3(std::string &);
  void setSite4(std::string &);
  void setSite(int &,std::string &);
  void setEnergy(double &);

  void setCount(int &);

  void incrementCount();
  void decrementCount();
  void increaseCountBy(int &);

  friend int operator ==(Quatuplet &,Quatuplet &);
  friend int operator !=(Quatuplet &,Quatuplet &);

  //void sortQuatuplet();

private:
  std::vector<double> dists;		      
  std::vector<std::string> elements;
  double distance1;
  double distance2;
  double distance3;
  double distance4;
  double distance5;
  double distance6;
  
  double LIMIT;
  std::string site1;
  std::string site2;
  std::string site3;
  std::string site4;

  int count;
  double energy;
  

};
