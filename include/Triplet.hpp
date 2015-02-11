#pragma once
#include <string>
#include <vector>


class Triplet
{
public:
  Triplet();
  Triplet(double &,double &,double &, std::string &, std::string &,std::string &);
  //  Triplet(double distance, std::string site1, std::string site2, int count);

  double getDistance1();
  double getDistance2();
  double getDistance3();

  std::string getSite1();
  std::string getSite2();
  std::string getSite3();

  double getEnergy();
  int getCount();
  void printTriplet();
  void setDistance1(double );
  void setDistance2(double );
  void setDistance3(double );

  void setSite1(std::string );
  void setSite2(std::string );
  void setSite3(std::string );
  void setEnergy(double );
  void setCount(int );

  void incrementCount();
  void decrementCount();
  void increaseCountBy(int &);

  friend int operator==(Triplet &,Triplet &);
  friend int operator!=(Triplet &,Triplet &);
  double LIMIT;
  void sortTriplet();
  void setAll(std::vector<double> &,std::vector<std::string> &);

private:
  double getLimit();
  double distance1;
  double distance2;
  double distance3;

  std::string site1;
  std::string site2;
  std::string site3;
  int count;
  double energy;
  
  
  void sortStrings(std::vector< std::string >&);
  
};
