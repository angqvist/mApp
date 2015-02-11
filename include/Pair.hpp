#pragma once
#include <string>


class Pair
{
public:
  Pair();
  Pair(double &, std::string &, std::string &);
  Pair(double &, std::string &, std::string &, int &);

  double getDistance();
  std::string getSite1();
  std::string getSite2();
  double getEnergy();
  int getCount();
  void printPair();
  void setDistance(double );
  void setSite1(std::string );
  void setSite2(std::string );
  void setEnergy(double );
  void setCount(int & );

  void incrementCount();
  void decrementCount();

  friend int operator==(Pair &,Pair &);
  friend int operator!=(Pair &,Pair &);
  
private:
  double getLimit();
  double distance;
  std::string site1;
  std::string site2;
  int count;
  double LIMIT;
  double energy;
};
