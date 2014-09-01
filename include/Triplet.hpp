#include <string>


class Triplet
{
public:
  Triplet();
  Triplet(double ,double,double, std::string , std::string,std::string);
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
  void setDistance1(double distance);
  void setDistance2(double distance);
  void setDistance3(double distance);

  void setSite1(std::string newSite);
  void setSite2(std::string newSite);
  void setSite3(std::string newSite);
  void setEnergy(double newEnergy);
  void setCount(int newCount);

  void incrementCount();
  void decrementCount();

  friend int operator==(Triplet,Triplet);
  friend int operator!=(Triplet,Triplet);
  double LIMIT;

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
};
