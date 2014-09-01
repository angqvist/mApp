#include <string>


class Pair
{
public:
  Pair();
  Pair(double distance, std::string site1, std::string site2);
  Pair(double distance, std::string site1, std::string site2, int count);

  double getDistance();
  std::string getSite1();
  std::string getSite2();
  double getEnergy();
  int getCount();
  void printPair();
  void setDistance(double distance);
  void setSite1(std::string newSite);
  void setSite2(std::string newSite);
  void setEnergy(double newEnergy);
  void setCount(int newCount);

  void incrementCount();
  void decrementCount();

  friend int operator==(Pair,Pair);
  friend int operator!=(Pair,Pair);
  
private:
  double getLimit();
  double distance;
  std::string site1;
  std::string site2;
  int count;
  double LIMIT;
  double energy;
};
