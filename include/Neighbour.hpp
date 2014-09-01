#include<string>

class Neighbour
{
public:
  Neighbour();
  Neighbour(std::string,std::string,double,double);
  double getEnergy();
  std::string getSite1();
  std::string getSite2();
  double getDistance();
  void setSite1(std::string);
  void setSite2(std::string);
  void setDistance(double);
  void print();
private:

  std::string site1;
  std::string site2;
  double energy;
  double distance;
  friend int operator==(Neighbour,Neighbour);
  friend int operator!=(Neighbour,Neighbour);
};
