#pragma once
#include <string>
#include <vector>

class Atom
{
public:
  Atom();
  Atom(std::string,double);
  //Atom(std::string,std::vector<double>);
  void setType(std::string);

  void setProperty(double);
  //void setProperty(int,double);
  
  std::string getType();
  double getProperty();


private:
  std::vector<double> properties;
  std::string type;
  double property;
};
