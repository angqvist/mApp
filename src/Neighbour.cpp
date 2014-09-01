#include <string>
#include <iostream>
#include "Neighbour.hpp"

Neighbour::Neighbour()
{
  site1="s1";
  site2="s2";
  energy=0.0;
  distance=0.0;
}

Neighbour::Neighbour(std::string newS1, std::string newS2,double newEnergy,double newDistance)
{   
  site1=newS1;
  site2=newS2;
  energy=newEnergy;
  distance=newDistance;
}

double Neighbour::getEnergy()
{
  return energy;
}

std::string Neighbour::getSite1()
{
  return site1;
}
std::string Neighbour::getSite2()
{
  return site2;
}
double Neighbour::getDistance()
{
  return distance;
}

int operator==(Neighbour nbr1, Neighbour nbr2)
{     
  if((nbr1.getSite1() != nbr2.getSite1() || nbr1.getSite2() != nbr2.getSite2()) && ((nbr1.getSite2() != nbr2.getSite1()) || nbr1.getSite1() != nbr2.getSite2()))
    {
      return false;
    }
  return true;
}

int operator !=(Neighbour nbr1, Neighbour nbr2)
{
  return !(nbr1==nbr2);
}

void Neighbour::print()
{
  std::cout<<getSite1()<< " "<<getSite2()<< " "<<getEnergy()<<" "<<getDistance()<<std::endl;
}

void Neighbour::setSite1(std::string newSite1)
{
  site1=newSite1;
}
void Neighbour::setSite2(std::string newSite2)
{
  site2=newSite2;
}

void Neighbour::setDistance(double newDistance)
{
  distance = newDistance;
}
