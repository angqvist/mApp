#include <iostream>
#include "Pair.hpp"
#include <cmath>

Pair::Pair()
{
  energy=0;
  LIMIT=1e-5;
  distance=0;
  site1="Ga";
  site2="Ga";
  count=0;
}

Pair::Pair(double newDistance, std::string newSite1, std::string newSite2, int newCount)
{
  energy=0;
  LIMIT=1e-4;
  distance=newDistance;
  site1=newSite1;
  site2=newSite2;
  count=newCount;
}

Pair::Pair(double newDistance, std::string newSite1, std::string newSite2)
{
  energy=0;
  LIMIT=1e-4;
  distance=newDistance;
  site1=newSite1;
  site2=newSite2;
  count=0;
}

//Getters
double Pair::getDistance()
{
  return distance;
}

double Pair::getEnergy()
{
  return energy;
}


std::string Pair::getSite1()
{
  return site1;
}
std::string Pair::getSite2()
{
  return site2;
}
int Pair::getCount()
{
  return count;
}

// setters
void Pair::setDistance(double newDistance)
{
   distance=newDistance;
}
void Pair::setSite1(std::string newSite)
{
  site1=newSite;
}
void Pair::setSite2(std::string newSite)
{
   site2=newSite;
}
void Pair::setCount(int newCount)
{
   count=newCount;
}

void Pair::setEnergy(double newEnergy)
{
  energy=newEnergy;
}

void Pair::incrementCount()
{
  count++;
}

void Pair::decrementCount()
{
  count--;
}

//this one feels stupid
double Pair::getLimit()
{
  return LIMIT;
}

int operator==(Pair p1, Pair p2)
{
  if(fabs(p1.getDistance() - p2.getDistance())>p1.getLimit())
    {
      return false;
    }
  if((p1.getSite1() != p2.getSite1() || p1.getSite2() != p2.getSite2()) && ((p1.getSite2() != p2.getSite1()) || p1.getSite1() != p2.getSite2()))
    {
      return false;
    }
  return true;
}

//Keep it simple
int operator!=(Pair p1, Pair p2)
{
  return !(p1==p2);
}

void Pair::printPair()
{
  std::cout<<getDistance()<<" "<<getSite1()<<" "<<getSite2()<< " "<<getEnergy()<< " "<<getCount()<<std::endl;
}
