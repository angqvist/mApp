#include "Triplet.hpp"
#include <cmath>
#include <iostream>


Triplet::Triplet()
{
  energy=0;
  LIMIT=1e-4;
  distance1=0;
  distance2=0;
  distance3=0;

  site1="A";
  site2="A";
  site3="A";
  count=0;
}

Triplet::Triplet(double dr1 ,double dr2,double dr3, std::string s1, std::string s2,std::string s3)
{
  energy=0;
  LIMIT=1e-4;
  distance1=dr1;
  distance2=dr2;
  distance3=dr3;
  site1=s1;
  site2=s2;
  site3=s3;
  count=0;
  sortTriplet();
}





//Getters
double Triplet::getDistance1()
{
  return distance1;
}
double Triplet::getDistance2()
{
  return distance2;
}

double Triplet::getDistance3()
{
  return distance3;
}


double Triplet::getEnergy()
{
  return energy;
}


std::string Triplet::getSite1()
{
  return site1;
}

std::string Triplet::getSite2()
{
  return site2;
}
std::string Triplet::getSite3()
{
  return site3;
}
int Triplet::getCount()
{
  return count;
}

void Triplet::setDistance1(double newDistance)
{
   distance1=newDistance;
}
void Triplet::setDistance2(double newDistance)
{
   distance2=newDistance;
}
void Triplet::setDistance3(double newDistance)
{
   distance3=newDistance;
}


void Triplet::setSite1(std::string newSite)
{
  site1=newSite;
}
void Triplet::setSite2(std::string newSite)
{
   site2=newSite;
}
void Triplet::setSite3(std::string newSite)
{
   site3=newSite;
}
void Triplet::setCount(int newCount)
{
   count=newCount;
}

void Triplet::setEnergy(double newEnergy)
{
  energy=newEnergy;
}

void Triplet::incrementCount()
{
  count++;
}

void Triplet::decrementCount()
{
  count--;
}
double Triplet::getLimit()
{
  return LIMIT;
}


int operator==(Triplet t1 ,Triplet t2)
{
  if(fabs(t1.getDistance1()-t2.getDistance1())> t1.getLimit())
    {
      return false;
    }
  if(fabs(t1.getDistance2()-t2.getDistance2())>t1.getLimit())
    {
      return false;
    }
  if(fabs(t1.getDistance3()-t2.getDistance3())>t1.getLimit())
    {
      return false;
    }

  //this is for and alphabetically ordered triplet

  if(t1.getSite1() !=t2.getSite1())
    {
      return false;
    }
  if(t1.getSite2() != t2.getSite2())
    {
      return false;
    }
  if(t1.getSite3() != t2.getSite3())
    {
      return false;
    }



      return true;
}

int operator!=(Triplet t1, Triplet t2)
{
  return !(t1==t2);
}

void Triplet::printTriplet()
{
  std::cout<<getDistance1()<< " "<<getDistance2()<< " "<<getDistance3()<<" "<<getSite1()<<" "<< getSite2()<< " "<<getSite3()<< " "<<getEnergy()<< " "<<getCount()<<std::endl;
}


void Triplet::sortTriplet()
{
  double diff1;
  double diff2;
  double diff3;
  
  if(this->distance1 > this->distance2 || this->distance1>this->distance3 || this->distance2 > this->distance3)
    {
      std::cout<<"Error distances not in order in sort triplet"<<std::endl;
      std::cout<<diff1<< " "<<diff2<<std::endl;
    }


  diff1= fabs(distance1-distance2);
  diff2= fabs(distance2-distance3);
  std::vector< std::string > types;
  if(diff1 < this->LIMIT && diff2 < this->LIMIT)
    {
      //all distances are considered equal
      //alphabetically sort site1, site2, site3
      types.push_back(site1);
      types.push_back(site2);
      types.push_back(site3);
      sortStrings(types);
      site1=types[0];
      site2=types[1];
      site3=types[2];
      
    }


  //if d2==d3, site1 kan vara på site2
  else if(diff2 < this->LIMIT) 
    {
      types.push_back(site1);
      types.push_back(site2);
      sortStrings(types);
      site1=types[0];
      site2=types[1];
      //alphabetically sort site1, site 2
    }
  //if d1==d2, site2 kan vara på site3
  else if(diff1 < this->LIMIT)
    {

      //alphatetically sort site2, site3

      types.push_back(site2);
      types.push_back(site3);
      sortStrings(types);
      site2=types[0];
      site3=types[1];

    }
}





void Triplet::sortStrings(std::vector< std::string> &type)
{
  bool swapped = true;
  std::string tempString;
  
  while(swapped)
    {
      swapped= false;
      for(int i=0; i<type.size()-1; i++)
	{
	  if(type[i].compare(type[i+1])>0)
	    {
	      tempString = type[i];
	      type[i]=type[i+1];
	      type[i+1]=tempString;
	      swapped=true;
	    }
	}
    }
}


void Triplet::setAll(std::vector<double> dists, std::vector<std::string> elements)
{


  if(dists.size() != 3 || elements.size() !=3)
    {
      std::cout<<"Error: size mismatch in setAll in triplet class"<<std::endl;
      return;
    }
  setDistance1(dists[0]);
  setDistance2(dists[1]);
  setDistance3(dists[2]);

  setSite1(elements[0]);
  setSite2(elements[1]);
  setSite3(elements[2]);
}
