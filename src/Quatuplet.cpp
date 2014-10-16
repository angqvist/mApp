#include "Quatuplet.hpp"
#include <cmath>
#include <iostream>

Quatuplet::Quatuplet()
{
  energy=0;
  LIMIT = 1e-4;
  distance1=0;
  distance2=0;
  distance3=0;
  distance4=0;
  distance5=0;
  distance6=0;

  site1="A";
  site2="A";
  site3="A";
  site4="A";
  count =0;
}

Quatuplet::Quatuplet(double d1,double d2, double d3, double d4, double d5,double d6,std::string s1,
		     std::string s2,std::string s3,std::string s4)
{
  distance1=d1;
  distance2=d2;
  distance3=d3;
  distance4=d4;
  distance5=d5;
  distance6=d6;

  site1=s1;
  site2=s2;
  site3=s3;
  site4=s4;

  count=0;
  energy=0;
  LIMIT= 1e-4;
}


Quatuplet::Quatuplet(std::vector<double> dist,std::vector<std::string> sites)
{
  distance1=dist[0];
  distance2=dist[1];
  distance3=dist[2];
  distance4=dist[3];
  distance5=dist[4];
  distance6=dist[5];

  site1=sites[0];
  site2=sites[1];
  site3=sites[2];
  site4=sites[3];

  count=0;
  energy=0;
  LIMIT= 1e-4;
}


double Quatuplet::getDistance1()
{
  return distance1;
}
double Quatuplet::getDistance2()
{
  return distance2;
}

double Quatuplet::getDistance3()
{
  return distance3;
}

double Quatuplet::getDistance4()
{
  return distance4;
}

double Quatuplet::getDistance5()
{
  return distance5;
}

double Quatuplet::getDistance6()
{
  return distance6;
}

double Quatuplet::getDistance(int index)
{
  if(index==0)
    {
      return distance1;
    }
  else if(index==1)
    {
      return distance2;
    }
  else if(index==2)
    {
      return distance3;
    }
  else if(index==3)
    {
      return distance4;
    }
  else if(index==4)
    {
      return distance5;
    }
  else if(index==5)
    {
      return distance6;
    }
  std::cout<<"index error in quatuplet getDistance, returnong 0 by default"<<std::endl;
  return 0;

}

std::string Quatuplet::getSite1()
{
  return site1;
}
std::string Quatuplet::getSite2()
{
  return site2;
}
std::string Quatuplet::getSite3()
{
  return site3;
}
std::string Quatuplet::getSite4()
{
  return site4;
}

std::string Quatuplet::getSite(int index)
{
  if(index==0)
    {
      return site1;
    }
  else if(index==1)
    {
      return site2;
    }
  else if(index==2)
    {
      return site3;
    }
  else if(index==3)
    {
      return site4;
    }
  std::cout<<"Error:  quatuplet getSite, index error, returning \"A\" by default"<<std::endl;
  return "A";
}

double Quatuplet::getEnergy()
{
  return energy;
}

int Quatuplet::getCount()
{
  return count;
}




void Quatuplet::setDistance1(double new_dist)
{

  distance1=new_dist;
}
void Quatuplet::setDistance2(double new_dist)
{

  distance2=new_dist;
}

void Quatuplet::setDistance3(double new_dist)
{

  distance3=new_dist;
}

void Quatuplet::setDistance4(double new_dist)
{

  distance4=new_dist;
}

void Quatuplet::setDistance5(double new_dist)
{

  distance5=new_dist;
}

void Quatuplet::setDistance6(double new_dist)
{

  distance6=new_dist;
}


void Quatuplet::setDistance(int index,double val)
{
  if(index==0)
    {
      distance1=val;
    }
  else if(index==1)
    {
      distance2=val;
    }
  else if(index==2)
    {
      distance3=val;
    }
  else if(index==3)
    {
      distance4=val;
    }
  else if(index==4)
    {
      distance5=val;
    }
  else if(index==5)
    {
      distance6=val;
    }
  std::cout<<"index error in quatuplet setDistancet"<<std::endl;
  

}




void Quatuplet::setSite1(std::string new_site)
{
  site1=new_site;
}
void Quatuplet::setSite2(std::string new_site)
{
  site2=new_site;
}
void Quatuplet::setSite3(std::string new_site)
{
  site3=new_site;
}
void Quatuplet::setSite4(std::string new_site)
{
site4=new_site;
}

void Quatuplet::setEnergy(double new_energy)
{
energy=new_energy;
}

void Quatuplet::setCount(int newCount)
{
  count=newCount;
}

void Quatuplet::incrementCount()
{
  count++;
}
void Quatuplet::decrementCount()
{
  count--;
}


int operator==(Quatuplet q1 ,Quatuplet q2)
{
  for(int i=0; i<6; i++)
    {
      if(fabs(q1.getDistance(i)-q2.getDistance(i))>1e-4)
	{
	  return false;
	}
    }
  for(int i=0; i<4; i++)
    {
      if(q1.getSite(i) !=q2.getSite(i))
	{
	  return false;
	}
    }
  return true;
}

int operator !=(Quatuplet q1, Quatuplet q2)
{
  return !(q1==q2);
}
