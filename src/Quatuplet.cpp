#include "Quatuplet.hpp"
#include <cmath>
#include <iostream>

Quatuplet::Quatuplet()
{
  energy=0;
  LIMIT = 1e-4;

  dists.resize(6);
  elements.resize(4);
  for(int i=0; i<6; i++)
    {
      dists[i]=0;
    }
  for(int i=0; i<4; i++)
    {
      elements[i]="A";
    }
  // distance1=0;
  // distance2=0;
  // distance3=0;
  // distance4=0;
  // distance5=0;
  // distance6=0;

  // site1="A";
  // site2="A";
  // site3="A";
  // site4="A";
  count =0;
}

Quatuplet::Quatuplet(double d1,double d2, double d3, double d4, double d5,double d6,std::string s1,
		     std::string s2,std::string s3,std::string s4)
{

  dists[0]=d1;
  dists[1]=d2;
  dists[2]=d3;
  dists[3]=d4;
  dists[4]=d5;
  dists[5]=d6;


  // distance1=d1;
  // distance2=d2;
  // distance3=d3;
  // distance4=d4;
  // distance5=d5;
  // distance6=d6;

  elements[0]=s1;
  elements[1]=s2;
  elements[2]=s3;
  elements[3]=s4;

  // site1=s1;
  // site2=s2;
  // site3=s3;
  // site4=s4;

  count=0;
  energy=0;
  LIMIT= 1e-4;
}


Quatuplet::Quatuplet(std::vector<double> new_dist,std::vector<std::string> new_sites)
{

  dists=new_dist;
  elements=new_sites;

  
  // distance1=dist[0];
  // distance2=dist[1];
  // distance3=dist[2];
  // distance4=dist[3];
  // distance5=dist[4];
  // distance6=dist[5];

  // site1=sites[0];
  // site2=sites[1];
  // site3=sites[2];
  // site4=sites[3];

  count=0;
  energy=0;
  LIMIT= 1e-4;
}


double Quatuplet::getDistance1()
{
  return dists[0];
  //return distance1;
}
double Quatuplet::getDistance2()
{
  return dists[1];
  //return distance2;
}

double Quatuplet::getDistance3()
{
  return dists[2];
  //return distance3;
}

double Quatuplet::getDistance4()
{
  return dists[3];
  //return distance4;
}

double Quatuplet::getDistance5()
{
  return dists[4];
  //return distance5;
}

double Quatuplet::getDistance6()
{
  return dists[5];
  //return distance6;
}

double Quatuplet::getDistance(int index)
{

  return dists[index];

  
  // if(index==0)
  //   {
  //     return distance1;
  //   }
  // else if(index==1)
  //   {
  //     return distance2;
  //   }
  // else if(index==2)
  //   {
  //     return distance3;
  //   }
  // else if(index==3)
  //   {
  //     return distance4;
  //   }
  // else if(index==4)
  //   {
  //     return distance5;
  //   }
  // else if(index==5)
  //   {
  //     return distance6;
  //   }
  // std::cout<<"index error in quatuplet getDistance, returnong 0 by default"<<std::endl;
  // return 0;

}

std::string Quatuplet::getSite1()
{
  return elements[0];
  //  return site1;
}
std::string Quatuplet::getSite2()
{
  return elements[1];
  //return site2;
}
std::string Quatuplet::getSite3()
{
  return elements[2];
  //return site3;
}
std::string Quatuplet::getSite4()
{
  return elements[3];
  //return site4;
}

std::string Quatuplet::getSite(int index)
{

  return elements[index];
  
  // if(index==0)
  //   {
  //     return site1;
  //   }
  // else if(index==1)
  //   {
  //     return site2;
  //   }
  // else if(index==2)
  //   {
  //     return site3;
  //   }
  // else if(index==3)
  //   {
  //     return site4;
  //   }
  // std::cout<<"Error:  quatuplet getSite, index error, returning \"A\" by default"<<std::endl;
  // return "A";
}

double Quatuplet::getEnergy()
{
  return energy;
}

int Quatuplet::getCount()
{
  return count;
}

std::vector<double> Quatuplet::getDists()
{
  return dists;
}

std::vector<std::string> Quatuplet::getElements()
{
  return elements;
}


void Quatuplet::setDistance1(double new_dist)
{
  dists[0]=new_dist;
  // distance1=new_dist;
}
void Quatuplet::setDistance2(double new_dist)
{
    dists[1]=new_dist;
    //  distance2=new_dist;
}

void Quatuplet::setDistance3(double new_dist)
{
  dists[2]=new_dist;
  // distance3=new_dist;
}

void Quatuplet::setDistance4(double new_dist)
{
  dists[3]=new_dist;
  // distance4=new_dist;
}

void Quatuplet::setDistance5(double new_dist)
{
  dists[4]=new_dist;
  //distance5=new_dist;
}

void Quatuplet::setDistance6(double new_dist)
{
  dists[5]=new_dist;
  // distance6=new_dist;
}


void Quatuplet::setDistance(int index,double val)
{

  dists[index]=val;

  // if(index==0)
  //   {
  //     distance1=val;
  //   }
  // else if(index==1)
  //   {
  //     distance2=val;
  //   }
  // else if(index==2)
  //   {
  //     distance3=val;
  //   }
  // else if(index==3)
  //   {
  //     distance4=val;
  //   }
  // else if(index==4)
  //   {
  //     distance5=val;
  //   }
  // else if(index==5)
  //   {
  //     distance6=val;
  //   }
  // std::cout<<"index error in quatuplet setDistancet"<<std::endl;
  

}




void Quatuplet::setSite1(std::string new_site)
{
  elements[0]=new_site;
  // site1=new_site;
}
void Quatuplet::setSite2(std::string new_site)
{
  elements[1]=new_site;
  // site2=new_site;
}
void Quatuplet::setSite3(std::string new_site)
{
  elements[2]=new_site;
  //site3=new_site;
}
void Quatuplet::setSite4(std::string new_site)
{
  elements[3]=new_site;
  // site4=new_site;
}


void Quatuplet::setSite(int index,std::string siteTYPE)
{

  elements[index]=siteTYPE;

  // if(index==0)
  //   {
  //     site1=siteType;
  //   }
  // else if(index==1)
  //   {
  //     site2=siteType;
  //   }
  
  // else if(index==2)
  //   {
  //     site3=siteType;
  //   }
  // else if(index==3)
  //   {
  //     site4=siteType;
  //   }
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
void Quatuplet::increaseCountBy(int addThisToCount)
{
  count += addThisToCount;
}
void Quatuplet::decrementCount()
{
  count--;
}

void Quatuplet::setAll(std::vector<double> new_dists, std::vector<std::string> new_elements)
{
  if(dists.size() != 6 || elements.size() != 4)
    {
      std::cout<<"error: size mismatch in function setAll in Quatuplet. "<<std::endl;
      std::cout<<"If it looks like a duck, and quacks like a duck, we have at least to consider the possibility that we have a small aquatic bird of the family anatidae on our hands."<<std::endl;
      return;
    }
  dists=new_dists;
  elements=new_elements;

  // for(int i=0; i<6; i++)
  //   {
  //     setDistance(i,dists[i]);
  //   }
  // for(int i=0; i<4; i++)
  //   {
  //     setSite(i,elements[i]);
  //   }
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

void Quatuplet::print()
{
  for(int i=0; i<6; i++)
    {
      std::cout<<dists[i]<< " ";
    }
  std::cout<<energy<< " ";
  std::cout<<count<<" ";
  for(int i=0; i<4; i++)
    {
      std::cout<<elements[i]<< " ";
    }
  std::cout<<std::endl;
}
