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

  site1="Ga";
  site2="Ga";
  site3="Ga";
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
  // if(fabs(t1.getDistance1()-t1.getDistance3())<t1.getLimit())
  //   {
  //     std::cout<<"=================================================="<<std::endl;

  //     t1.printTriplet();
  //     t2.printTriplet();
  //   }

    //   std::cout<<"=================================================="<<std::endl;

   //    t1.printTriplet();
    //  t2.printTriplet();

  if( (t1.getSite1() != t2.getSite1() || t1.getSite2() != t2.getSite2() || t1.getSite3() != t2.getSite3()) && 
      (fabs(t1.getDistance1()-t2.getDistance2())>t1.getLimit() || (t1.getSite1() != t2.getSite3() || t1.getSite3() != t2.getSite1() || t1.getSite2() != t2.getSite2() ))&&
((fabs(t1.getDistance2()-t2.getDistance3())>t1.getLimit() || (t1.getSite1() != t2.getSite2() || t1.getSite2() != t2.getSite1() || t1.getSite3() != t2.getSite3() ))) &&
      (fabs(t1.getDistance1()-t2.getDistance3())>t1.getLimit() ||!(
  								  ((t1.getSite1() == t2.getSite1()) && ((t1.getSite2() == t2.getSite2() && t1.getSite3() == t2.getSite3() ) || (t1.getSite3() == t2.getSite2() && t1.getSite2() == t2.getSite3() ))) ||
  								  ((t1.getSite3() == t2.getSite1()) && ((t1.getSite1() == t2.getSite2() && t1.getSite2() == t2.getSite3() ) || (t1.getSite2() == t2.getSite2() && t1.getSite1() == t2.getSite3() ))) ||
  								  ((t1.getSite2() == t2.getSite1()) && ((t1.getSite3() == t2.getSite2() && t1.getSite1() == t2.getSite3() ) || (t1.getSite1() == t2.getSite2() && t1.getSite3() == t2.getSite3() ))))))




//  if( (t1.getSite1() != t2.getSite1() || t1.getSite2() != t2.getSite2() || t1.getSite3() != t2.getSite3()) && 
//       ((fabs(t1.getDistance1()-t2.getDistance2())>t1.getLimit() || (t1.getSite2() != t2.getSite3() || t1.getSite3() != t2.getSite2() || t1.getSite1() != t2.getSite1() )))&&
// ((fabs(t1.getDistance2()-t2.getDistance3())>t1.getLimit() || (t1.getSite1() != t2.getSite2() || t1.getSite2() != t2.getSite1() || t1.getSite3() != t2.getSite3() ))) &&
//       (fabs(t1.getDistance1()-t2.getDistance3())>t1.getLimit() ||(
//   								  ((t1.getSite1() != t2.getSite1()) || ((t1.getSite2() != t2.getSite2() || t1.getSite3() != t2.getSite3() ) || (t1.getSite3() != t2.getSite2() || t1.getSite2() != t2.getSite3() ))) ||
//   								  ((t1.getSite3() != t2.getSite1()) || ((t1.getSite1() != t2.getSite2() || t1.getSite2() != t2.getSite3() ) || (t1.getSite2() != t2.getSite2() || t1.getSite1() != t2.getSite3() ))) ||
//   								  ((t1.getSite2() != t2.getSite1()) || ((t1.getSite3() != t2.getSite2() || t1.getSite1() != t2.getSite3() ) || (t1.getSite1() != t2.getSite2() || t1.getSite3() != t2.getSite3() ))))))





  // if( (t1.getSite1() != t2.getSite1() || t1.getSite2() != t2.getSite2() || t1.getSite3() != t2.getSite3()) || 
  //     ((fabs(t1.getDistance1()-t2.getDistance2())>t1.getLimit() || (t1.getSite2() != t2.getSite3() || t1.getSite3() != t2.getSite2() || t1.getSite1() != t2.getSite1() )))||
  //     (fabs(t1.getDistance1()-t2.getDistance3())>t1.getLimit()))

    {
      // if(fabs(t1.getDistance1()-t1.getDistance3())<t1.getLimit())
      //  	{
      //  	  std::cout<<"False"<<std::endl;
      // 	}
      return false;
    }
  // if( !((fabs(t1.getDistance2()-t2.getDistance3())>t1.getLimit())))
  //   {
  // std::cout<<"==============================="<<std::endl;
  // t1.printTriplet();
  // t2.printTriplet();
  // std::cout<<(t1.getSite1() != t2.getSite1() || t1.getSite2() != t2.getSite2() || t1.getSite3() != t2.getSite3())<< " "<<((fabs(t1.getDistance1()-t2.getDistance2())>t1.getLimit() || (t1.getSite2() != t2.getSite3() || t1.getSite3() != t2.getSite2() || t1.getSite1() != t2.getSite1() )))<<" "<<((fabs(t1.getDistance2()-t2.getDistance3())>t1.getLimit() || (t1.getSite1() != t2.getSite2() || t1.getSite2() != t2.getSite1() || t1.getSite3() != t2.getSite3() ))) << " "<< (fabs(t1.getDistance1()-t2.getDistance3())>t1.getLimit() ||(
  // 								  ((t1.getSite1() != t2.getSite1()) || ((t1.getSite2() != t2.getSite2() || t1.getSite3() != t2.getSite3() ) || (t1.getSite3() != t2.getSite2() || t1.getSite2() != t2.getSite3() ))) ||
  // 								  ((t1.getSite3() != t2.getSite1()) || ((t1.getSite1() != t2.getSite2() || t1.getSite2() != t2.getSite3() ) || (t1.getSite2() != t2.getSite2() || t1.getSite1() != t2.getSite3() ))) ||
  // 								  ((t1.getSite2() != t2.getSite1()) || ((t1.getSite3() != t2.getSite2() || t1.getSite1() != t2.getSite3() ) || (t1.getSite1() != t2.getSite2() || t1.getSite3() != t2.getSite3() )))))<<std::endl;
  //   }  
      // if(fabs(t1.getDistance1()-t1.getDistance3())<t1.getLimit())
      // 	{
      // 	  std::cout<<"True"<<std::endl;
      // 	}
      

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
