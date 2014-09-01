#include "NeighbourList.hpp"
#include "Neighbour.hpp"
#include "LatticeList.hpp"
#include "ParameterList.hpp"
#include "Pair.hpp"
#include <vector>
#include <iostream>
#include <cmath>


// should really not be used with empty constructor
NeighbourList::NeighbourList()
{
  thisIndex=0;
  distanceLimit=0.001;
  // empty
}
//check up if you should send references to latticeList or something like that
NeighbourList::NeighbourList(int i, LatticeList ll, ParameterList pl)
{
  thisIndex=i;
  distanceLimit=0.001;
  findNeighbours(ll,pl);
}
void NeighbourList::addNeighbour(Neighbour nbr1)
{
  if(newNeighbour(nbr1))
    {
      nbrList.push_back(nbr1);
    }
}
int NeighbourList::newNeighbour(Neighbour nbr1)
{
  for(size_t i=0; i < nbrList.size(); i++)
    {
      if(nbr1==nbrList[i])
      {
	// std::cout<<nbrList.size()<<std::endl;
	// std::cout<<"Tried to add neighbour which already existed in neighbourList."<<std::endl<<"Neighbour tried to add: ";
	// nbr1.print();
	// std::cout<<"Index: "<<i<<std::endl;
	return false;
      }
    }
  return true;
}
void NeighbourList::printList()
{
  for(size_t i=0; i < manyNbrLists.size(); i++)
    {
      std::cout<<indexList[i]<<": ";
      printManyNbr(manyNbrLists[i]);
      std::cout<<"--------------------------"<<std::endl;
      std::cout<<" "<<std::endl;
    }
}

void NeighbourList::printManyNbr(std::vector<Neighbour> nl)
{
  for(size_t i=0; i < nl.size(); i++)
    {
      nl[i].print();
    }
}
void NeighbourList::findNeighbours(LatticeList ll, ParameterList pl)
{
  double tempDist;
  bool paramAtIndex;
  int singletIndex=0;
  while(pl.getPair(singletIndex).getDistance()<distanceLimit)
    {
      //std::cout<<"index: "<<singletIndex<<" dist: "<<pl.getPair(singletIndex).getDistance()<< " type "<<pl.getPair(singletIndex).getSite1()<< " energy "<<pl.getPair(singletIndex).getEnergy()<<std::endl;
      singletType.push_back(pl.getPair(singletIndex).getSite1());
      singletEnergy.push_back(pl.getPair(singletIndex).getEnergy());
      singletIndex++;
    }
  //std::cout<<"singlet index="<<singletIndex<<std::endl;
 


  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      paramAtIndex=false;
      if(i == thisIndex)
	{
	  continue;
	}
      tempDist=ll.getDistance(i,thisIndex);
      for(int j=0; j<pl.getNbrOfParams(); j++)
	{
	  if(fabs(tempDist-pl.getPair(j).getDistance())<distanceLimit && pl.getPair(j).getDistance()>0.00001)
	    {
	      Neighbour tempNbr = Neighbour(pl.getPair(j).getSite1(),pl.getPair(j).getSite2(),pl.getPair(j).getEnergy(),pl.getPair(j).getDistance());
	      if(newNeighbour(tempNbr) )
		{
		  nbrList.push_back(tempNbr);
		  paramAtIndex=true;		      
		}
	    }
	}
	
      if(!(nbrList.empty()))
	{
	  manyNbrLists.push_back(nbrList);
	}
      nbrList.clear();
      if(paramAtIndex)
	{
	  indexList.push_back(i);
	}
    }// end i loop
}







double  NeighbourList::getMatchingEnergy(Neighbour nbr, std::vector<Neighbour> nl)
{
  for(size_t j=0; j<nl.size(); j++)
    {
      //std::cout<<"tempNbr: "<<std::endl;
      //nbr.print();
      //std::cout<<"neigbourlist  "<<std::endl;
      //nl[j].print();

      if(nbr==nl[j])
	{
	  return nl[j].getEnergy();
	}
    }
  return 0;
}

double NeighbourList::getLocalEnergy(LatticeList ll)
{

  localEnergy=0.0;
  //add singlet Energy times two?

  for(int i=0; i<singletEnergy.size(); i++)
    {
      if(ll.getSite(thisIndex)==singletType[i])
  	{
  	  localEnergy+=singletEnergy[i]*2.0;
  	}
    }
  Neighbour tempNbr= Neighbour();
  tempNbr.setSite1(ll.getSite(thisIndex));  
  for(size_t i=0; i < indexList.size(); i++)
    {
      tempNbr.setSite2(ll.getSite(indexList[i]));
      localEnergy += getMatchingEnergy(tempNbr,manyNbrLists[i]);
    }
  return localEnergy;
}




int NeighbourList::isThisMatchingNeighbour(LatticeList ll, Neighbour n1)
{
  Neighbour tempNbr= Neighbour();
  tempNbr.setSite1(ll.getSite(thisIndex));
  int returnNumber=0;
  for(size_t i=0; i<indexList.size(); i++)
    {
      tempNbr.setDistance(ll.getDistance(thisIndex,indexList[i]));
      tempNbr.setSite2(ll.getSite(indexList[i]));			  
      if(tempNbr==n1 && fabs(tempNbr.getDistance()-n1.getDistance())<1e-3)
	{
	  // std::cout<<"==================== "<< thisIndex<< " "<<indexList[i]<<std::endl;
	  // tempNbr.print();
	  // n1.print();
	  returnNumber++;
	}
    }
  return returnNumber;
}
      
      
double NeighbourList::getCurrentLocalEnergy()
{
  return currentLocalEnergy;
}

void NeighbourList::setCurrentLocalEnergy(double newEnergy)
{
 currentLocalEnergy=newEnergy;
}

void NeighbourList::calcCurrentLocalEnergy(LatticeList ll)
{
  currentLocalEnergy = getLocalEnergy(ll);
}
