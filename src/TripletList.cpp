#include "TripletList.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "Triplet.hpp"
#include "LatticeList.hpp"
#include <cmath>

TripletList::TripletList()
{
  nbrOfTriplets=0;
}

int TripletList::getNbrOfTriplets()
{
  return nbrOfTriplets;
}

int TripletList::updateTriplet(Triplet newTriplet,bool add)
{
    for (size_t i=0; i< tripletList.size(); i++)
	 {	
	    if(tripletList[i]==newTriplet) 
		{
			tripletList[i].incrementCount();
			 return false;
		}
    }
  if(add)
    {
      tripletList.push_back(newTriplet);
      nbrOfTriplets++;
    }
  return true;
}

int TripletList::isAtomInSubElements(std::string atom, std::vector<std::string> subelements)
{
  for( size_t i=0; i< subelements.size(); i++)
    {
      if(atom==subelements[i])
	{
	  return true;
	}
    }
  std::cout<<"Atom type: "<<atom<< " not in subelements(in triplet class)"<<std::endl;
  return false;
}





void TripletList::initializeTriplets(LatticeList ll, std::vector<std::string> subelements, double cutOff)
{
  double dr1;
  double dr2;
  double dr3;
  Triplet tempTriplet=Triplet();
  std::vector<double> orderDr;
  std::vector<int> orderIndex;
  orderDr.resize(3);
  orderIndex.resize(3);
  for(size_t i=0; i< ll.getNbrOfSites(); i++)//first loop
    {
      if(!(isAtomInSubElements(ll.getSite(i),subelements)))
	{
	  std::cout<<"Atom not in subelements in TripletList"<<std::endl;
	  continue;
	}
      for(size_t j=i; j<ll.getNbrOfSites(); j++) // second loop
	{
	  if(j==i)
	    {
	      continue;
	    }
	  if(!(isAtomInSubElements(ll.getSite(j),subelements)))
	    {
	      std::cout<<"Atom not in subelements in TripletList"<<std::endl;
	      continue;
	    }
	  for(size_t j2=j; j2<ll.getNbrOfSites(); j2++) // third loop
	    {

	      if(!(isAtomInSubElements(ll.getSite(j2),subelements)))
		{
		  std::cout<<"Atom not in subelements in TripletList"<<std::endl;
		  continue;
		}
	      if(i==j2 || j == j2)
		{
		  //std::cout<<"ofta eller?"<<std::endl;
		  continue;
		}
	      
	      dr1=ll.getDistance(i,j);
	      dr2=ll.getDistance(j,j2);
	      dr3=ll.getDistance(i,j2);
	 
	      if(dr1>cutOff || dr2>cutOff)
		{
		  continue;
		}
	      //std::cout<<"over cutoff "<<dr1<<" "<<dr2<< " "<<dr3<<std::endl;

	      orderDr[0]=dr1;
	      orderDr[1]=dr2;
	      orderDr[2]=dr3;
	      orderIndex[0]=i;
	      orderIndex[1]=j;
	      orderIndex[2]=j2;
	       //std::cout<<"Before "<<orderDr[0]<< " "<<orderDr[1]<< " "<<orderDr[2]<<std::endl;
	      sortOrder(orderDr,orderIndex);
	       // std::cout<<"After "<<orderDr[0]<< " "<<orderDr[1]<< " "<<orderDr[2]<<std::endl;

	      for(int k=0; k< subelements.size(); k++)
		{
		  for(int l=0; l<subelements.size(); l++)
		    {
		      for(int p=0; p<subelements.size(); p++)
			{
			  tempTriplet.setDistance1(orderDr[0]);
			  tempTriplet.setDistance2(orderDr[1]);
			  tempTriplet.setDistance3(orderDr[2]);
			  tempTriplet.setSite1(subelements[k]);
			  tempTriplet.setSite2(subelements[l]);
			  tempTriplet.setSite3(subelements[p]);
			  tempTriplet.sortTriplet();
			  updateTriplet(tempTriplet,true);		  
			}//end p loop
		    }//end l loop	  
		}//end for k loop
	    }//end j2 loop
	}//end for j loop
    }//end i loop
}

void TripletList::printList()
{
  for(size_t i=0; i<tripletList.size(); i++)
    {
      tripletList[i].printTriplet();
    }
}

void TripletList::resetCounts()
{
  for(size_t i=0; i< tripletList.size(); i++)
    {
      tripletList[i].setCount(0);
    }
}

Triplet& TripletList::getTriplet(int index)
{
  return tripletList[index];
}

void TripletList::sortOrder(std::vector<double>  &orderDr, std::vector<int>  &orderIndex)
{
  bool swapped=true;
  double tempDr;
  int tempIndex;
  while(swapped)
    {
      swapped=false;
      for(int i=0; i<2; i++)
	{
	  if(orderDr[i]>orderDr[i+1])
	    {
	      tempDr=orderDr[i];
	      tempIndex=orderIndex[i];
	      orderDr[i]=orderDr[i+1];
	      orderDr[i+1]=tempDr;
	      orderIndex[i]=orderIndex[i+1];
	      orderIndex[i+1]=tempIndex;
	      swapped=true;
	    }
	}
    }
}


std::vector<std::vector<double> > TripletList::getUniqueDistances()
{
  
  std::vector<std::vector<double> > distances;
  
  for(int i=0; i<tripletList.size(); i++)
    {
      // tripletList[i].printTriplet();
      isTripletUnique(tripletList[i],distances,true);

    }
  


  return distances;
}

bool TripletList::isTripletUnique(Triplet t1, std::vector<std::vector<double> > &distances, bool add)
{
  
  if(distances.size()==0)
    {
      if(add)
	{
	  std::vector<double> addThisToDistances;
	  addThisToDistances.push_back(t1.getDistance1());
	  addThisToDistances.push_back(t1.getDistance2());
	  addThisToDistances.push_back(t1.getDistance3());
	  distances.push_back(addThisToDistances);
	}
      return true;
    }

  for(int i=0; i<distances.size(); i++)
    {
      if((fabs(distances[i][0]-t1.getDistance1())<1e-3)
	 && (fabs(distances[i][1]-t1.getDistance2())<1e-3)
	 && (fabs(distances[i][2]-t1.getDistance3())<1e-3))
	{
	  return false;
	} 
  
    }

    if(add)
	{
	  std::vector<double> addThisToDistances;
	  addThisToDistances.push_back(t1.getDistance1());
	  addThisToDistances.push_back(t1.getDistance2());
	  addThisToDistances.push_back(t1.getDistance3());
	  distances.push_back(addThisToDistances);
	}
      return true;

}


  
  
