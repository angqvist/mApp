#include "PairList.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "Pair.hpp"
#include "LatticeList.hpp"


PairList::PairList()
{
  nbrOfPairs=0;
}

int PairList::getNbrOfPairs()
{
  return nbrOfPairs;
}

int PairList::updatePair(Pair newPair,bool add)
{
  for (size_t i=0; i< pairList.size(); i++)
    {
      if(pairList[i]==newPair)
	{
	  pairList[i].incrementCount();
	  return false;
	}
    }
  if(add)
    {
      pairList.push_back(newPair);
      nbrOfPairs++;
    }
  return true;
}
  

int PairList::isAtomInSubElements(std::string atom, std::vector<std::string> subelements)
{
  for( size_t i=0; i< subelements.size(); i++)
    {
      if(atom==subelements[i])
	{
	  return true;
	}
    }
  std::cout<<"Atom type: "<<atom<< " not in subelements"<<std::endl;
  return false;
}



void PairList::initializePairs(LatticeList ll, std::vector<std::string> subelements, double cutOff)
{
  double dr;
  Pair tempPair=Pair();
  for(size_t i=0; i< ll.getNbrOfSites(); i++)
    {
      if(!(isAtomInSubElements(ll.getSite(i),subelements)))
	{
	  continue;
	}
      for(size_t j=i; j<ll.getNbrOfSites(); j++)
	{
	  if(i==j)
	    {
	      continue;
	    }
	  if(!(isAtomInSubElements(ll.getSite(j),subelements)))
	    {
	      continue;
	    }
	  dr=ll.getDistance(i,j);
	  if(dr>cutOff)
	    {
	      continue;
	    }
	  for(size_t k=0; k< subelements.size(); k++)
	    {
	      for(size_t l=0; l<subelements.size(); l++)
		{
		  tempPair.setDistance(dr);
		  tempPair.setSite1(subelements[k]);
		  tempPair.setSite2(subelements[l]);
		  updatePair(tempPair,true);		  
		}//end l loop	  
	    }//end for k loop
	}//end for j loop
    }//end i loop
}

void PairList::printList()
{
  for(size_t i=0; i<pairList.size(); i++)
    {
      pairList[i].printPair();
    }
}

void PairList::resetCounts()
{
  for(size_t i=0; i< pairList.size(); i++)
    {
      pairList[i].setCount(0);
    }
}

Pair& PairList::getPair(int index)
{
  return pairList[index];
}

void PairList::divideCountByTwo()
{
    for(size_t i=0; i< pairList.size(); i++)
    {
      pairList[i].setCount(pairList[i].getCount()/2);
    }
}
