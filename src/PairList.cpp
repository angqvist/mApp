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

std::vector<double> PairList::getUniqueDistances(double cutoff)
{
  std::vector<double> uniq_dist;
  bool addDistance;
  for(int i=0; i<pairList.size(); i++)
    {
      addDistance=true;
      if(pairList[i].getDistance()>cutoff)
	{
	  continue;
	}
      for(int j=0; j<uniq_dist.size(); j++)
	{
	  if(fabs(pairList[i].getDistance()-uniq_dist[j])<1e-4)
	    {
	      addDistance=false;
	      break;
	    }
	}
      if(addDistance)
	uniq_dist.push_back(pairList[i].getDistance());
    }
  return uniq_dist;
}
int s1;
int s2;


std::vector<double> PairList::getClusterVector(std::vector<std::string> subElements, double cutoff,bool doAverage)
{
  const double PI = 3.1415926535897932384626;

  std::vector<double> uniq_dist = getUniqueDistances(cutoff);  
  int Mi=subElements.size();
  std::vector<double> clusterVector;
  
  //first do singlets---  no! do it somewhere else... you need latticeList for it

  int tempT=0;
  double tempVal=0;
  int totalPairs;
  double tempAverage;

      
  for(int i=0; i<uniq_dist.size(); i++)
    {
      for(int m=2; m<=Mi; m++)
	{	  
	  for(int t=0; t<m-1; t++)
	    {
	      tempAverage=0;
	      totalPairs=0;
	      for(int j=0; j<pairList.size(); j++)
		{		  
		  if(fabs(pairList[j].getDistance()-uniq_dist[i])<1e-4)
		    {
		      
		      for(int ii=0; ii<Mi; ii++)
			{
			  if(subElements[ii]==pairList[j].getSite1())
			    {
			      s1=ii;
			    }
		   
			  if(subElements[ii]==pairList[j].getSite2())
			    {
			      s2=ii;
			    }
			}
		      tempT=(m/2); //round down aye
		      if(((m-2)%2==0))
			{
			  tempVal=-cos(2*PI*s1*tempT/(subElements.size()));
			}
		      else
			{
			  tempVal=-sin(2*PI*s1*tempT/(subElements.size()));
			}
		      tempT=((t+2)/2); //round down aye
				  
		      if((t%2==0))
			{
			  tempVal*=-cos(2*PI*s2*tempT/(subElements.size()));
			}
		      else
			{
			  tempVal*=-sin(2*PI*s2*tempT/(subElements.size()));
			}	
		      totalPairs += pairList[j].getCount();
		      tempAverage +=pairList[j].getCount()*tempVal;
		    }
		}
	      if( !doAverage )
		{
		  clusterVector.push_back(tempAverage);
		}
	      else
		{
		  clusterVector.push_back(tempAverage/(double)totalPairs);
		}
	    }
	}
    }
  return clusterVector;
}
