#include "PairList.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "Pair.hpp"
#include "LatticeList.hpp"
#include "clust.hpp"
#include "helperFunctions.hpp"

PairList::PairList()
{
  nbrOfPairs=0;
}

int PairList::getNbrOfPairs()
{
  return nbrOfPairs;
}

int PairList::updatePair(Pair &newPair,bool &add)
{
  for (int i=0; i< pairList.size(); i++)
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
  

int PairList::isAtomInSubElements(std::string atom, std::vector<std::string> &subelements)
{
  for(int i=0; i< subelements.size(); i++)
    {
      if(atom==subelements[i])
	{
	  return true;
	}
    }
  std::cout<<"Atom type: "<<atom<< " not in subelements"<<std::endl;
  return false;
}



void PairList::initializePairs(LatticeList &ll, std::vector<std::string> &subelements, double &cutOff)
{
  ll.calculate_lookup_table();
  std::cout<<"elements size:"<<subelements.size()<<std::endl;
  std::cout<<"number of sites "<< ll.getNbrOfSites()<<std::endl;

  for(int k=0;k<subelements.size(); k++)
    {
      std::cout<<subelements[k]<<std::endl;
    }
  double dr;
  bool addPairWhileUpdating=true;
  Pair tempPair=Pair();
  std::vector<double> dist;
  dist.resize(1);
  std::vector<std::vector< std::string> > all_element_combinations;
  std::vector<double> dists;
  for(int i=0; i < ll.get_original_atoms_count(); i++)
    {
      if(!(isAtomInSubElements(ll.getSite(i),subelements)))
	{
	  continue;
	}
      for(int j=i+1; j<ll.getNbrOfSites(); j++)
	{
	  
	  if(!(isAtomInSubElements(ll.getSite(j),subelements)))
	    {
	      continue;
	    }
	  //dr=ll.getDistance(i,j);
	  dr=ll.fast_distance(i,j);
	  if(dr>cutOff)
	    {
	      continue;
	    }
	  //dists=ll.getPeriodicDistance(i,j,cutOff);
	 
	  dist[0]=dr;
	  
	  all_element_combinations = symmetric_cluster_function(dist,subelements);
	  

	  // std::cout<<"d= "<<d<< " distance= "<<dists[d]<<std::endl;
	  // dist[0]=dists[d];
	  tempPair.setDistance(dist[0]);
	  tempPair.setCount(1);
	  for(int k=0; k<all_element_combinations.size(); k++)
	    {
	      tempPair.setSite1(all_element_combinations[k][0]);
	      tempPair.setSite2(all_element_combinations[k][1]);
	      updatePair(tempPair,addPairWhileUpdating);
	    }

	}//end for j loop
    }//end i loop
  sortPairs();
}


void PairList::sortPairs()
{
  Pair tempPair;
  bool switched=true;

  if(pairList.size()==0)
    {
      switched=false;
    }
  while(switched)
    {
      switched=false;
      for(int i=0; i<pairList.size()-1; i++)
	{
	  if(pairList[i].getDistance()>pairList[i+1].getDistance())
	    {
	      tempPair=pairList[i];
	      pairList[i]=pairList[i+1];
	      pairList[i+1]=tempPair;
	      switched=true;
	    }
	}
    }
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
  //death by reference
  int zero = 0;
  for(size_t i=0; i< pairList.size(); i++)
    {
      pairList[i].setCount(zero);
    }
}

Pair& PairList::getPair(int &index)
{
  return pairList[index];
}

void PairList::divideCountByTwo()
{
  int count;
    for(int i=0; i< pairList.size(); i++)
    {
      count = pairList[i].getCount()/2;
      pairList[i].setCount( count);
    }
}

std::vector<double> PairList::getUniqueDistances(double &cutoff)
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
  std::sort (uniq_dist.begin(),uniq_dist.end());
  return uniq_dist;
}
//int s1;
//int s2;


std::vector<double> PairList::getClusterVector(std::vector<std::string> &subElements, double &cutoff, bool &doAverage)
{
  const double PI = 3.1415926535897932384626;

  std::vector<double> uniq_dist = getUniqueDistances(cutoff);  
  // std::sort (uniq_dist.begin(),uniq_dist.end());
  int Mi=subElements.size();
  std::vector<double> clusterVector;
  

  int tempT=0;
  double tempVal=0;
  int totalPairs;
  double tempAverage;
  bool lexicalSort=true;


  std::vector<double> dist;
  dist.resize(1);
  std::vector<std::vector<int> > cluster_functions;
  dist[0]=uniq_dist[0];
  //true for reverse sort
  //same cluster_function for all pairs so can take it out of loop, note that it is not true for triplets and so on
  int atom1;
  int atom2;
  cluster_functions = symmetric_cluster_function(dist,Mi,lexicalSort);
  std::vector<std::vector<std::string > > all_element_combinations= symmetric_cluster_function(dist,subElements);
  
  for(int i=0; i<uniq_dist.size(); i++)
    {
      
      for(int j=0; j<cluster_functions.size(); j++)
	{
	  tempAverage=0;
	  totalPairs=0;
	  for(int jj2=0; jj2<pairList.size(); jj2++)
	    {
	      
	      if(fabs(pairList[jj2].getDistance()-uniq_dist[i])<1e-4)
		{
		  for(int jj=0; jj<subElements.size(); jj++)
		    {
		      if(pairList[jj2].getSite1()==subElements[jj])
			{
			  atom1=jj;
			}
		      if(pairList[jj2].getSite2()==subElements[jj])
			{
			  atom2=jj;
			}

		    }
		  
		  
		  //	  pairList[jj2].printPair();

		  //		  std::cout<<"s1,s2,weight: S1 "<<atom1<< "  s2 "<<atom2<< " clusterFunction1 "<<cluster_functions[j][0]<< " clustFunc2 "<<cluster_functions[j][1]<< " weight "<< (clusterFunction(Mi,atom1,cluster_functions[j][0]))*(clusterFunction(Mi,atom2,cluster_functions[j][1]))<<std::endl;;

		  //  std::cout<<cluster_functions[j][0]<< " "<< Mi<< " "<<atom1<< " "<< clusterFunction(Mi,atom1,cluster_functions[j][0])<<std::endl;
		  tempAverage += clusterFunction(Mi,atom1,cluster_functions[j][0])*
		    clusterFunction(Mi,atom2,cluster_functions[j][1])*pairList[jj2].getCount();

		  //std::cout<<"COUNT "<<pairList[jj2].getCount()<<std::endl;
		  //std::cout<<"hoho"<<std::endl;
		  //	  std::cout<<clusterFunction(Mi,atom1,cluster_functions[j][0])*
		  // clusterFunction(Mi,atom2,cluster_functions[j][1])<<std::endl;
		  
		  //  std::cout<<std::endl;
		  totalPairs += pairList[jj2].getCount();
		}
	    }
	  if(!doAverage)
	    {
	      //  std::cout<< "j= "<< j<< " tempAverage = "<< tempAverage<<std::endl;
	      clusterVector.push_back(tempAverage);
	    }
	  else
	    {
	      if(totalPairs==0)
		{
		  std::cout<<" tempAverage " <<tempAverage<< "  division: "<< tempAverage/(double)totalPairs<<std::endl;
		  pairList[i*3].printPair();
		  pairList[i*3+1].printPair();
		  pairList[i*3+2].printPair();
		  
		}
	      clusterVector.push_back(tempAverage/(double)totalPairs);
	    }
	}
    
    }


  return clusterVector;
}
