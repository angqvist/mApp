#include "TripletList.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "Triplet.hpp"
#include "LatticeList.hpp"
#include <cmath>
#include "clust.hpp"
TripletList::TripletList()
{
  nbrOfTriplets=0;
}

TripletList::TripletList(LatticeList &ll, std::vector<std::string> &subelements, double &cutOff)
{
  std::cout<<"initialize_trips"<<std::endl;

  nbrOfTriplets=0;
  std::cout<<"initialize_trips"<<std::endl;
  initializeTriplets(ll, subelements,  cutOff);

}



int TripletList::getNbrOfTriplets()
{
  return nbrOfTriplets;
}

int TripletList::updateTriplet(Triplet &newTriplet,bool &add,int &multiplicity)
{
  //newTriplet.sortTriplet();
  for (size_t i=0; i< tripletList.size(); i++)
    {	
      if(tripletList[i]==newTriplet) 
	{
	  tripletList[i].increaseCountBy(multiplicity);
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

void TripletList::reset()
{
  tripletList.clear();
  nbrOfTriplets=0;
}

int TripletList::isAtomInSubElements(std::string atom, std::vector<std::string> &subelements)
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





void TripletList::initializeTriplets(LatticeList &ll, std::vector<std::string> &subelements, double &cutOff)
{
  std::cout<<"initialize lookuptable in triplets"<<std::endl;
  ll.calculate_lookup_table();
  std::cout<<"done"<<std::endl;
  reset();
  bool addTripletWhenUpdating=true;
  int countForTripletWhenUpdating=1;
  bool lexicalSort=false;
  double dr1;
  double dr2;
  double dr3;
  Triplet tempTriplet=Triplet();
  std::vector<double> orderDr;
  std::vector<int> orderIndex;
  orderDr.resize(3);
  orderIndex.resize(3);
  std::vector<std::string> elements;
  elements.resize(3);
  std::vector<std::vector< std::string> > all_element_combinations;

  for(int i=0; i< ll.get_original_atoms_count(); i++)//first loop
    {
      if(!(isAtomInSubElements(ll.getSite(i),subelements)))
	{
	  std::cout<<"Atom not in subelements in TripletList"<<std::endl;
	  continue;
	}
      for(int j=i+1; j<ll.getNbrOfSites(); j++) // second loop
	{
	  if(!(isAtomInSubElements(ll.getSite(j),subelements)))
	    {
	      std::cout<<"Atom not in subelements in TripletList"<<std::endl;
	      continue;
	    }
	  for(int j2=j+1; j2<ll.getNbrOfSites(); j2++) // third loop
	    {

	      if(!(isAtomInSubElements(ll.getSite(j2),subelements)))
		{
		  std::cout<<"Atom not in subelements in TripletList"<<std::endl;
		  continue;
		}
	      
	      dr1=ll.fast_distance(i,j);
	      dr2=ll.fast_distance(i,j2);
	      dr3=ll.fast_distance(j,j2);
	 

	      //std::cout<<"over cutoff "<<dr1<<" "<<dr2<< " "<<dr3<<std::endl;

	      orderDr[0]=dr1;
	      orderDr[1]=dr2;
	      orderDr[2]=dr3;
	      // orderIndex[0]=i;
	      // orderIndex[1]=j;
	      // orderIndex[2]=j2;

	      elements[0]=ll.getSite(i);
	      elements[1]=ll.getSite(j);
	      elements[2]=ll.getSite(j2);
	      
	      // all_element_combinations=symmetric_cluster_function
	      //false for sorting alphabetically

	      
	      //tuple_remodulator(orderDr,elements,false);
	      


	      
	       //std::cout<<"Before "<<orderDr[0]<< " "<<orderDr[1]<< " "<<orderDr[2]<<std::endl;
	      // sortOrder(orderDr,orderIndex);
	      if(orderDr[0]>cutOff || orderDr[1]>cutOff || orderDr[2]>cutOff)
		{
		  continue;
		}
	       // std::cout<<"After "<<orderDr[0]<< " "<<orderDr[1]<< " "<<orderDr[2]<<std::endl;

	      //gives all remodulated, symmetric distinct triplets for a given dists
	      all_element_combinations= symmetric_cluster_function(orderDr,subelements); //DOESNT CHANGE DIST ORDER 

	      
	      std::vector<std::string> ghost_atoms;
	      ghost_atoms.push_back("A");
	      ghost_atoms.push_back("A");
	      ghost_atoms.push_back("A");
	      //   std::cout<<orderDr[0]<< " "<<orderDr[1]<< " "<<orderDr[2]<<std::endl;
	      tuple_remodulator(orderDr,ghost_atoms,lexicalSort);
	      //    std::cout<<orderDr[0]<< " "<<orderDr[1]<< " "<<orderDr[2]<<std::endl;
	      
	      // for(int kk=0; kk<3; kk++)
	      // 	{
	      // 	  std::cout<<orderDr[kk]<< " ";
		
	      // 	}
	      //std::cout<<std::endl;
	    
	      tempTriplet.setDistance1(orderDr[0]);
	      tempTriplet.setDistance2(orderDr[1]);
	      tempTriplet.setDistance3(orderDr[2]);

	      
	      for(int k=0; k<all_element_combinations.size(); k++)
		{

		  // std::cout<<all_element_combinations[k][0]<< " "<<all_element_combinations[k][1] << " "<<all_element_combinations[k][2]<<std::endl;
		  tempTriplet.setSite1(all_element_combinations[k][0]);
		  tempTriplet.setSite2(all_element_combinations[k][1]);
		  tempTriplet.setSite3(all_element_combinations[k][2]);
		  updateTriplet(tempTriplet,addTripletWhenUpdating,countForTripletWhenUpdating);		  

		}
		  


	    }//end j2 loop
	}//end for j loop
    }//end i loop
}

void TripletList::printList()
{
  for(int i=0; i<tripletList.size(); i++)
    {
      tripletList[i].printTriplet();
    }
}

void TripletList::resetCounts()
{
  int zero=0;
  for(int i=0; i< tripletList.size(); i++)
    {
      tripletList[i].setCount(zero);
    }
}

Triplet& TripletList::getTriplet(int &index)
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
	      if(i==0)
		{
		  //swap atoms 2 and 3
		  tempIndex=orderIndex[1];
		  orderIndex[1]=orderIndex[2];
		  orderIndex[2]=tempIndex;
		}
	      if(i==1)
		{
		  //swap atoms 1 and 2
		  tempIndex=orderIndex[0];
		  orderIndex[0]=orderIndex[1];
		  orderIndex[1]=tempIndex;
		}

		
	      tempDr=orderDr[i];
	      orderDr[i]=orderDr[i+1];
	      orderDr[i+1]=tempDr;	      
	      swapped=true;
	    }
	}
    }
}


std::vector<std::vector<double> > TripletList::getUniqueDistances(double &cutoff)
{
  
  std::vector<std::vector<double> > distances;
  bool add=true;
  for(int i=0; i<tripletList.size(); i++)
    {
      if(tripletList[i].getDistance2()>cutoff)
	{
	  continue;
	}
      // tripletList[i].printTriplet();
      isTripletUnique(tripletList[i],distances,add);
    }
  clust_sort_dists(distances);
  return distances;
}



bool TripletList::isTripletUnique(Triplet &t1, std::vector<std::vector<double> > &distances, bool &add)
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
      if((fabs(distances[i][0]-t1.getDistance1())<1e-4)
	 && (fabs(distances[i][1]-t1.getDistance2())<1e-4)
	 && (fabs(distances[i][2]-t1.getDistance3())<1e-4))
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


  
  
void TripletList::sortTripletList()
{
  bool swapped=true;
  Triplet tempTriplet;
  
  while(swapped)
    {
      swapped=false;
      for(int i=0; i<tripletList.size()-1; i++)
	{
	  if(tripletList[i].getDistance1()>tripletList[i+1].getDistance1())
	    {
	      tempTriplet=tripletList[i];
	      tripletList[i]=tripletList[i+1];
	      tripletList[i+1]=tempTriplet;
	      swapped=true;
	    }
	}
    }
  swapped=true;

  while(swapped)
    {
      swapped=false;
      for(int i=0; i<tripletList.size()-1; i++)
	{

	   if(fabs(tripletList[i].getDistance1()-tripletList[i+1].getDistance1())<1e-4
	  	  && tripletList[i].getDistance2()>tripletList[i+1].getDistance2())
	    {
	      tempTriplet=tripletList[i];
	      tripletList[i]=tripletList[i+1];
	      tripletList[i+1]=tempTriplet;
	      swapped=true;
	    }

	}

    }
  swapped=true;
  while(swapped)
    {
      swapped=false;
      for(int i=0; i<tripletList.size()-1; i++)
	{

	  if(fabs(tripletList[i].getDistance1()-tripletList[i+1].getDistance1())<1e-4
	  	  && fabs(tripletList[i].getDistance2()-tripletList[i+1].getDistance2())<1e-4
	  	  && tripletList[i].getDistance3()>tripletList[i+1].getDistance3())
	    {
	      tempTriplet=tripletList[i];
	      tripletList[i]=tripletList[i+1];
	      tripletList[i+1]=tempTriplet;
	      swapped=true;
	    }

	}
    }
}


std::vector<double> TripletList::getClusterVector(std::vector<std::string > &subElements,double &cutoff,bool &average)
{

  const double PI = 3.1415926535897932384626;
  bool lexicalSort=true;
  sortTripletList(); //importante
  std::vector<std::vector<double> > uniq_dists = getUniqueDistances(cutoff);
  clust_sort_dists(uniq_dists);

  std::cout<<"--------------------------"<<std::endl;
  for(int i=0; i<uniq_dists.size(); i++)
    {
      for(int j=0;j<3; j++)
	{
	  std::cout<< uniq_dists[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"--------------------------"<<std::endl;

  std::vector<double> clusterVector;

  int  Mi=subElements.size();
  int tripletCount;
  double tempAverage;
  int tempT=0;
  double tempVal=0;
  int tempTripletCount;

  int s1;
  int s2;
  int s3;



  std::vector<double> dist;
  dist.resize(3);
  std::vector<std::vector<int> > clusterFunctions;


  //new stuff, shiny and exciting
  for(int i=0; i<uniq_dists.size(); i++)
    {
      clusterFunctions = symmetric_cluster_function(uniq_dists[i],Mi,lexicalSort);
      

      for(int ii2=0; ii2<clusterFunctions.size(); ii2++)
	{
	  tempAverage=0.0;
	  tempTripletCount=0;

	  for(int j=0; j<tripletList.size(); j++)
	    {
	      
	  
	      if( fabs(tripletList[j].getDistance1()-uniq_dists[i][0])<1e-4
		  && fabs(tripletList[j].getDistance2()-uniq_dists[i][1])<1e-4
		  && fabs(tripletList[j].getDistance3()-uniq_dists[i][2])<1e-4)
		{
		  for(int ii=0; ii<subElements.size(); ii++)
		    {
		      if(subElements[ii]==tripletList[j].getSite1())
			{
			  s1=ii;
			}
		      if(subElements[ii]==tripletList[j].getSite2())
			{
			  s2=ii;
			}
		      if(subElements[ii]==tripletList[j].getSite3())
			{
			  s3=ii;
			}
		    }
		  //  std::cout<<"TRIP count "<<tripletList[j].getCount()<<std::endl;
		  //	  std::cout<<  clusterFunction(Mi,s1,clusterFunctions[ii2][0])*
		  //  clusterFunction(Mi,s2,clusterFunctions[ii2][1])*
		  //    clusterFunction(Mi,s3,clusterFunctions[ii2][2])*tripletList[j].getCount()<<std::endl;

		  tempAverage += clusterFunction(Mi,s1,clusterFunctions[ii2][0])*
		    clusterFunction(Mi,s2,clusterFunctions[ii2][1])*
		    clusterFunction(Mi,s3,clusterFunctions[ii2][2])*tripletList[j].getCount();
		  tempTripletCount += tripletList[j].getCount();
		}
	    }
	  if(average)
	    {		  
	      clusterVector.push_back(tempAverage/(double)tempTripletCount);
	    }
	  else
	    {
	      clusterVector.push_back(tempAverage);
	    }
	}
    

    }//end new stuff  


  return clusterVector;
}		  
		      
		  
      

