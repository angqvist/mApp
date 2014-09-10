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

TripletList::TripletList(LatticeList ll, std::vector<std::string> subelements, double cutOff)
{
  nbrOfTriplets=0;

  initializeTriplets(ll, subelements,  cutOff);

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
	      dr2=ll.getDistance(i,j2);
	      dr3=ll.getDistance(j,j2);
	 
	      if(dr1>cutOff || dr2>cutOff || dr3>cutOff)
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


std::vector<std::vector<double> > TripletList::getUniqueDistances(double cutoff)
{
  
  std::vector<std::vector<double> > distances;
  
  for(int i=0; i<tripletList.size(); i++)
    {
      if(tripletList[i].getDistance3()>cutoff)
	{
	  continue;
	}
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


std::vector<double> TripletList::getClusterVector(std::vector<std::string > subElements,double cutoff)
{

  const double PI = 3.1415926535897932384626;

  sortTripletList(); //importante
  std::vector<std::vector<double> > uniq_dists =getUniqueDistances(cutoff);
  std::vector<double> clusterVector;

  int  Mi=subElements.size();
  // if(Mi=2)
  //   {
  //     clusterVector.resize(4*uniq_dists.size());
  //   }
  // else if(Mi=3)
  //   {
  //     clusterVector.resize(10*uniq_dists.size());
  //   }
  // else if(Mi=4)
  //   {
  //     clusterVector.resize(20*uniq_dists.size());
  //   }
  // else
  //   {
  //     std::cout<<"Error, size of subelements is bad. Size is: "<<Mi<< " expects value between [2,4] "<<std::endl;
  //   }
  int tripletCount;
  double average;
  double tempAverage;
  int tempT=0;
  double tempVal=0;

  int s1;
  int s2;
  int s3;

  // for(int i=0; i<uniq_dists.size(); i++)
  //   {
  //     std::cout<<uniq_dists[i][0]<<" "<<uniq_dists[i][1]<< " "<<uniq_dists[i][2]<<std::endl;
  //   }

  for(int i=0; i<uniq_dists.size(); i++)
    {
      // std::cout<<"========="<<std::endl;
      //    std::cout<<uniq_dists[i][0]<<" "<<uniq_dists[i][1]<< " "<<uniq_dists[i][2]<<std::endl;
      if((uniq_dists[i][2]+1e-4)>cutoff)
      	{
	  // continue;
      	}
      for(int m=2; m<=Mi; m++)
	{
	  for(int t=0; t<m-1; t++)
	    {
	      for(int l=0; l<=t; l++)
		{
		  tempVal=0.0;
		  tempAverage=0.0;


		  for(int j=0; j<tripletList.size(); j++)
		    {
		      if(tripletList[j].getDistance1()<(uniq_dists[i][0]-1e-4))
		      	{
		      	  //continue;
		      	}
		      // if(tripletList[j].getDistance1()>uniq_dists[i][0])
		      // 	{
		      // 	  break; //speeds things up
		      // 	}

		      //      std::cout<<fabs(tripletList[j].getDistance1()-uniq_dists[i][0])<< " "<<fabs(tripletList[j].getDistance2()-uniq_dists[i][1])<< " "<< fabs(tripletList[j].getDistance3()-uniq_dists[i][2])<<std::endl;
		      
		      if( fabs(tripletList[j].getDistance1()-uniq_dists[i][0])<1e-4
			  && fabs(tripletList[j].getDistance2()-uniq_dists[i][1])<1e-4
			  && fabs(tripletList[j].getDistance3()-uniq_dists[i][2])<1e-4)
			{
			  //  tripletList[j].printTriplet();
			  tempVal=0.0;
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

			  tempT=(m/2);
			  if(((m-2)%2==0))
			    {
			      tempVal=-cos(2*PI*s1*tempT/(Mi));
			    }
			  else
			    {
			      tempVal=-sin(2*PI*s1*tempT/(Mi));
			    }
			  tempT=((t+2)/2); //round down aye
				  
			  if((t%2==0))
			    {
			      tempVal*=-cos(2*PI*s2*tempT/(Mi));
			    }
			  else
			    {
			      tempVal*=-sin(2*PI*s2*tempT/(Mi));
			    }	
			  tempT=((l+2)/2);
			  if((l%2==0))
			    {
			      tempVal*=-cos(2*PI*s3*tempT/(Mi));
			    }
			  else
			    {
			      tempVal*=-sin(2*PI*s3*tempT/(Mi));
			    }	
			  //	  std::cout<<tempVal<<" "<<tripletList[i].std::endl;
			  
			  tempAverage +=tripletList[j].getCount()*tempVal;
			}
		    }//end tripletlist loop
		  clusterVector.push_back(tempAverage);
		}//end l loop
	    }//en t loop
	}//end m loop
    }//end uniq dist loop

  return clusterVector;
}		  
		      
		  
      
