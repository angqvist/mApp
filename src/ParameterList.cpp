#include "ParameterList.hpp"
#include <iostream>
#include <fstream>
#include "Pair.hpp"
#include "PairList.hpp"
#include "clust.hpp"
#include <cmath>

// // // // // #include "Atom.hpp"
// // #include "Triplet.hpp"
// #include "Quatuplet.hpp"

ParameterList::ParameterList()
{
  fileName="params.data";
  nbrOfParams=-1;
  eCutOff=0.001;
  readParams();
}

ParameterList::ParameterList(std::string newFileName,double newCutOff)
{
  fileName=newFileName;
  nbrOfParams=-1;
  eCutOff=newCutOff;
  readParams();
}

ParameterList::ParameterList(PairList pl)
{
  fileName="";
  nbrOfParams=0;
  eCutOff=-0.01;
  readParamsWithPL(pl);
}
ParameterList::ParameterList(std::string newFileName,double newCutOff,std::vector<std::string > subElements)
{
  fileName=newFileName;
  eCutOff=newCutOff;
  nbrOfParams=0;

  //readParamsATATStyle(subElements);
  readParams_new(subElements);
  unwrapSinglets(paramList_0,subElements);
  unwrapPairs(paramList,subElements);
  unwrapTriplets(paramList_3,subElements);
  unwrapQuatuplets(paramList_4,subElements);
  std::cout<<"Done unwrapping"<<std::endl;
  std::cout<<" Number of pairs: "<< paramList.size()<<" Triplets: "<<paramList_3.size()<< " Quatuplets: "<< paramList_4.size()<<std::endl;
  for(int i=0; i<paramList_4.size(); i++)
    {
      paramList_4[i].print();
    }

}


void ParameterList::readParamsWithPL(PairList pl)
{
  for(int i=0; i<pl.getNbrOfPairs(); i++)
    {
      paramList.push_back(pl.getPair(i));
      nbrOfParams++;
    }
}

int ParameterList::getNbrOfPairs()
{
  // if(nbrOfParams==-1)
  //   {
  //     std::cout<<" No params has been read, nbr of params = -1"<<std::endl;
  //   }
  return paramList.size();
}


void ParameterList::readParamsATATStyle(std::vector<std::string > subElements)
{

const double PI = 3.1415926535897932384626;

  std::ifstream in(fileName.c_str());  
  if (!in)
    {
      std::cout<< "Cannot open file: "<<fileName<<std::endl;
      return;
    }

  int Mi=subElements.size();

  int nbrSinglet;
  int nbrPairs;
  if(Mi==2)
    {
      nbrSinglet=1;
      nbrPairs=1;
    }

  else if(Mi==3)
    {
      nbrSinglet=2;
      nbrPairs=3;
    }

  else if(Mi==4)
    {
      nbrSinglet=3;
      nbrPairs=4;
    }

  else if(Mi==5)
    {
      nbrSinglet=4;
      nbrPairs=5;
    }
  
  

  std::vector<double> distances;
  std::vector<double> energies;
  double tempValue;
  while(!(in.eof()))
    {
      in >> tempValue;
      distances.push_back(tempValue);
      in >> tempValue;
      energies.push_back(tempValue);
      // std::cout<<tempValue<<std::endl;
    }
  int count=0;
  Pair tempPair = Pair();
  double tempEnergy=0;
  int tempT;
  for(int i=0; i<Mi; i++)
    {
      count=0;
      tempPair.setSite1(subElements[i]);
      tempPair.setSite2(subElements[i]);
      tempPair.setDistance(0.0);
      for(int m=2; m<=Mi; m++)
	{
	  tempT=(m/2); //round down aye
	  if((((m-2)%2)==0))
	    {
	      tempEnergy += -cos(2*PI*i*tempT/Mi)*energies[count];
	    }

	  else
	    {
	      tempEnergy += -sin(2*PI*i*tempT/Mi)*energies[count];
	    }
	  count++;
	}
      tempPair.setEnergy(tempEnergy);
      paramList.push_back(tempPair);	      
      nbrOfParams++;
      tempEnergy=0;
    }
      
				

 
  for(int i=nbrSinglet; i<distances.size(); i +=nbrPairs)
    {
      //std::cout<< i<< " "<<distances[i]<<std::endl;
      for(int ii=0; ii<Mi; ii++)
	{
	  for(int jj=ii; jj<Mi; jj++)
	    {
	      // std::cout<<distances[i]<<" "<< subElements[ii]<< " "<<subElements[jj]<<std::endl;
	      
	      tempPair.setSite1(subElements[ii]);
	      tempPair.setSite2(subElements[jj]);
	      tempPair.setDistance(distances[i]);
	      tempPair.setEnergy(0);
	      tempEnergy=0;
	      count=0;
	      for(int m=2; m<=Mi; m++)
		{
		  for(int t=0; t<m-1; t++)
		    {
		      //std::cout<<count<<std::endl;
		      tempT=(m/2); //round down aye
		      if((((m-2)%2)==0))
			{
			  tempEnergy = -cos(2*PI*ii*tempT/Mi);
			}
		      else
			{
			  tempEnergy = -sin(2*PI*ii*tempT/Mi);
			}

		      tempT=((t+2)/2); //round down aye
				  
		      if(((t%2)==0))
			{
			  tempEnergy *= -cos(2*PI*jj*tempT/Mi);
			}
		      else
			{
			  tempEnergy *= -sin(2*PI*jj*tempT/Mi);
			}	
		      //std::cout<<distances[i]<< " "<<distances[i-1]<< " "<<energies[i+count]<<std::endl;
		      tempEnergy = tempPair.getEnergy()+tempEnergy*energies[i+count];
		      tempPair.setEnergy(tempEnergy);
		      count++;


		    }
		}// end m
	      if(fabs(tempPair.getEnergy())>eCutOff)
		{
		  paramList.push_back(tempPair);
		  nbrOfParams++;
		}

		  
	    }
	}
    }
   
}

void ParameterList::readParams()
{
  std::ifstream in(fileName.c_str());  
  if (!in)
    {
      std::cout<< "Cannot open file: "<<fileName<<std::endl;
      return;
    }
  Pair tempPair= Pair();
  std::string tempS1;
  std::string tempS2;
  double tempDist;
  double tempEnergy;
  int tempNbrParams=0;
  while(!(in.eof()))
    {
      //if( in.eof() ) break;

      in>> tempDist;
      in>> tempS1;
      in>> tempS2;
      in>> tempEnergy;
      tempPair.setDistance(tempDist);
      tempPair.setSite1(tempS1);
      tempPair.setSite2(tempS2);
      tempPair.setEnergy(tempEnergy);      
      if(fabs(tempEnergy)>eCutOff)
	{
	  tempNbrParams +=1;
	  paramList.push_back(tempPair);
	}
      nbrOfParams=tempNbrParams;

    } 
  in.close();
}

void ParameterList::readParams_new(std::vector<std::string > subElements)
{
  offset_value=0;
  std::cout<<"starting to read params.."<<std::endl;

  int Mi=subElements.size();
  std::ifstream in(fileName.c_str());  
  if (!in)
    {
      std::cout<< "Cannot open file: "<<fileName<<std::endl;
      return;
    }

  Atom tempAtom = Atom();
  Pair tempPair= Pair();
  Triplet tempTriplet = Triplet();
  Quatuplet tempQ = Quatuplet();
  std::string tempS1;
  std::string tempS2;
  
  double tempDist;
  double tempEnergy;
  int tempNbrParams=0;
  
  int tuple_order;
  
  while(!(in.eof()))
    {
      in >> tuple_order;
      
      
      //offset, zeroth cluster
      if(tuple_order==0)
	{
	  in >> offset_value;
	}
      
      //singlets
      if(tuple_order==1)
	{
	  in >> tempEnergy;
	  tempAtom.setProperty(tempEnergy);	  
	  paramList_0.push_back(tempAtom); 	    
	}
      
      
      //pairs
      if(tuple_order==2)
	{	  
	  in>> tempDist;
	  in>> tempEnergy;
	  tempPair.setDistance(tempDist);
	  tempPair.setEnergy(tempEnergy);   
	  paramList.push_back(tempPair);
	}
      
      //triplets
      if(tuple_order==3)
	{

	  in >> tempDist;
	  tempTriplet.setDistance1(tempDist);

	  in >> tempDist;
	  tempTriplet.setDistance2(tempDist);
	  
	  in >> tempDist;
	  tempTriplet.setDistance3(tempDist);


	  in >> tempEnergy;
	  tempTriplet.setEnergy(tempEnergy);
	  
	  std::cout<<"=========+====+=++=+============++==+++=++"<<std::endl;
	  tempTriplet.printTriplet();
	  paramList_3.push_back(tempTriplet);	  
	}

      if(tuple_order==4)
	{
	  //need quatuplet code
	  
	  for(int i=0; i<6; i++)
	    {
	      in>>tempDist;
	      tempQ.setDistance(i,tempDist);
	    }
	  in >> tempEnergy;
	  std::cout<<"read param 4 energy "<<tempEnergy<<std::endl;
	  tempQ.setEnergy(tempEnergy);

	  paramList_4.push_back(tempQ);	  
	}
      

    } 

  std::cout<<"Done reading params. Number of pairs: "<< paramList.size()<<" Triplets: "<<paramList_3.size()<< " Quatuplets: "<< paramList_4.size()<<std::endl;
  in.close();
}


void ParameterList::unwrapSinglets(std::vector<Atom> singlet_list, std::vector<std::string> subelements)
{
  std::vector<Atom> new_list;
  
  Atom temp_atom = Atom();

  std::vector<std::vector<std::string> > unWrapped_elements;
  std::vector<std::vector<int> > clusterFunctions;
  double energy;
  for(int i=0; i<singlet_list.size(); i++)
    {
      std::vector<double> tomt;
      unWrapped_elements = symmetric_cluster_function(tomt,subelements);
      clusterFunctions = symmetric_cluster_function(tomt,subelements.size(), true);

      for(int j=0; j<unWrapped_elements.size(); j++)
	{
	  temp_atom.setType(unWrapped_elements[j][0]);
	  
	  int atom1;
	  for(int jj=0; jj<subelements.size(); jj++)
	    {
	      if(unWrapped_elements[j][0]==subelements[jj])
		{
		  atom1=jj;
		}
	    }
	  energy=0;
	  for(int k=0; k<clusterFunctions.size(); k++)
	    {
	      energy += clusterFunction(subelements.size(),atom1,clusterFunctions[k][0])*singlet_list[i].getProperty();
	    }

	  temp_atom.setProperty(energy);
	  new_list.push_back(temp_atom);
	}
    }
  paramList_0=new_list; //booyah
}

void ParameterList::unwrapPairs(std::vector<Pair> parList,std::vector<std::string> subelements)
{
  std::vector<Pair> new_list;
  Pair temp_pair = Pair();
  std::vector<std::vector<std::string> > unWrapped_elements;
  std::vector<std::vector<int> > clusterFunctions;
  std::vector<double> dists;
  dists.resize(1);
  double energy;
  for(int i=0; i<parList.size(); i++)
    {
      temp_pair.setDistance(parList[i].getDistance());
      dists[0]=parList[i].getDistance();
      unWrapped_elements = symmetric_cluster_function(dists,subelements);
      //true for reverse_sort
      clusterFunctions = symmetric_cluster_function(dists,subelements.size(),true);

      for(int j=0; j<unWrapped_elements.size(); j++)
	{

	  
	  int atom1;
	  int atom2;
	  for(int jj=0; jj<subelements.size(); jj++)
	    {
	      if(unWrapped_elements[j][0]==subelements[jj])
		{
		  atom1=jj;
		}
	      if(unWrapped_elements[j][1]==subelements[jj])
		{
		  atom2=jj;
		}

	    }
	  std::cout<<j<< " "<<unWrapped_elements[j][0]<< " "<<unWrapped_elements[j][1]<< " "<< atom1<< " "<<atom2<<std::endl;
	  energy=0;
	  for(int k=0; k<clusterFunctions.size(); k++)
	    {
	      energy += clusterFunction(subelements.size(),atom1,clusterFunctions[k][0])*
		clusterFunction(subelements.size(),atom2,clusterFunctions[k][1])*parList[i].getEnergy();
	      
	    }
	  temp_pair.setSite1(unWrapped_elements[j][0]);
	  temp_pair.setSite2(unWrapped_elements[j][1]);
	  temp_pair.setEnergy(energy);
	  new_list.push_back(temp_pair);
	}
    }
  paramList=new_list; //REALLY?!?!?! yes
}
	  
      

void ParameterList::unwrapTriplets(std::vector<Triplet> trip_list,std::vector<std::string> subelements)
{
  std::vector<Triplet> new_list;
  Triplet temp_trip = Triplet();
  std::vector<std::vector<std::string> > unWrapped_elements;
  std::vector<std::vector<int> > clusterFunctions;
  std::vector<double> dists;
  dists.resize(3);
  double energy;

  for(int i=0; i<trip_list.size(); i++)
    {
      temp_trip.setDistance1(trip_list[i].getDistance1());
      temp_trip.setDistance2(trip_list[i].getDistance2());
      temp_trip.setDistance3(trip_list[i].getDistance3());
      dists[0]=trip_list[i].getDistance1();
      dists[1]=trip_list[i].getDistance2();
      dists[2]=trip_list[i].getDistance3();
      
      unWrapped_elements = symmetric_cluster_function(dists,subelements);
      clusterFunctions = symmetric_cluster_function(dists,subelements.size(),true);
      
      for(int j=0; j<unWrapped_elements.size(); j++)
	{
	  int atom1;
	  int atom2;
	  int atom3;
	  for(int jj=0; jj<subelements.size(); jj++)
	    {
	      if(unWrapped_elements[j][0]==subelements[jj])
		{
		  atom1=jj;
		}
	      if(unWrapped_elements[j][1]==subelements[jj])
		{
		  atom2=jj;
		}

	      if(unWrapped_elements[j][2]==subelements[jj])
		{
		  atom3=jj;
		}
	    }
	  energy=0;

	  for(int k=0; k<clusterFunctions.size(); k++)
	    {
	      energy += clusterFunction(subelements.size(),atom1,clusterFunctions[k][0])*
		clusterFunction(subelements.size(),atom2,clusterFunctions[k][1])*
		clusterFunction(subelements.size(),atom3,clusterFunctions[k][2])*trip_list[i].getEnergy();
	    }
	  temp_trip.setSite1(unWrapped_elements[j][0]);
	  temp_trip.setSite2(unWrapped_elements[j][1]);
	  temp_trip.setSite3(unWrapped_elements[j][2]);
	  temp_trip.setEnergy(energy);
	  new_list.push_back(temp_trip);
	}
    }
  paramList_3=new_list;
}

void ParameterList::unwrapQuatuplets(std::vector<Quatuplet> quat_list, std::vector<std::string> subelements)
{
  std::vector<Quatuplet> new_list;
  Quatuplet temp_quat = Quatuplet();
    std::vector<std::vector<std::string> > unWrapped_elements;
  std::vector<std::vector<int> > clusterFunctions;
  std::vector<double> dists;
  dists.resize(6);
  double energy;
  
  for(int i=0; i<quat_list.size(); i++)
    {
      temp_quat.setDistance1(quat_list[i].getDistance1());
      temp_quat.setDistance2(quat_list[i].getDistance2());
      temp_quat.setDistance3(quat_list[i].getDistance3());
      temp_quat.setDistance4(quat_list[i].getDistance4());
      temp_quat.setDistance5(quat_list[i].getDistance5());
      temp_quat.setDistance6(quat_list[i].getDistance6());
      dists[0]=quat_list[i].getDistance1();
      dists[1]=quat_list[i].getDistance2();
      dists[2]=quat_list[i].getDistance3();
      dists[3]=quat_list[i].getDistance4();
      dists[4]=quat_list[i].getDistance5();
      dists[5]=quat_list[i].getDistance6();

      unWrapped_elements = symmetric_cluster_function(dists,subelements);
      clusterFunctions = symmetric_cluster_function(dists,subelements.size(),true);
      
      for(int j=0; j<unWrapped_elements.size(); j++)
	{
	  int atom1;
	  int atom2;
	  int atom3;
	  int atom4;
	  for(int jj=0; jj<subelements.size(); jj++)
	    {
	      if(unWrapped_elements[j][0]==subelements[jj])
		{
		  atom1=jj;
		}
	      if(unWrapped_elements[j][1]==subelements[jj])
		{
		  atom2=jj;
		}

	      if(unWrapped_elements[j][2]==subelements[jj])
		{
		  atom3=jj;
		}
	      if(unWrapped_elements[j][3]==subelements[jj])
		{
		  atom4=jj;
		}
	    }
	  
	  energy=0;

	  for(int k=0; k<clusterFunctions.size(); k++)
	    {
	      energy += clusterFunction(subelements.size(),atom1,clusterFunctions[k][0])*
		clusterFunction(subelements.size(),atom2,clusterFunctions[k][1])*
		clusterFunction(subelements.size(),atom3,clusterFunctions[k][2])*
		clusterFunction(subelements.size(),atom4,clusterFunctions[k][3])*quat_list[i].getEnergy();
	    }
	  temp_quat.setSite1(unWrapped_elements[j][0]);
	  temp_quat.setSite2(unWrapped_elements[j][1]);
	  temp_quat.setSite3(unWrapped_elements[j][2]);
	  temp_quat.setSite4(unWrapped_elements[j][3]);
	  std::cout<<"unwrap 4 energy "<<energy<<std::endl;
	  temp_quat.setEnergy(energy);
	  new_list.push_back(temp_quat);
	}
    }
  paramList_4=new_list;
}


  


Pair& ParameterList::getPair(int i)
{
  return paramList[i];
}

void ParameterList::printList()
{
  for (size_t i=0; i< paramList.size(); i++)
    {
      printPair(i);
    }
}

void ParameterList::printPair(int i)
{
  getPair(i).printPair();
}



Triplet& ParameterList::getTriplet(int i)
{
  return paramList_3[i];
}

int ParameterList::getNbrOfTriplets()
{
  return paramList_3.size();
}


Quatuplet& ParameterList::getQuatuplet(int i)
{
  return paramList_4[i];
}
int ParameterList::getNbrOfQuatuplets()
{
  return paramList_4.size();
}


Atom& ParameterList::getSinglet(int i)
{
  return paramList_0[i];
}

int ParameterList::getNbrOfSinglets()
{
  return paramList_0.size();
}

double ParameterList::getOffsetValue()
{
  return offset_value;
}

std::vector<Atom> ParameterList::returnSingletVector()
{
  return paramList_0;
}
