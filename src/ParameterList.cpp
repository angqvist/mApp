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

  readParamsATATStyle(subElements);
}


void ParameterList::readParamsWithPL(PairList pl)
{
  for(int i=0; i<pl.getNbrOfPairs(); i++)
    {
      paramList.push_back(pl.getPair(i));
      nbrOfParams++;
    }
}

int ParameterList::getNbrOfParams()
{
  if(nbrOfParams==-1)
    {
      std::cout<<" No params has been read, nbr of params = -1"<<std::endl;
    }
  return nbrOfParams;
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
  
  std::string tempS1;
  std::string tempS2;
  
  double tempDist;
  double tempEnergy;
  int tempNbrParams=0;
  
  int tuple_order;
  
  while(!(in.eof()))
    {
      //if( in.eof() ) break;
      in >> tuple_order;

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
      
      if(tuple_order==3)
	{
	  in >> tempEnergy;
	  tempTriplet.setEnergy(tempEnergy);
	  
	  in >> tempDist;
	  tempTriplet.setDistance1(tempDist);

	  in >> tempDist;
	  tempTriplet.setDistance2(tempDist);
	  
	  in >> tempDist;
	  tempTriplet.setDistance3(tempDist);
	  
	  paramList_3.push_back(tempTriplet);	  
	}

      if(tuple_order==4)
	{
	  //need quatuplet code

	  
	}
      

    } 
  in.close();


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


