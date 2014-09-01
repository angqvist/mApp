#include "MC_methods.hpp"
#include "LatticeList.hpp"
#include "Neighbour.hpp"
#include "NeighbourList.hpp"
#include <vector>
#include "MC.hpp"
#include <string>
double MC_totalEnergy(LatticeList ll, std::vector<NeighbourList> nl)
{
  double ret=0;
  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      ret += nl[i].getLocalEnergy(ll);
    }
  return ret/2.0;
}

double MC_orderValue(LatticeList ll, double orderLimit)
{
  int nbrGa=0;
  int nbrGe=0;
  int nbrNbrs=0;
  double ret=0;
  
  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      if(ll.getSite(i)!="Ga")
	{
	  continue;
	}
      nbrGa++;
      nbrNbrs=0;
      nbrGe=0;
      for(int j=0; j<ll.getNbrOfSites(); j++)
	{
	  if(j==i)
	    {
	      continue;
	    }
	  if(ll.getDistance(i,j)<orderLimit)
	    {
	      nbrNbrs++;
	      if(ll.getSite(j)=="Ge")
		{
		  nbrGe++;
		}
	    }
	}// end j loop
      ret += nbrGe/((double)nbrNbrs);
    }
  ret =1-ret/((double)nbrGa*(30.0/46.0)); //// NOTE SHOULD CHANGE IF GO AWAY FROM BaGaGe
  return ret;
}
	      


void MC_occupations(LatticeList ll, double &c6, double &k24, double &i16,std::string type)
{
  c6=0;
  k24=0;
  i16=0;

  for(int i=0; i<ll.getNbrOfSites(); i++)
    {

      if(ll.getSite(i)==type)
	{
	  if(i%46<16)
	    {
	      i16++;
	    }
	  if(i%46>15 && i%46<22)
	    {
	      c6++;
	    }
	  if(i%46>21)
	    {
	      k24++;
	    }
	}
    }

  int nbrOfCells=ll.getNbrOfSites()/46;
  c6 = c6/(6.0*nbrOfCells);
  i16 = i16/(16.0*nbrOfCells);
  k24 = k24/(24.0*nbrOfCells);
}


std::vector<int> neighbourCount(LatticeList ll, std::vector<std::string> subElements, double distMin,double distMax)
{
  std::vector<int> neighbourCount;
  neighbourCount.resize(ll.getNbrOfSites()*subElements.size());
  std::vector<int> tempCount;
  tempCount.resize(subElements.size());
  double tempDistance;
  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      for(int ii=0; ii<tempCount.size(); ii++)
      {
	tempCount[ii]=0;
	neighbourCount[i*subElements.size()+ii]=0;
      }
      for(int j=0; j<ll.getNbrOfSites(); j++)
	{
	  if(i==j)
	    {
	      continue;
	    }
	  
	  tempDistance=ll.getDistance(i,j);
	  
	  if(tempDistance>distMin && tempDistance<distMax)
	    {
	      for(int ii=0; ii<tempCount.size(); ii++)
		{
		  if(subElements[ii]==ll.getSite(j))
		    {
		      neighbourCount[i*subElements.size()+ii]++;
		    }
		}
	    }
	}//end j
    }//end i

  return  neighbourCount;
}
	  
// aborted since I might want to keep track of higher order neighbours, might need many intervvals

// void printPosAndNbrInfo(LatticeList ll, std::vector<std::string> subElements, double distMin,double distMax,std::string)
// {
//   std::vector<int> = neighbourCount(ll,subElements,distMin,distMax);

//   std::ofstream outF;
//   outF.open(filename.c_str());

//   outF<<"#Atom type X-Position Y-Position Z-position,
//   for(int i=0; i<ll.getNbrOfSites(); i++)
//     {
//       for(int j=0; j<subElements; j++)
