#include "MC.hpp"
#include "Neighbour.hpp"
#include "NeighbourList.hpp"
#include "LatticeList.hpp"
#include <math.h> 
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "MC_methods.hpp"

MC::MC()
{
  T=0;  
  ensemble="SGC";
  nbrOfSuccesfulSwaps=0;
  totalSwaps=0;
  kbeta=8.6173324e-5;
}

double MC::step(LatticeList &ll,std::vector<NeighbourList> nl)
{
  double energyDiff;

  energyDiff=stepSGC(ll,nl);
  //if(ensemble=="SGC")
  // {
  // }
  return energyDiff;
}


double MC::step(int N,LatticeList &ll,std::vector<NeighbourList> nl)
{
  double energyDiff=0; 
  for(int i=0; i<N; i++)
    {
      energyDiff += step(ll,nl);
    }   
  return energyDiff;
}

//double MC::averageStep(int mcsteps,int averageStep,LatticeList &ll,std::vector<std::vector<NeighbourList> > nlVectors,double &stdE,double &orderP,double &c6,double &k24, double &i16,double &stdc6,double &stdk24, double &stdi16 )
double MC::averageStep(int mcsteps,int averageStep,LatticeList &ll,std::vector<std::vector<NeighbourList> > nlVectors,std::vector<double> &data,std::vector<std::string> Wtype)
{
  std::string fileName = "dataFiles/mcConfs/binary/conf_new_bagage";
  for(int i=0; i<data.size(); i++)
    {
      data[i]=0;
    }

  //data:  E-0, std E-1,BG-2, stdBG-3,VOL-4,stdVOL-5,Lat-6,stdLat-7, c6-8, k24-9, i16-10, stdC6-11, stdk24-12, stdi16-13


  // std::vector<double> avg_prop;
  // std::vector<double> std_prop;
  int properties=nlVectors.size();

  // avg_prop.resize(properties);
  // std_prop.resize(properties);
  
  double averageNrg=0;
  double averageBG=0;
  double averageVOL=0;
  double averageLAT=0;

  double nrgSquared=0;
  double squareBG;
  double squareVOL;
  double squareLAT;

  double tempEnergy;

  double tempBG;
  double tempVOL;
  double tempLAT;
  double tempC6;
  double tempK24;
  double tempI16;
  double fastEnergy = MC_totalEnergy(ll,nlVectors[0]);
  
  std::vector<double> occupationVector;
  std::vector<std::string> wykSites;
  wykSites.push_back("6c");
  wykSites.push_back("24k");
  wykSites.push_back("16i");
  double tempVal;
  // i16 = 0;
  //  k24 = 0;
  //  c6  = 0;
  //  stdc6=0;
  //  stdk24 =0;
  //  stdi16 =0;
  //orderP=0;
  for(int k=0; k<averageStep; k++)
    {
      double energyDiff=0;  
      for(int i=0; i<mcsteps; i++)
	{
	  energyDiff += step(ll,nlVectors[0]);
	} 

      fastEnergy += energyDiff;      
      tempEnergy = fastEnergy;

      averageNrg += tempEnergy;
      nrgSquared += pow(tempEnergy,2.0);


      
      for(int i=0; i<Wtype.size(); i++)
	{
	  occupationVector=ll.getWyckoffOccupancy(wykSites,Wtype[i]);
	  tempC6=occupationVector[0];
	  tempK24=occupationVector[1];
	  tempI16=occupationVector[2];

	  
	  //MC_occupations(ll,tempC6,tempK24,tempI16,Wtype[i]);
	  data[2*properties+i*5+i] +=tempC6;
	  data[2*properties+1+i*5+i] +=pow(tempC6,2.0);

	  data[2*properties+2+i*5+i] +=tempK24;
	  data[2*properties+3+i*5+i] +=pow(tempK24,2.0);

	  data[2*properties+4+i*5+i] +=tempI16;
	  data[2*properties+5+i*5+i] +=pow(tempI16,2.0);
	}
 
      

      
      for( int i=0; i< properties-1; i++)
	{
	  tempVal=MC_totalEnergy(ll,nlVectors[i+1]);
	  data[2*(i+1)+0] += tempVal;
	  data[2*(i+1)+1] += pow(tempVal,2.0);
	}
            
      

    }
  
  
  

  data[0]=averageNrg/(double)averageStep;
  data[1] = sqrt(nrgSquared/(double)averageStep-pow(data[0],2.0));
  for(int i=1; i<data.size()/2; i++)
    {
      data[2*i]=data[2*i]/(double)averageStep;
      data[2*i+1]=sqrt(data[2*i+1]/(double)averageStep-pow(data[2*i],2.0));
    }
  




  // data[10] =data[10]/(double)averageStep;
  // data[9] = data[9]/(double)averageStep;
  // data[8]  = data[8]/(double)averageStep;
           
  // data[12] = sqrt(data[12]/(double)averageStep-pow(data[9],2.0));
  // data[13] = sqrt(data[13]/(double)averageStep-pow(data[10],2.0));
  // data[11] = sqrt(data[11]/(double)averageStep-pow(data[8],2.0));
  

  // data[2]=data[2]/(double)averageStep;
  // data[4]=data[4]/(double)averageStep;
  // data[6]=data[6]/(double)averageStep;

  // data[3]=sqrt(data[3]/(double)averageStep-pow(data[2],2.0));
  // data[5]=sqrt(data[5]/(double)averageStep-pow(data[4],2.0));
  // data[7]=sqrt(data[7]/(double)averageStep-pow(data[6],2.0));

  // orderP = orderP/(double)averageStep;
  return data[0];
}



double MC::stepSGC(LatticeList &ll,std::vector<NeighbourList> nl)
{
  totalSwaps++;
  double energyDiff;
  double nrgyBefore;
  double nrgyAfter;
  // double energyBefore1;
  //  double energyBefore2;
  double energyAfter1;
  double energyAfter2;


  int k=rand()%(ll.get_original_atoms_count());
  int l=rand()%(ll.get_original_atoms_count());
  while(ll.getSite(l) == ll.getSite(k))
    {
      l=rand()%(ll.get_original_atoms_count());
    }
  // energyBefore1=nl[k].getCurrentLocalEnergy();
  // energyBefore2=nl[l].getCurrentLocalEnergy();

  nrgyBefore=nl[k].getLocalEnergyForMC(ll)+nl[l].getLocalEnergyForMC(ll);
  //nrgyBefore= nl[k].getCurrentLocalEnergy() + nl[l].getCurrentLocalEnergy();
  std::string tempSite1= ll.getSite(k);
  std::string tempSite2= ll.getSite(l);
  ll.setSite(k,tempSite2);
  ll.setSite(l,tempSite1);
  // energyAfter1=nl[k].getLocalEnergy(ll);
  // energyAfter2=nl[l].getLocalEnergy(ll);
  // nrgyAfter=energyAfter1 + energyAfter2;

  nrgyAfter=nl[k].getLocalEnergyForMC(ll)+nl[l].getLocalEnergyForMC(ll);
  energyDiff=nrgyAfter-nrgyBefore;
  //std::cout<<energyDiff<<std::endl;

  if(energyDiff>0)
    {
      if(exp(-energyDiff/(kbeta*T))<(rand()/(double)RAND_MAX)) //revert
	{
	  ll.setSite(k,tempSite1);
	  ll.setSite(l,tempSite2);
	  return 0;
	}
    }
  nbrOfSuccesfulSwaps++;
  return energyDiff;  
}

void MC::setTemperature(double newTemp)
{
  T=newTemp;
}

void MC::setEnsemble(std::string newEnsemble)
{
  ensemble=newEnsemble;
}

void MC::setKBeta(double newKB)
{
  kbeta=newKB;
}

void MC::resetSwapCounts()
{
  nbrOfSuccesfulSwaps=0;
  totalSwaps=0;
}

double MC::getSwapRatio()
{
  return (double) nbrOfSuccesfulSwaps/totalSwaps;
}


void MC::printInfo(std::string fileName,LatticeList ll, double minDist, double maxDist)
{
  std::vector<std::string > subE; // NEED TO FIX THIS!!
  subE.push_back("Al");
  subE.push_back("Ga");
  subE.push_back("Si");

  //std::vector<int> nc = neighbourCount(ll,subE,0.5,3.0);

  
  std::ofstream outF;
  outF.open(fileName.c_str());
  outF<<ll.getNbrOfSites()<<std::endl;
  outF<<"#Temperature: "<<T<< " Atom type X-Position Y-Position Z-position "<<std::endl;
  // for(int i=0; i<subE.size(); i++)
  //   {
  //     outF<<subE[i]<< " ";
  //   }
  // outF<<std::endl;

  std::string type;
  double xpos;
  double ypos;
  double zpos;
  


  outF <<"    "<<ll.getLx()<<"    0.0000000000000000   0.0000000000000000"<<std::endl;
  outF<<"    0.0000000000000000    "<<ll.getLy()<<"   0.0000000000000000"<<std::endl;
  outF<<"    0.0000000000000000    0.0000000000000000   "<<ll.getLz()<<std::endl;

 

   for(int i=0; i<ll.get_original_atoms_count(); i++)
     {
       type=ll.getAtomInfo(i,xpos,ypos,zpos);
       outF<<type<< " "<<xpos<< " "<<ypos<< " "<<zpos<<" ";
       // for(int j=0; j<subE.size(); j++)
       // 	 {
       // 	   outF<<nc[i*subE.size()+j]<< " ";
       // 	 }
       outF<<std::endl;
     }
}
