#include <iostream>
#include <string>
#include <fstream>
#include "LatticeList.hpp"
#include <stdlib.h>
#include <cmath>
#include <string>
using namespace std;

// (A bit inelegant, nbrOfAtoms isnt really nbr of atoms but the number of atoms in the unit cell times 3)-- I think not
// i.e the number of numbers required to get the positions of all the atoms
// (Number of sites, is the number of total sites in the lattice times three) -- I think not

// default constructor just makes a unit cell lattice
LatticeList::LatticeList() 
{
  cellSizeX=1;
  cellSizeY=1;
  cellSizeZ=1;
  nbrOfAtoms=46; //46*3
  nbrOfSites=nbrOfAtoms*cellSizeX*cellSizeY*cellSizeZ;
  posList2.resize(nbrOfSites*3);
  atomTypeList.resize(nbrOfSites);
  latticeConstant=10.9868239608;
  Lx=cellSizeX*latticeConstant;
  Ly=cellSizeY*latticeConstant;
  Lz=cellSizeZ*latticeConstant;
  fileName="configs/confs4/config_0";
  readIdealPos2();

}

LatticeList::LatticeList(int sizeX, int sizeY, int sizeZ)
{
  cellSizeX=sizeX;
  cellSizeY=sizeY;
  cellSizeZ=sizeZ;
  nbrOfAtoms=46; //46
  nbrOfSites=nbrOfAtoms*cellSizeX*cellSizeY*cellSizeZ;
  posList2.resize((nbrOfSites*3));
  atomTypeList.resize(nbrOfSites);
  latticeConstant=10.9868239608;
  Lx=cellSizeX*latticeConstant;
  Ly=cellSizeY*latticeConstant;
  Lz=cellSizeZ*latticeConstant;
  fileName="configs/confs4/config_0";
  readIdealPos2();
}
LatticeList::LatticeList(int sizeX, int sizeY, int sizeZ, int newNbrOfAtoms, string newFileName) 
{
  cellSizeX=sizeX;
  cellSizeY=sizeY;
  cellSizeZ=sizeZ;
  nbrOfAtoms=newNbrOfAtoms;
  nbrOfSites=nbrOfAtoms*cellSizeX*cellSizeY*cellSizeZ;
  posList2.resize((nbrOfSites*3));
  atomTypeList.resize(nbrOfSites);
  latticeConstant=10.9868239608;  
  Lx=cellSizeX*latticeConstant;
  Ly=cellSizeY*latticeConstant;
  Lz=cellSizeZ*latticeConstant;

  fileName=newFileName;
  readIdealPos2();
}
LatticeList::LatticeList(int sizeX, int sizeY, int sizeZ, int newNbrOfAtoms, int nbrOfProperties, string newFileName,std::vector<std::string> subElements) 
{
  cellSizeX=sizeX;
  cellSizeY=sizeY;
  cellSizeZ=sizeZ;
  nbrOfAtoms=newNbrOfAtoms;
  nbrOfSites=nbrOfAtoms*cellSizeX*cellSizeY*cellSizeZ;
  posList2.resize((nbrOfSites*3));
  atomTypeList.resize(nbrOfSites);
  latticeConstant=10.9868239608;  
  Lx=cellSizeX*latticeConstant;
  Ly=cellSizeY*latticeConstant;
  Lz=cellSizeZ*latticeConstant;
  properties.resize(nbrOfProperties);
  fileName=newFileName;
  readIdealPos3(subElements);

}




void LatticeList::readIdealPos()
{
  ifstream in(fileName.c_str());  
  if (!in)
    {
      cout<< "Cannot open file.\n";
      return;
    }
  for(int i=0; i<(nbrOfAtoms*3); i++)
    {
      in >>posList2[i];
    }
  for(int i=(nbrOfAtoms*3); i<cellSizeX*(nbrOfAtoms*3); i+=3)
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)]+latticeConstant;
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)];
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)];
    }
  for(int i=(nbrOfAtoms*3)*cellSizeX; i<cellSizeX*cellSizeY*(nbrOfAtoms*3); i+=3)
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)*cellSizeX];
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)*cellSizeX]+latticeConstant;
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)*cellSizeX];
    }
  for(int i=(nbrOfAtoms*3)*cellSizeX*cellSizeY; i<cellSizeX*cellSizeY*cellSizeZ*(nbrOfAtoms*3); i+=3 )
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)*cellSizeX*cellSizeY];
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)*cellSizeX*cellSizeY];
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)*cellSizeX*cellSizeY]+latticeConstant;
    }
  in.close(); 
}

// reads the new type of configs -keeping the old one, don't burn bridges.
void LatticeList::readIdealPos2()
{
  ifstream in(fileName.c_str());  
  if (!in)
    {
      cout<< "Cannot open file.\n";
      return;
    }
  double tempLatticeConstant;
  for(int i=0; i<3; i++)
    {
      for( int j=0; j<3; j++)
	{
	  in >> tempLatticeConstant;
	  if( j != i)
	    {
	      continue;
	    }
	  if(i==0)
	    {
	      latticeConstant=tempLatticeConstant; // only got one latticeconstant
	    }
	}
    }
  Lx=cellSizeX*latticeConstant;
  Ly=cellSizeY*latticeConstant;
  Lz=cellSizeZ*latticeConstant;
  
  
  double tempEnergy;
  in >> tempEnergy;
  setEnergy(tempEnergy);
  in >> tempEnergy; //bandgap now
  setBandGap(tempEnergy);
  in >> tempEnergy; // volume
  setVolume(tempEnergy);
  in >> tempEnergy; // averageLatticeConstant
  setAverageLatticeConstant(tempEnergy);
  std::string tempSite;
  double tempPos1;
  double tempPos2;
  double tempPos3;
  int loopIndex=0;
  while(!(in.eof()))
    {
      if(loopIndex+1>nbrOfAtoms)
	{
	  break;
	}
      in >>tempSite;
      in >>tempPos1;
      in >>tempPos2;
      in >>tempPos3;
      bool newAtom=true;
      if(elements.size()==0)
	{
	  elements.push_back(tempSite);
	  elementCounts.push_back(1);
	}
      else
	{
	  for(int eleIndex=0; eleIndex<elements.size(); eleIndex++)
	    {
	      if(tempSite==elements[eleIndex])
		{
		  elementCounts[eleIndex]++;
		  newAtom=false;
		  break;
		}
	      
	    }//end for
	  if(newAtom)
	    {
	      elementCounts.push_back(1);
	      elements.push_back(tempSite);
	    }
	}// end else

      if(tempSite=="Ba")
	{
	  continue;
	}

      atomTypeList[loopIndex]=tempSite;
      posList2[3*loopIndex]=tempPos1;
      posList2[3*loopIndex+1]=tempPos2;
      posList2[3*loopIndex+2]=tempPos3;
      loopIndex++;      
    }
  
  for(int i=(nbrOfAtoms*3); i<cellSizeX*(nbrOfAtoms*3); i+=3)
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)]+latticeConstant;
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)];
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)];
    }
  for(int i=(nbrOfAtoms*3)*cellSizeX; i<cellSizeX*cellSizeY*(nbrOfAtoms*3); i+=3)
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)*cellSizeX];
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)*cellSizeX]+latticeConstant;
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)*cellSizeX];
    }
  for(int i=(nbrOfAtoms*3)*cellSizeX*cellSizeY; i<cellSizeX*cellSizeY*cellSizeZ*(nbrOfAtoms*3); i+=3 )
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)*cellSizeX*cellSizeY];
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)*cellSizeX*cellSizeY];
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)*cellSizeX*cellSizeY]+latticeConstant;
    }
  in.close(); 
}




void LatticeList::readIdealPos3(std::vector<std::string> subElements)
{
  ifstream in(fileName.c_str());  
  if (!in)
    {
      cout<< "Cannot open file.\n";
      return;
    }
  double tempLatticeConstant;
  for(int i=0; i<3; i++)
    {
      for( int j=0; j<3; j++)
	{
	  in >> tempLatticeConstant;
	  if( j != i)
	    {
	      continue;
	    }
	  if(i==0)
	    {
	      latticeConstant=tempLatticeConstant; // only got one latticeconstant
	    }
	}
    }
  Lx=cellSizeX*latticeConstant;
  Ly=cellSizeY*latticeConstant;
  Lz=cellSizeZ*latticeConstant;
  
  for(int i=0; i<properties.size(); i++)
    {
      in >> properties[i];
    }

  std::string tempSite;
  double tempPos1;
  double tempPos2;
  double tempPos3;
  int loopIndex=0;
  while(!(in.eof()))
    {
      if(loopIndex+1>nbrOfAtoms)
	{
	  break;
	}
      in >>tempSite;
      in >>tempPos1;
      in >>tempPos2;
      in >>tempPos3;
      bool newAtom=true;


  
      
      if(elements.size()==0)
	{
	  elements.push_back(tempSite);
	  elementCounts.push_back(1);
	}
      else
	{
	  for(int eleIndex=0; eleIndex<elements.size(); eleIndex++)
	    {
	      if(tempSite==elements[eleIndex])
		{
		  elementCounts[eleIndex]++;
		  newAtom=false;
		  break;
		}
	      
	    }//end for
	  if(newAtom)
	    {
	      elementCounts.push_back(1);
	      elements.push_back(tempSite);
	    }
	}// end else



      bool illegalAtom=true;
      for(int i=0; i<subElements.size(); i++)
      	{
      	  if(subElements[i]==tempSite)
      	    {
      	      illegalAtom=false;
      	    }
      	}
      if(illegalAtom)
      	{
      	  continue;
      	}

 

      atomTypeList[loopIndex]=tempSite;
      posList2[3*loopIndex]=tempPos1;
      posList2[3*loopIndex+1]=tempPos2;
      posList2[3*loopIndex+2]=tempPos3;
      loopIndex++;      
    }
  
  for(int i=(nbrOfAtoms*3); i<cellSizeX*(nbrOfAtoms*3); i+=3)
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)]+latticeConstant;
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)];
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)];
    }
  for(int i=(nbrOfAtoms*3)*cellSizeX; i<cellSizeX*cellSizeY*(nbrOfAtoms*3); i+=3)
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)*cellSizeX];
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)*cellSizeX]+latticeConstant;
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)*cellSizeX];
    }
  for(int i=(nbrOfAtoms*3)*cellSizeX*cellSizeY; i<cellSizeX*cellSizeY*cellSizeZ*(nbrOfAtoms*3); i+=3 )
    {
      posList2[i]=posList2[i-(nbrOfAtoms*3)*cellSizeX*cellSizeY];
      posList2[i+1]=posList2[i+1-(nbrOfAtoms*3)*cellSizeX*cellSizeY];
      posList2[i+2]=posList2[i+2-(nbrOfAtoms*3)*cellSizeX*cellSizeY]+latticeConstant;
    }
  in.close(); 
}



void LatticeList::printList()
{
  cout<<"PasCarBen01.cif"<<endl;
  cout<<"   "<<"10.9868239607846370"<<endl;
  cout<<"    1.0000000000000000    0.0000000000000000   0.0000000000000000"<<endl;
  cout<<"    0.0000000000000000    1.0000000000000000   0.0000000000000000"<<endl;
  cout<<"    0.0000000000000000    0.0000000000000000   1.0000000000000000"<<endl;
  
  //cout<<"   Ba   Ga   Ge" <<endl;
  if(elements.size()==1)
    {
      cout<<"   "<<elements[0]<<endl;
      cout<<"     "<<elementCounts[0]<<endl;
    }
  else if(elements.size()==2)
    {
      cout<<"   "<<elements[0]<<"   "<<elements[1] <<endl;
      cout<<"     "<<elementCounts[0]<<"    "<<elementCounts[1]<<endl;
    }
  else if(elements.size()==3)
    {
      cout<<"   "<<elements[0]<<"   "<<elements[1] <<"   "<<elements[2]<<endl;
      cout<<"     "<<elementCounts[0]<<"    "<<elementCounts[1]<<"   "<<elementCounts[2]<<endl;

    }
  else if(elements.size()==4)
    {
      cout<<"   "<<elements[0]<<"   "<<elements[1] <<"   "<<elements[2]<<"   "<<elements[3]<<endl;
      cout<<"     "<<elementCounts[0]<<"    "<<elementCounts[1]<<"   "<<elementCounts[2]<<"   "<<elementCounts[3]<<endl;
    }
  //cout<<"     8    16    30"<<endl;
  cout<<"Direct"<<endl;
  cout<<"  0.0000000000000000  0.0000000000000000  0.0000000000000000"<<endl;
  cout<<"  0.5000000000000000  0.5000000000000000  0.5000000000000000"<<endl;
  cout<<"  0.2500000000000000  0.5000000000000000  0.0000000000000000"<<endl;
  cout<<"  0.7500000000000000  0.5000000000000000  0.0000000000000000"<<endl;
  cout<<"  0.5000000000000000  0.0000000000000000  0.2500000000000000"<<endl;
  cout<<"  0.0000000000000000  0.2500000000000000  0.5000000000000000"<<endl;
  cout<<"  0.5000000000000000  0.0000000000000000  0.7500000000000000"<<endl;
  cout<<"  0.0000000000000000  0.7500000000000000  0.5000000000000000"<<endl;


  for(int eleIndex=0; eleIndex<elements.size(); eleIndex++)
    {
  
      for(int i=0; i<nbrOfSites; i++)
	{
	  if(atomTypeList[i]==elements[eleIndex])
	    {
	      cout<<"  ";
	      for(int j=0; j<3;j++)
		{
		  cout<<posList2[i*3+j]/10.9868239607846370<< "  ";
		}
	      cout<<endl;
	    }      
	}
    }
}


void LatticeList::printList(std::string filename)
{

  std::ofstream outF;
  outF.open(filename.c_str());
  outF<<"PasCarBen01.cif"<<endl;
  outF<<"   "<<"10.9868239607846370"<<endl;
  outF<<"    1.0000000000000000    0.0000000000000000   0.0000000000000000"<<endl;
  outF<<"    0.0000000000000000    1.0000000000000000   0.0000000000000000"<<endl;
  outF<<"    0.0000000000000000    0.0000000000000000   1.0000000000000000"<<endl;
  outF<<"   Ba   Ga   Ge" <<endl;
  outF<<"     8    16    30"<<endl;
  outF<<"Direct"<<endl;
  outF<<"  0.0000000000000000  0.0000000000000000  0.0000000000000000"<<endl;
  outF<<"  0.5000000000000000  0.5000000000000000  0.5000000000000000"<<endl;
  outF<<"  0.2500000000000000  0.5000000000000000  0.0000000000000000"<<endl;
  outF<<"  0.7500000000000000  0.5000000000000000  0.0000000000000000"<<endl;
  outF<<"  0.5000000000000000  0.0000000000000000  0.2500000000000000"<<endl;
  outF<<"  0.0000000000000000  0.2500000000000000  0.5000000000000000"<<endl;
  outF<<"  0.5000000000000000  0.0000000000000000  0.7500000000000000"<<endl;
  outF<<"  0.0000000000000000  0.7500000000000000  0.5000000000000000"<<endl;


  for(int i=0; i<nbrOfSites; i++)
    {
      if(atomTypeList[i]=="Ga")
	{
	  outF<<"  ";
	  for(int j=0; j<3;j++)
	    {
	      outF<<posList2[i*3+j]/10.9868239607846370<< "  ";
	    }
	  outF<<endl;
	}      
    }
  for(int i=0; i<nbrOfSites; i++)
    {
      if(atomTypeList[i]=="Ge")
	{
	  outF<<"  ";
	  for(int j=0; j<3;j++)
	    {
	      outF<<posList2[i*3+j]/10.9868239607846370<< "  ";
	    }
	  outF<<endl;
	}      
    }




}





std::string LatticeList::getAtomInfo(int i,double &xpos, double &ypos, double &zpos)
{
  xpos=posList2[i*3];
  ypos=posList2[i*3+1];
  zpos=posList2[i*3+2];
  return getSite(i);
}






double LatticeList::getDistance(int i, int j)
{
  double dx=fabs(posList2[i*3]-posList2[j*3]);
  double dy=fabs(posList2[i*3+1]-posList2[j*3+1]);
  double dz=fabs(posList2[i*3+2]-posList2[j*3+2]);
  // std::cout<<dx<<" "<< Lx-dx<<" "<< std::min(dx,Lx-dx) <<std::endl;
  dx=std::min(dx,Lx-dx);
  dy=std::min(dy,Ly-dy);
  dz=std::min(dz,Lz-dz);
  //if(dx<0 ||  dy<0 || dz<0)
    //  {
  //std::cout<<"dx= : "<<pow(dx,2.0)<< " "<<"dy= : "<<dy<<" dz = :"<<dz<<std::endl;
  //  }
  return sqrt(pow(dx,2.0)+pow(dy,2.0)+pow(dz,2.0));
}

int LatticeList::getNbrOfSites()
{
  return nbrOfSites;
}

void LatticeList::setRandomSites(int nbrA, std::string type1, std::string type2)
{
  int k;
  int countA=0;
  for(size_t i=0; i<nbrOfSites; i++)
    {
      atomTypeList[i]=type2;
    }

  while(nbrA*cellSizeX*cellSizeY*cellSizeZ>countA)
    {
      k = rand()%nbrOfSites;
      if(atomTypeList[k]==type2)
	{
	  atomTypeList[k]=type1;
	  countA++;
	}
    }
}

void LatticeList::setRandomSites(int nbrA,int nbrB, std::string type1, std::string type2,std::string type3)
{
  int k;
  int countA=0;
  int countB=0;
  for(size_t i=0; i<nbrOfSites; i++)
    {
      atomTypeList[i]=type3;
    }
  while(nbrA*cellSizeX*cellSizeY*cellSizeZ>countA)
    {
      k = rand()%nbrOfSites;
      if(atomTypeList[k]==type3)
	{
	  atomTypeList[k]=type1;
	  countA++;
	}
    }
  while(nbrB*cellSizeX*cellSizeY*cellSizeZ>countB)
    {
      k = rand()%nbrOfSites;
      if(atomTypeList[k]==type3)
	{
	  atomTypeList[k]=type2;
	  countB++;
	}
    }
}


void LatticeList::setRandomSites(int nbrA,int nbrB, int nbrC, std::string type1, std::string type2,std::string type3,std::string type4)
{
  int k;
  int countA=0;
  int countB=0;
  int countC=0;
  for(size_t i=0; i<nbrOfSites; i++)
    {
      atomTypeList[i]=type4;
    }
  while(nbrA*cellSizeX*cellSizeY*cellSizeZ>countA)
    {
      k = rand()%nbrOfSites;
      if(atomTypeList[k]==type4)
	{
	  atomTypeList[k]=type1;
	  countA++;
	}
    }
  while(nbrB*cellSizeX*cellSizeY*cellSizeZ>countB)
    {
      k = rand()%nbrOfSites;
      if(atomTypeList[k]==type4)
	{
	  atomTypeList[k]=type2;
	  countB++;
	}
    }
  while(nbrC*cellSizeX*cellSizeY*cellSizeZ>countC)
    {
      k = rand()%nbrOfSites;
      if(atomTypeList[k]==type4)
	{
	  atomTypeList[k]=type3;
	  countC++;
	}
    }
}







std::string LatticeList::getSite(int i)
{
  return atomTypeList[i];
}


void LatticeList::setSite(int index,std::string newSite)
{
  atomTypeList[index]=newSite;
}

void LatticeList::setEnergy(double newEnergy)
{
  energy=newEnergy;
}


double LatticeList::getEnergy()
{
  return energy;
}


void LatticeList::setBandGap(double newBandGap)
{
  bandGap=newBandGap;
}

double LatticeList::getBandGap()
{
  return bandGap;
}


void LatticeList::setVolume(double newVolume)
{
  volume=newVolume;
}


double LatticeList::getVolume()
{
  return volume;
}

void LatticeList::setAverageLatticeConstant(double newAverageLatticeConstant)
{
  averageLatticeConstant = newAverageLatticeConstant;
}



double LatticeList::getAverageLatticeConstant()
{
  return averageLatticeConstant;
}


double LatticeList::getLx()
{
  return Lx;
}
double LatticeList::getLy()
{
  return Ly;
}
double LatticeList::getLz()
{
  return Lz;
}

int LatticeList::getNumberOfProperties()
{
  return properties.size();
}


double LatticeList::getProperty(int index)
{
  if(index >= properties.size())
    {
      std::cout<<"ERROR: trying to fetch property from lattice beyond the limit of property vector"<<std::endl;
      std::cout<<"ERROR: trying to get property "<<index<<" length of property: "<<properties.size()<<std::endl;
      return 0;
    }
  return properties[index];
}


double LatticeList::getConcentration(std::string type)
{
  for(int i=0; i<elements.size(); i++)
    {
      if(type==elements[i])
	{
	  return (double)elementCounts[i]/(double)nbrOfAtoms;
	}
    }
  return 0;
}
