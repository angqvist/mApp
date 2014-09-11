#include "helperFunctions.hpp"
#include <iostream>
#include <fstream>
#include "LatticeList.hpp"
#include "PairList.hpp"
#include "Pair.hpp"
#include "TripletList.hpp"
#include "Triplet.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <gsl/gsl_multimin.h>
#include "NeighbourList.hpp"
#include "ParameterList.hpp"
#include "Neighbour.hpp"
#include <math.h> 
#include <algorithm> 
#include <sstream>
#include<omp.h>
  struct minPar
  {
    gsl_matrix * parA;
    gsl_matrix * parAtA;
    gsl_vector * parEnergy;
    gsl_vector * parftA;
    gsl_vector * parD;
    gsl_vector * parB;
    double parMu;
    double parLambda;
    double parAlpha; // maybe not need alpha
    int parRows;
    int parColumns;
  };
//calc norm of error for valid and training as function of how many configs in training - do this for maybe three different cutoffs
// calc error for valid&training for all cutoffs
// calc number of nonzero parameters as function of cutoff

void printData(int configStart, int configStop,int configStep, int nbrOfConfigs, std::string posFile, std::string energyFile, double paramLimit, std::string nrgyFiles,std::string paramDataFile,bool doNrgNorm,int nbrFitParams,bool doTriplet,bool doParams,double doubleCO,double tripletCO,bool doBandGap,bool doVolume,bool doLattice,bool printTrainingAndValidationEnergy,std::string energyTrainFile,std::string energyValidFile,std::string confFileName,std::string confDirectory,std::string ECIParamOutFile,std::vector<std::string> subElements,bool printECIParams,bool doATAT,int printParamAtIndex,bool doCV)
{
  //initial stuff
  std::cout<<"start print data"<<std::endl;
  std::vector<LatticeList> dftPos = readConfig(confDirectory,nbrOfConfigs,46,4,subElements);
  std::vector<double> energies;
  for(int i=0; i<nbrOfConfigs; i++)
    {
       std::cout<<i<< " "<<dftPos[i].getBandGap()<< " "<<dftPos[i].getVolume()<< " "<<dftPos[i].getAverageLatticeConstant()<< " "<<dftPos[i].getEnergy()<<std::endl;
      
      if(doBandGap)
	{
	  energies.push_back(dftPos[i].getBandGap());
	}
      else if(doVolume)
	{
	  energies.push_back(dftPos[i].getVolume());
	}
      else if(doLattice)
	{
	  energies.push_back(dftPos[i].getAverageLatticeConstant());
	}      
      else
	{
	  energies.push_back(dftPos[i].getEnergy());
	}
    }
  std::cout<<"read data: "<< energies.size()<< " entries"<<std::endl;


  //shuffleLists(energies,dftPos);

  
  LatticeList lista = LatticeList(1,1,1);
  if(subElements.size()==2)
    {
  lista.setRandomSites(16,subElements[0],subElements[1]);
    }
  else if(subElements.size()==3)
    {
      lista.setRandomSites(16,10,subElements[0],subElements[1],subElements[2]);
    }
  else if(subElements.size()==4)
    {
      lista.setRandomSites(16,10,10,subElements[0],subElements[1],subElements[2],subElements[3]);
    }
  else
    {
      std::cout<<"Error in subelements, size is: "<<subElements.size()<<std::endl;
    }
  PairList pl = PairList();
  std::vector<std::string> subelements;
  std::vector<double> nrgyValid;
  std::vector<double> nrgyTrain;
  nrgyValid.resize(round((configStop-configStart)/(double)configStep));
  nrgyTrain.resize(round((configStop-configStart)/(double)configStep));

 
  pl.initializePairs(lista,subElements,doubleCO);  
  TripletList tl = TripletList();

  if(doTriplet)
    {
      tl.initializeTriplets(lista,subElements,tripletCO);
      std::cout<<"nbr of triplets: "<<tl.getNbrOfTriplets()<<std::endl;
    }


  std::vector<double> X;
  std::vector<double> X2;
  std::vector<double> params;
  std::vector<double> calcEnergies;
  const double alpha=0.1;
  const double mu=0.65;
  //const double mu=0.65;
  const double lambda=100;
  const bool doSB=true;
  const int sbIters=1000;
  const double sbTol=1e-6;
  const int bfgsIters=5000;
  const double bfgsTol=1e-5;
  const bool verbal=false;
  double tempNrgy;
  double tempParamNorm;
  std::vector<double> cutOffVector;
  std::vector<double> distances; //for ATAT
  // ===================== del me
  // X2=getAwithATAT(dftPos,nbrOfConfigs,subElements,doubleCO,distances,false);

  // int numbConfs=70;
  // X = getAMatrixFromExistingOne(X2,numbConfs,distances.size());

  // //  X=getAwithATAT(dftPos,numbConfs,subElements,doubleCO,distances,false);
  // double tempPNorm=0;
  // double tempErrorT=0;
  // double tempErrorV=0;
  // double nonZero=0;

  // for(double mu2=0.000001; mu2<1e14; mu2 *=1.05)
  //   {
  //     nonZero=0;
  //     const double mu3=mu2;
  //     tempPNorm=0;
  //     tempErrorT=0;
  //     tempErrorV=0;
  //     params = doMinimize(X,distances.size(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);
  //     std::vector<double> tempEnergyAll = energyFromParams(params,X2);
  //     for(int jj2=0; jj2<numbConfs; jj2++)
  // 	{
  // 	  tempErrorT +=pow(((tempEnergyAll[jj2]-energies[jj2])),2.0);
  // 	}

  //     for(int jj2=numbConfs; jj2<nbrOfConfigs; jj2++)
  // 	{
  // 	  tempErrorV +=pow(((tempEnergyAll[jj2]-energies[jj2])),2.0);
  // 	}
  //     double medel=0;
  //     double std=0;

  //     for(int jj2=4; jj2<params.size(); jj2++)
  // 	{
  // 	  medel += fabs(params[jj2]);
  // 	  std += pow(params[jj2],2.0);
  // 	}
  //     medel = medel/((double)params.size()-4);
  //     std=sqrt(std/((double)params.size()-4)-pow(medel,2.0));

  //     // std::cout<<medel<< " "<<std<<" "<<medel-0.5*medel<<std::endl;
  //     for(int jj2=0; jj2<params.size(); jj2++)
  // 	{
	 
  // 	  //  std::cout<<params[jj2]<<std::endl;
	    
  // 	  tempPNorm += fabs(params[jj2]);
  // 	  // if(fabs(params[jj2])>1e-4)
  // 	  //   {
  // 	  //     nonZero++;
  // 	  //   }
	  
  // 	  if(fabs(params[jj2])>(medel-0.5*medel))
  // 	    {
  // 	      nonZero++;
  // 	    }
  // 	}
  //     std::cout<<mu2<< " "<<std::sqrt(tempErrorT/numbConfs)<< " "<<std::sqrt(tempErrorV/((double)nbrOfConfigs-numbConfs))<<" "<<tempPNorm<<" "<<nonZero<<std::endl;
  //  }


      //===================del me

  cutOffVector.push_back(3.0);
  cutOffVector.push_back(4.0);
  cutOffVector.push_back(5.0);
  cutOffVector.push_back(6.0);
  cutOffVector.push_back(7.0);
  cutOffVector.push_back(8.0);
  cutOffVector.push_back(9.0);
  cutOffVector.push_back(10.0);
  if(doNrgNorm)
    {

      std::cout<<"start doing energy norm"<<std::endl;
  //     for( int CO=0;CO<cutOffVector.size();CO++)
  // 	{ 
  if(doTriplet)
    {
      std::cout<<"generating X2 matrix with triplets. Number of triplets: "<<tl.getNbrOfTriplets()<< ". Number of pairs "<<pl.getNbrOfPairs()<<std::endl;
      X2=getAMatrixWTriplets(pl,tl,dftPos,nbrOfConfigs);
    }
  else if(doATAT)
    {  
      if(doCV)
	{
	  int cvStructs=40;
	  X2=getAwithATAT(dftPos,nbrOfConfigs,subElements,doubleCO,distances,false);
	  
	  int testSize=20;
	  std::cout<<" #cutoff: "<<doubleCO<< ". number of averaging at every Nconf: "<< cvStructs<< ". testSize: "<<testSize<<std::endl;
	  int fixTrainSize=160;
	  // std::cout<<"#CV as function of numper of pairs calculated for train size of: "<<fixTrainSize<<  "number of averaging at every pair size: "<< cvStructs<< ". testSize: "<<testSize<<std::endl;
	  std::cout<<"#Npairs, cvAverage, variance CV, mseTest, mseTrain, nonzeroParams"<<std::endl;
	  for(int jj=distances.size()+1; jj<=nbrOfConfigs-testSize; jj+=5)
	  //while(distances.size()>=1)
	    {
	      //int jj = fixTrainSize;
	      
	      double cv=0.0;
	      double cvSquare=0.0;
	      double mseTest=0.0;
	      double mseTrain=0.0;
	      double nonZeroParams=0.0;
	       
	      for(int averageCV=0; averageCV<cvStructs; averageCV++)
		{
		  double tempCV=0.0;		  
		  shuffleLists(energies,dftPos);

		  X2=getAwithATAT(dftPos,jj,subElements,doubleCO,distances,false);
		  X=getAwithATAT(dftPos,jj+testSize,subElements,doubleCO,distances,false);
		  params = doMinimize(X2,distances.size(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);
		  std::vector<double> tempEnergy = energyFromParams(params,X2);
		  std::vector<double> cvCorr = getAwithATAT(dftPos,jj,subElements,doubleCO,distances,true);
		  double medel=0;
		  double std=0;

		  for(int jj2=4; jj2<params.size(); jj2++)
		    {
		      medel += fabs(params[jj2]);
		      std += pow(params[jj2],2.0);
		    }
		  medel = medel/((double)params.size()-4);
		  std=sqrt(std/((double)params.size()-4)-pow(medel,2.0));

		
		  for(int jj2=0; jj2<params.size(); jj2++)
		    {		   
		      if(fabs(params[jj2])>(0.5*medel))
			{
			  nonZeroParams++;
			}
		    }

		  for(int jj2=0; jj2<cvCorr.size(); jj2++)
		    {
		      tempCV += pow(((tempEnergy[jj2]-energies[jj2])/(1-cvCorr[jj2])),2.0);
		      mseTrain += pow(tempEnergy[jj2]-energies[jj2],2.0);
		    }
		  tempEnergy= energyFromParams(params,X);
		  for(int jj2=cvCorr.size(); jj2<cvCorr.size()+testSize; jj2++)
		    {
		      mseTest +=pow(tempEnergy[jj2]-energies[jj2],2.0);
		    }		  
		  cv += tempCV/((double)cvCorr.size());
		  cvSquare += pow(tempCV/((double)cvCorr.size()),2.0);

		}
	      //std::cout<<distances.size()<<" "<<cv/cvStructs<< " "<<cvSquare/cvStructs-pow(cv/((double)cvStructs),2.0)<< " "<<mseTest/(testSize*cvStructs)<< " "<<mseTrain/(cvStructs*jj)<<" "<<nonZeroParams/((double)cvStructs)<< std::endl;
	      	      std::cout<<jj<<" "<<cv/cvStructs<< " "<<cvSquare/cvStructs-pow(cv/((double)cvStructs),2.0)<< " "<<mseTest/(testSize*cvStructs)<< " "<<mseTrain/(cvStructs*jj)<<" "<<nonZeroParams/((double)cvStructs)<< std::endl;


		      //doubleCO=distances[distances.size()-1]-1e-4;


	    }      
	}
      X2=getAwithATAT(dftPos,nbrOfConfigs,subElements,doubleCO,distances,false);      

    }
  else
    {
      X2=getAMatrix(pl,dftPos,nbrOfConfigs);
    }
  //PairList pl2 = PairList();
  // pl2.initializePairs(lista,subelements,cutOffVector[CO]);  
  //#pragma omp parallel for private(X,params,calcEnergies,tempNrgy) shared(nrgyValid,nrgyTrain)  
  for(int i=configStart; i<configStop; i+=configStep)
    {
      std::cout<<"Number in training set: "<<i<<" "<<"Number in validation: "<<nbrOfConfigs-i<<std::endl;

      if(doTriplet)
	{
	  //std::cout<<"generating X matrix with triplets. Number of triplets: "<<tl.getNbrOfTriplets()<<std::endl;
	  //X =getAMatrixWTriplets(pl,tl,dftPos,i);
	  X = getAMatrixFromExistingOne(X2,i,pl.getNbrOfPairs()+tl.getNbrOfTriplets());
	  // std::cout<<"Do the minimizing for triplets.."<<std::endl;
	  params = doMinimize(X,pl.getNbrOfPairs()+tl.getNbrOfTriplets(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);
	  //  std::cout<<"after minimzing. Size of params: "<<params.size()<<std::endl;
	  //	  std::cout<<"Number of nonzero params:"<<std::endl;
	}
      else if(doATAT)
	{
	  std::cout<<"Generating x matrix for ATAT "<<distances.size()<<" "<<std::endl;
	  X = getAMatrixFromExistingOne(X2,i,distances.size());
	  std::cout<<"size div by distSize= columns = "<<X2.size()/((double)distances.size())<<std::endl;
	  params = doMinimize(X,distances.size(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);
	  // std::cout<<params.size()<<std::endl;
	  // for(int k2=0; k2<10; k2++)
	  //   {
	  //     std::cout<<"config ======= "<<k2<< std::endl;
	  //     for(int kk=0; kk<params.size(); kk++)
	  // 	{
	  // 	  std::cout<<"distance: "<<distances[kk]<< " energy contribution "<<params[kk]*X[k2*params.size()+kk]<<std::endl;
	  // 	}
	  //   }	  
	}
      else
	{
	  std::cout<<"generating X matrix with doublets. Number of doublets: "<<pl.getNbrOfPairs()<<std::endl;
	  //X =getAMatrix(pl,dftPos,i);
	  X = getAMatrixFromExistingOne(X2,i,pl.getNbrOfPairs());
	  params = doMinimize(X,pl.getNbrOfPairs(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);
	}

      if(printECIParams && i==printParamAtIndex)
	{
	  std::cout<<"Printing parameters to file "<<ECIParamOutFile<<std::endl;
	  std::ofstream outECIParams;
	  outECIParams.open(ECIParamOutFile.c_str());

	  if(!doTriplet && !doATAT)
	    {
	      outECIParams<<pl.getPair(0).getDistance()<< " "<< pl.getPair(0).getSite1()<< " "<<pl.getPair(0).getSite2()<< " "<<params[0];

	      for(int h=1; h<params.size(); h++) // begings with one because the first is outside loop #no empty line at bottom
		{
		  outECIParams<<std::endl;
		  outECIParams<<pl.getPair(h).getDistance()<< " "<< pl.getPair(h).getSite1()<< " "<<pl.getPair(h).getSite2()<< " "<<params[h];
		}
	    }
	  else if(doATAT)
	    {
	      outECIParams<<distances[0]<< " "<<params[0];
	      for(int h=1; h<params.size(); h++) // begings with one because the first is outside loop #no empty line at bottom
		{
		  outECIParams<<std::endl;
		  outECIParams<<distances[h]<< " "<<params[h];
		}
	    }
	  else
	    {
	  outECIParams<<pl.getPair(0).getDistance()<< " "<< pl.getPair(0).getSite1()<< " "<<pl.getPair(0).getSite2()<< " "<<params[0];

	      for(int h=0; h<pl.getNbrOfPairs(); h++)
		{
		  outECIParams<<std::endl;
		  outECIParams<<pl.getPair(h).getDistance()<< " "<< pl.getPair(h).getSite1()<< " "<<pl.getPair(h).getSite2()<< " "<<params[h];
		}
	      for(int h=0; h<tl.getNbrOfTriplets(); h++)
		{
		  outECIParams<<std::endl;
		  outECIParams<<tl.getTriplet(h).getDistance1()<< " " <<tl.getTriplet(h).getDistance2()<< " "  <<tl.getTriplet(h).getDistance3()<< 
		    " "<<tl.getTriplet(h).getSite1()<<" "<<tl.getTriplet(h).getSite2()<<" " <<tl.getTriplet(h).getSite3()<<" " <<params[h+pl.getNbrOfPairs()];
		}

	    }
	}
      std::cout<<"calculating energies from params..."<<std::endl;
      calcEnergies = energyFromParams(params,X2);
      tempNrgy=0;
      
      std::ofstream outTraining;
      if(printTrainingAndValidationEnergy && i==printParamAtIndex)
	{
	  outTraining.open(energyTrainFile.c_str());  
	}
      tempParamNorm=0;
      for(int j=0; j<params.size(); j++)
	{
	  tempParamNorm += fabs(params[j]);
	}
      for(int j=0; j<i; j++)
	{
	  tempNrgy += fabs(calcEnergies[j]-energies[j]);
	  if(printTrainingAndValidationEnergy && i==printParamAtIndex )
	    {
	      outTraining<< energies[j]<<" "<<calcEnergies[j]<<std::endl;
	    }
	}
       nrgyTrain[(i-configStart)/configStep]=(tempNrgy/(double)i);
      //nrgyTrain[(i-configStart)/configStep]=i;

      std::ofstream outValidation;
      if(printTrainingAndValidationEnergy&& i==printParamAtIndex )
	{
	  outValidation.open(energyValidFile.c_str());  
	}
      tempNrgy=0;
      for(int j=i; j<nbrOfConfigs; j++)
	{
	  tempNrgy += pow(calcEnergies[j]-energies[j],2.0)/((double)(nbrOfConfigs-i));
	  if(printTrainingAndValidationEnergy && i==printParamAtIndex)
	    {
	      outValidation<<energies[j]<< " "<<calcEnergies[j]<<std::endl;
	    }
	}
        nrgyValid[(i-configStart)/configStep]=(tempNrgy/(double)(nbrOfConfigs-i));
	std::cout<<"Error for validation set: "<<std::sqrt(tempNrgy)<<" norm of params: "<<tempParamNorm<<" last value of param "<<params[params.size()-1]<< std::endl;

	// nrgyValid[(i-configStart)/configStep]=calcEnergies[0];

    }
      
  // std::ostringstream ss;
  // ss << cutOffVector[CO];
  // std::string test=".\/dataFiles\/" + nrgyFiles + ss.str();
  std::string test="dataFiles/" + nrgyFiles;

  std::cout<<"filnamnet: "<<test<<std::endl;
  //nrgyTrainFile=".\\dataFiles\\" +nrgyTrainFile;

  std::ofstream nrgyV;
  nrgyV.open(test.c_str());  
  for(int i =0; i<nrgyValid.size();i++)
    {
      nrgyV<<i*configStep+configStart<<" "<<nrgyValid[i]<<" "<<nrgyTrain[i]<<std::endl;
    }
  nrgyValid.clear();
  nrgyTrain.clear();
}
// }

  if(doParams)
    {
      std::cout<<"Starting doing parameters"<<std::endl;

      std::string paramFile="dataFiles/" + paramDataFile;
      std::ofstream parFile;
      parFile.open(paramFile.c_str());

      int nonZeroParams;
      int nonZeroParamsTrip;

      for(double i=2.45; i<10; i+=1)
	{

	  nonZeroParams=0;
	  nonZeroParamsTrip=0;
	  PairList pl2 = PairList();
	  pl2.initializePairs(lista,subelements,i);  
	  TripletList tl2 = TripletList();

	  if(doTriplet)
	    {
	      tl2.initializeTriplets(lista,subelements,i);
	  
	      X2=getAMatrixWTriplets(pl2,tl2,dftPos,nbrFitParams);
	      X = getAMatrixWTriplets(pl2,tl2,dftPos,nbrOfConfigs);
	      params = doMinimize(X,pl2.getNbrOfPairs()+tl2.getNbrOfTriplets(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);

	    }
	  else
	    {
	      X2=getAMatrix(pl2,dftPos,nbrFitParams);
	      X = getAMatrix(pl2,dftPos,nbrOfConfigs);
	      params = doMinimize(X,pl2.getNbrOfPairs(),energies,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);
	    }
	  std::cout<<i<<std::endl;
	  for(int k=0; k<params.size(); k++)
	    {
	      if(fabs(params[k])>paramLimit)
		{
		  if(doTriplet)
		    {
		      if(k>pl2.getNbrOfPairs())
			{
			  nonZeroParamsTrip++;
			}
		    }
		  else
		    {
		      nonZeroParams++;
		    }
		}
	    }
	  if(doTriplet)
	    {

	      parFile<<i<< " "<<nonZeroParams<<" "<<nonZeroParamsTrip<< " "<<pl2.getNbrOfPairs()<<" "<<tl2.getNbrOfTriplets()<<std::endl;
	    }
	  else
	    {
	      parFile<<i<< " "<<nonZeroParams<<" "<< " "<<pl2.getNbrOfPairs()<<std::endl;
	    }
	}	  
    }
}










std::vector<double> readEnergies(std::string fileName,int linesToRead)
{
  std::vector<double> ret;  
  std::ifstream in(fileName.c_str());
  if(!in)
    {
      std::cout<<"cannot open file! in: readEnergies"<<std::endl;
    }
  double temp;
  int linesRead=0;
  while(!(in.eof()))
    {
      linesRead++;
      in>>temp;
      ret.push_back(temp);
      if(linesRead>linesToRead)
	{
	  break;
	}
    }
  in.close();
  return ret;
}

  

void shuffleLists(std::vector<double> &energies,std::vector<LatticeList> &dftPos)
{
  //shuffle list
  std::vector<int> ordering;
  for(int i=0; i<energies.size(); i++)
    ordering.push_back(i);

  for(int i=0; i<50*energies.size(); i++)
    {
      int temp =rand()%(energies.size());
      int temp2 = rand()%(energies.size());
      LatticeList templl;
      double tempE;
      while(temp == temp2)
	temp2=rand()%(energies.size());
      int temp3;
      templl=dftPos[temp];
      dftPos[temp]=dftPos[temp2];
      dftPos[temp2]=templl;
      tempE=energies[temp];
      energies[temp]=energies[temp2];
      energies[temp2]=tempE;
    }
}

void shuffleLists(std::vector<double> &energies,std::vector<double> &energies2, std::vector<LatticeList> &dftPos)
{
  //shuffle list
  std::vector<int> ordering;
  for(int i=0; i<energies.size(); i++)
    ordering.push_back(i);

  for(int i=0; i<50*energies.size(); i++)
    {
      int temp =rand()%(energies.size());
      int temp2 = rand()%(energies.size());
      LatticeList templl;
      double tempE;
      double tempE2;
      while(temp == temp2)
	temp2=rand()%(energies.size());
      int temp3;
      templl=dftPos[temp];
      dftPos[temp]=dftPos[temp2];
      dftPos[temp2]=templl;
      tempE=energies[temp];
      tempE2=energies2[temp];
      energies[temp]=energies[temp2];
      energies[temp2]=tempE;
      energies2[temp]=energies2[temp2];
      energies2[temp2]=tempE2;
    }
}




std::vector<LatticeList> readConfig(std::string posFileName,int configs,int nbrOfAtoms,int nbrOfProperties,std::vector<std::string> subElements)
{
 std::vector<LatticeList> ret;
 for(int i=0; i<configs; i++)
   {
     std::ostringstream ss;
     //ss << "confs4/config_"<<i;
     ss << posFileName<<i;
     //LatticeList ll = LatticeList(1,1,1,nbrOfAtoms,nbrOfProperties,ss.str(),subElements);
      LatticeList ll = LatticeList(1,1,1,nbrOfAtoms,ss.str());
     
     ret.push_back(ll);
   }
 return ret; 
}



PairList countPairs(LatticeList ll, PairList pl)
{
  double dr;
  pl.resetCounts();
  Pair tempPair;
  for(size_t i=0; i<ll.getNbrOfSites(); i++)
    {
      for(size_t j=i; j<ll.getNbrOfSites(); j++)
	{
	  if(i==j)
	    {
	      continue;
	    }
	  dr=ll.getDistance(i,j);
	  tempPair=Pair(dr,ll.getSite(i),ll.getSite(j));
	  pl.updatePair(tempPair,false);
	}
    }
  // pl.divideCountByTwo();
  return pl;
}

TripletList countTriplets(LatticeList ll, TripletList tl)
{
  double dr1;
  double dr2;
  double dr3;
  std::vector<double> orderDr;
  std::vector<int> orderIndex;
  orderDr.resize(3);
  orderIndex.resize(3);
  tl.resetCounts(); 

  Triplet tempTriplet;
  //#pragma omp parallel for private(dr1,dr2,dr3,orderDr,orderIndex,tempTriplet,tl) shared(ll)
  for(size_t i=0; i<ll.getNbrOfSites(); i++)
    {
      for(size_t j=i+1; j<ll.getNbrOfSites(); j++)
	{
	  for(size_t k=j+1; k<ll.getNbrOfSites(); k++)
	    {
	      if(i==k || j == k || i==j)
		{
		  continue;
		}
	      dr1=ll.getDistance(i,j);
	      dr2=ll.getDistance(i,k);
	      dr3=ll.getDistance(j,k);
	      orderDr[0]=dr1;
	      orderDr[1]=dr2;
	      orderDr[2]=dr3;
	      orderIndex[0]=i;
	      orderIndex[1]=j;
	      orderIndex[2]=k;
	      sortOrder(orderDr,orderIndex);
	      
	      tempTriplet=Triplet(orderDr[0],orderDr[1],orderDr[2],ll.getSite(orderIndex[0]),ll.getSite(orderIndex[1]),ll.getSite(orderIndex[2]));
	      tl.updateTriplet(tempTriplet,false);
		}		
	    
	}
    }
  return tl;
}

void sortOrder(std::vector<double>  &orderDr, std::vector<int>  &orderIndex)
{
  bool swapped=true;
  double tempDr;
  int tempIndex;
  if(orderDr.size() != orderIndex.size())
    {
      std::cout<<"Error size mismatch between orderDr and ORderIndex"<<std::endl;
      std::cout<<"Aborting sorting.... luser"<<std::endl;
    }
  else{
    while(swapped)
      {
	swapped=false;
	for(int i=0; i<orderDr.size()-1; i++)
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
}
std::vector<int> getPairCounts(PairList pl)
{
  int length=pl.getNbrOfPairs();
  std::vector<int> ret;
  ret.resize(length);
  
  for(int i=0; i<length; i++)
    {
      ret[i]=pl.getPair(i).getCount();
    }
  return ret;
}




std::vector<double> getAMatrix(PairList pl ,std::vector<class LatticeList> dftPos,int nbrOfConfigs)
{
  //std::vector<std::vector<int> > ret;
  std::vector<double> ret2;
  ret2.resize(nbrOfConfigs*pl.getNbrOfPairs());

  if(nbrOfConfigs > dftPos.size())
    {
      std::cout<<"Error: nbrOfConfigs to read is larger that configs input in dftPos"<<std::endl;
    }

  for(int i=0; i< nbrOfConfigs; i++)
    {
      pl=countPairs(dftPos[i],pl);
      for(int j=0; j<pl.getNbrOfPairs(); j++)
	{
	  ret2[i*pl.getNbrOfPairs()+j]=(pl.getPair(j).getCount());
	}
    }
  return ret2;
}


std::vector<double> getAMatrixFromExistingOne(std::vector<double> X2,int nbrOfConfigs,int columns)
{
  std::vector<double> ret;
  ret.resize(nbrOfConfigs*columns);
  for(int i=0; i< ret.size(); i++)
    {
      ret[i]=X2[i];
    }
  return ret;
}



std::vector<double> getAMatrixWith1(PairList pl ,std::vector<class LatticeList> dftPos,int nbrOfConfigs)
{
  //std::vector<std::vector<int> > ret;
  std::vector<double> ret2;
  ret2.resize(nbrOfConfigs*(pl.getNbrOfPairs()+1));

  if(nbrOfConfigs > dftPos.size())
    {
      std::cout<<"Error: nbrOfConfigs to read is larger that configs input in dftPos"<<std::endl;
    }

  for(int i=0; i< nbrOfConfigs; i++)
    {
      pl=countPairs(dftPos[i],pl);
      ret2[i*(pl.getNbrOfPairs()+1)]=1.0;

      for(int j=0; j<pl.getNbrOfPairs(); j++)
	{
	  ret2[i*(pl.getNbrOfPairs()+1)+j+1]=(pl.getPair(j).getCount());
	}
    }
  return ret2;
}





std::vector<double> getAMatrixWTriplets(PairList pl ,TripletList tl, std::vector<class LatticeList> dftPos,int nbrOfConfigs)
{
  //std::vector<std::vector<int> > ret;
  std::vector<double> ret2;
  ret2.resize(nbrOfConfigs*(pl.getNbrOfPairs()+tl.getNbrOfTriplets()));

  if(nbrOfConfigs> dftPos.size())
    {
      std::cout<<"Error: nbrOfConfigs to read is larger that configs input in dftPos"<<std::endl;
    }
   PairList plCopy;
   TripletList tlCopy;

  int index;
  // int i;
  // int j;
  // int k;
  //#pragma omp parallel for shared(ret2,dftPos) //private(tl,pl) shared(ret2)
  for(int i=0; i<nbrOfConfigs; i++)
    {
      //  #pragma omp nowait
       plCopy=countPairs(dftPos[i],pl);
      //    std::cout<<i<<" "<<plCopy.getNbrOfPairs()<<" "<<i+plCopy.getNbrOfPairs()<<std::endl;
       tlCopy=countTriplets(dftPos[i],tl);
      for(int j=0; j<plCopy.getNbrOfPairs(); j++)
	{
	  //std::cout<<i*(plCopy.getNbrOfPairs()+tlCopy.getNbrOfTriplets())+j<<" ";
	  ret2[i*(plCopy.getNbrOfPairs()+tlCopy.getNbrOfTriplets())+j]=(plCopy.getPair(j).getCount());
	}
      //std::cout<<"\n";
      for(int k=0; k<tlCopy.getNbrOfTriplets(); k++)
	{
	  //std::cout<<i*(tlCopy.getNbrOfTriplets()+plCopy.getNbrOfPairs())+k+plCopy.getNbrOfPairs()<<" ";
	  ret2[i*(tlCopy.getNbrOfTriplets()+plCopy.getNbrOfPairs())+k+plCopy.getNbrOfPairs()]=(tlCopy.getTriplet(k).getCount());
	}
      //std::cout<<"\n";

    }
  return ret2;
}



std::vector<double> getA2Matrix(std::vector<double> A, int nbrOfPairs,int nbrOfConfigs)
{
  std::vector<double> ret;

  ret.resize(A.size()/3);

  // 300*nbrofPairs/3
  //add the first and third and subtract the second
  for(int i=0; i<nbrOfConfigs*nbrOfPairs/3; i++)
    {
       ret[i]=A[3*i] - A[3*i+1] + A[3*i+2];
    }
  return ret;
}




void transformToGSLMatrix(std::vector<double> mat,gsl_matrix * ret,int rows, int columns)
{
  if(rows*columns != mat.size())
    {
      std::cout<<"error, rows*columns not equal to size of vector sent in transformToGSLMatrix"<<std::endl;
      std::cout<<"rows: "<<rows<< " columns: "<<columns<< " matrixSize " << mat.size()<<std::endl;
    }
  //gsl_matrix * ret = gsl_matrix_alloc (10, 3);  
  //gsl_matrix * m = gsl_matrix_alloc(10,10);
  //gsl_matrix ret = gsl_matrix_calloc (10, 3);
  for(int i=0; i<rows; i++)
    {
      for(int j=0; j<columns; j++)
  	{
  	  gsl_matrix_set(ret,i,j,mat[i*columns+j]); 
  	}
    }
}

void printTheCV()
{
  std::vector<double> mat;
  mat.push_back(1);
  mat.push_back(2);
  mat.push_back(3);
  mat.push_back(4);
  gsl_matrix * test = gsl_matrix_alloc(2,2);
  
  int rows=2;
  int columns=2;
  transformToGSLMatrix(mat,test,rows,columns);

  for(int i=0; i<rows; i++)
    {
      for(int j=0; j<rows; j++)
	{
	  std::cout<<gsl_matrix_get (test,i,j)<<" ";
	}
      std::cout<<std::endl;
    }
}

std::vector<double> getCVCorrection(std::vector<double> X,int rows, int columns)
{

  std::vector<double> ret;

  if(X.size() != rows*columns)
    {
      std::cout<< "error in getCVCorrection, rows and columns size doesnt match matrix size"<<std::endl;
    }
  //gsl_matrix * test = gsl_matrix_alloc(rows,columns);
  gsl_matrix * C = gsl_matrix_alloc(columns,columns);
  gsl_matrix_set_zero(C);
  //std::vector<double> xdouble(X.begin(),X.end());
  //transformToGSLMatrix(xdouble,test,rows,columns);
  const gsl_matrix_view B = gsl_matrix_view_array(&X[0],rows,columns);
  // gsl_matrix C = gsl_matrix_calloc(columns,columns); // for xTx  
  //int gsl_blas_dgemm(CBLAS_TRANSPOSE_t, CBLAS_TRANSPOSE_t, double, const gsl_matrix*, const gsl_matrix*, double, gsl_matrix*)’
  gsl_blas_dgemm(CblasTrans, CblasNoTrans,1.0, &B.matrix, &B.matrix,0.0, C);
  gsl_matrix * inverse = gsl_matrix_alloc(columns,columns);
  gsl_permutation * perm =gsl_permutation_alloc(columns);  
  int s; //signum for LU decompos (yeeees....)  
  gsl_linalg_LU_decomp(C,perm,&s);
  gsl_linalg_LU_invert(C,perm,inverse);
  gsl_vector * Xi = gsl_vector_alloc(columns);
  gsl_vector * XiT = gsl_vector_alloc(columns);
  // Y= (XTX)^-1 * XiT
  gsl_vector * Y = gsl_vector_alloc(columns);
  // vvResult=Xi * Y = double = the droi... err the correction I am looking for
  double vvResult;
  for(int k=0; k<rows; k++)
    {
      gsl_matrix_get_row(Xi,&B.matrix,k);
      gsl_blas_dgemv(CblasNoTrans,1.0,inverse,Xi,0.0,Y);
      gsl_blas_ddot(Xi,Y,&vvResult);
      //std::cout<<"result: "<<vvResult<<std::endl;
      ret.push_back(vvResult);
    }

  gsl_matrix_free(C);
  gsl_matrix_free(inverse);
  gsl_permutation_free(perm);
  gsl_vector_free(Xi);
  gsl_vector_free(XiT);
  return ret;
}



std::vector<double> doMinimize(std::vector<double> inA, int columns,std::vector<double> nrgy,double alpha,double lambda,double mu,bool doSB,int sbIterMax,double sbTol,int bfgsIters,bool verbal,double bfgsTol)
{  
  int rows = inA.size()/columns;
  gsl_matrix * AtA = gsl_matrix_alloc(columns,columns);
  gsl_matrix_set_zero(AtA); // not sure if necessary
  const gsl_matrix_view B = gsl_matrix_view_array(&inA[0],rows,columns);
  gsl_matrix * A = gsl_matrix_alloc(rows, columns);
  gsl_vector * energies = gsl_vector_alloc(rows);
  gsl_vector * ftA = gsl_vector_alloc(columns);
  gsl_vector * ceParams = gsl_vector_alloc(columns);
  gsl_vector * d = gsl_vector_calloc(columns);
  gsl_vector * b = gsl_vector_calloc(columns);
  gsl_vector * res = gsl_vector_alloc(columns); //

  for(int i=0; i<columns; i++) // set random parameters
    {
      //double tempRand = (rand()/(double)RAND_MAX-0.5)*1e-1;
      gsl_vector_set(ceParams,i,0.0);
    }
  for(int i=0; i<rows; i++) // set random parameters
    {
      gsl_vector_set(energies,i,nrgy[i]);
      for(int j=0; j<columns; j++)
	{
	  gsl_matrix_set(A,i,j,inA[i*columns+j]);
	}
    }  

  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&B.matrix,&B.matrix,0.0,AtA);
  gsl_blas_dgemv(CblasTrans,1.0,&B.matrix,energies,0.0,ftA);  

  // gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&A,&A,0.0,AtA);
  //  gsl_blas_dgemv(CblasTrans,1.0,&A,energies,0.0,ftA);    

 minPar inParams;
  inParams.parA=A;
  inParams.parAtA=AtA;
  inParams.parEnergy=energies;
  inParams.parftA=ftA;
  inParams.parRows=rows;
  inParams.parColumns=columns;
  inParams.parB=b;
  inParams.parD=d;
  inParams.parMu=mu;
  inParams.parLambda=lambda;
  inParams.parAlpha=alpha;

  size_t iter = 0;
  int status;
 
  gsl_vector * df = gsl_vector_alloc(columns); // keep the df vector
  gsl_vector_set_zero(df);
  gsl_multimin_function_fdf my_func;
  gsl_multimin_function_fdf my_funcTest;

      my_funcTest.n = columns;
      my_funcTest.f = my_f;
      //my_funcTest.df = my_df;
      //my_funcTest.fdf = my_fdf;
      my_funcTest.params =(void*) &inParams;

  if(doSB)
    {
      my_func.n = columns;
      my_func.f = my_fSB;
      my_func.df = my_dfSB;
      my_func.fdf = my_fdfSB;
      my_func.params =(void*) &inParams;
    }
  else
    {
      my_func.n = columns;
      my_func.f = my_f;
      my_func.df = my_df;
      my_func.fdf = my_fdf;
      my_func.params =(void*) &inParams;

    }
  double initialValue;// = my_f(ceParams,my_func.params);
  double initialValueSB;// = my_f(ceParams,my_func.params);
  
  const gsl_multimin_fdfminimizer_type *T;
  
  gsl_multimin_fdfminimizer *s;

  T= gsl_multimin_fdfminimizer_vector_bfgs2;
  //T= gsl_multimin_fdfminimizer_steepest_descent;
  //T= gsl_multimin_fdfminimizer_conjugate_fr;
  //T= gsl_multimin_fminimizer_nmsimplex2;

  s = gsl_multimin_fdfminimizer_alloc(T,columns);
  //ss = gsl_multimin_fminimizer_alloc(Ttest,columns);
  double newNorm;
  double oldNorm=0.0;

  if(!doSB)
    {
      sbIterMax=1;
    }

  for(int sbIter=0; sbIter<sbIterMax; sbIter++)
    {
      if(doSB)
	{
      initialValue = my_fSB(ceParams,my_func.params);      
	}
      else
	{
	  initialValue = my_f(ceParams,my_func.params);      
	}
      gsl_multimin_fdfminimizer_set(s,&my_func,ceParams,0.10,bfgsTol); //1e-7
      do
	{  
	  iter++;
	  status = gsl_multimin_fdfminimizer_iterate(s);    
	  if(status)
	    {	      
	      break;
	    }
	  status = gsl_multimin_test_gradient(s->gradient,bfgsTol);  //1e-7
	  if(status == GSL_SUCCESS)
	    {
	      if(verbal)
		{
		  printf("Minimum found at: \n");
		}
	    }
	  //std::cout<<"Current approx min(BFGS) "<<gsl_multimin_fdfminimizer_minimum(s)<<std::endl;
	}
      while(status == GSL_CONTINUE && iter < bfgsIters);
  
      if(!doSB)
	{
	  gsl_vector_memcpy(res,gsl_multimin_fdfminimizer_x(s));
		  
	  //  res = gsl_multimin_fdfminimizer_x(s);
	  double endVal=my_f(res,my_func.params);
	  //std::cout<<"The result: "<<std::endl;
	  if(verbal)
	    {
	      std::cout<<" BFGS iters "<<iter<<" SBiter: "<<sbIter<<" initial value: :"<<initialValue<<" end value "<< endVal<<std::endl;
	      std::cout<<"---------------------------------"<<std::endl;
	    }
	}
      else
	{
	  //status= GSL_CONTINUE;
	  gsl_vector_memcpy(res,gsl_multimin_fdfminimizer_x(s));
	  //  res = gsl_multimin_fdfminimizer_x(s);
	  gsl_vector_memcpy(ceParams,gsl_multimin_fdfminimizer_x(s));	  // u^k+1=u^k
	  
	  
	  double endVal=my_fSB(res,my_func.params);
	  //	  std::cout<<std::endl;
	  shrink(ceParams,my_func.params);
	  //  std::cout<<"after "<<gsl_vector_get(inParams.parD,0)<<" "<<gsl_vector_get(inParams.parD,1)<< std::endl;
	  //	  std::cout<<std::endl;
	  //	  std::cout<<std::endl;
	  updateB(ceParams,my_func.params);
	  //	  std::cout<<"efter(1&2) "<<gsl_vector_get(inParams.parB,0)<<" "<< gsl_vector_get(inParams.parB,1)<<std::endl;
	  //	  std::cout<<std::endl;
	  newNorm=gsl_blas_dasum(ceParams);
	  if(verbal)
	    {
	      std::cout<<"SBiter: "<<sbIter<<" initial value: :"<<initialValue<<" end value "<< endVal<<" Difference: "<< initialValue-endVal<<" New norm "<<newNorm<< " normDiff "<<fabs(newNorm-oldNorm)<<" BFGS iters "<<iter<<std::endl;
	    }
	  if(fabs(newNorm-oldNorm)<sbTol)
	    {
	      break;
	    }
	  oldNorm=newNorm;
	  iter=0;

	  
	  // do d and b iter and change ceParams into  gsl_multimin_fdfminimizer_x(s); I guess

	}
    }// end sbIter


  std::vector<double> returnValues;
  for(int i =0; i<columns; i++)
    {
      returnValues.push_back(gsl_vector_get(res,i));
    }

  // ======================== free memory === dont leak
  gsl_matrix_free(AtA);
  gsl_matrix_free(A);
  gsl_vector_free(energies);
  gsl_vector_free(ftA );
  gsl_vector_free(d);
  gsl_vector_free(b);
  gsl_vector_free(res);
  gsl_vector_free(ceParams);
  gsl_multimin_fdfminimizer_free(s);


  // ====================== freedom
  return returnValues;   
}

// void updateX(const gsl_vector * ceParams, void *params)
// {
//   minPar *inParams = (minPar *) params;
//   gsl_vector * tmpX = gsl_vector_alloc(inParams->parColumns);
//   gsl_vector_memcpy(tmpX,ceParams);

// }

void updateB(const gsl_vector * ceParams, void *params)
{
  minPar *inParams = (minPar *) params;
  gsl_vector * tmpV = gsl_vector_alloc(inParams->parColumns);
  gsl_vector_memcpy(tmpV,ceParams); //x
  gsl_vector_scale(tmpV,inParams->parMu); //mu*x
  gsl_vector_add(inParams->parB,tmpV); //b+mu*x
  gsl_vector_sub(inParams->parB,inParams->parD); //b+mu*x-d
  gsl_vector_free(tmpV);
}

void shrink(const gsl_vector * ceParams, void *params)
{
  minPar *inParams = (minPar *) params;
  double mu = inParams->parMu;
  double lambdaInv = 1.0/inParams->parLambda;
  gsl_vector * tmpX = gsl_vector_alloc(inParams->parColumns);
  gsl_vector * tmpY = gsl_vector_alloc(inParams->parColumns);
  gsl_vector_memcpy(tmpX,ceParams);//x
  gsl_vector_scale(tmpX,mu);//mu*x
  gsl_vector_add(tmpX,inParams->parB); // tmpX= mu*x+b
  double temp;
  for(int i=0; i<inParams->parColumns; i++)
    {
      temp =std::max((fabs(gsl_vector_get(tmpX,i))-lambdaInv),0.0);   

      if(gsl_vector_get(tmpX,i)>0)
	{
	  gsl_vector_set(tmpY,i,temp);
	}
      else
	{
	  gsl_vector_set(tmpY,i,-1.0*temp);
	}
    }  
gsl_vector_memcpy(inParams->parD,tmpY); // does this change d?  yes.

  gsl_vector_free(tmpX);
  gsl_vector_free(tmpY);


}


double my_fSB(const gsl_vector * ceParams, void *params)
{
  
  minPar *inParams = (minPar *) params;
  gsl_vector * Ax = gsl_vector_alloc(inParams->parRows);
  gsl_blas_dgemv(CblasNoTrans,1.0,inParams->parA,ceParams,0.0,Ax);
  gsl_vector_sub(Ax,inParams->parEnergy); //Ax-y
  double term1;
  gsl_blas_ddot(Ax,Ax,&term1);
  term1=0.5*term1;
  double term2;
  gsl_vector * tempV = gsl_vector_alloc(inParams->parColumns);
  gsl_vector * tempV2= gsl_vector_alloc(inParams->parColumns);
  gsl_vector_memcpy(tempV,inParams->parD);
  gsl_vector_memcpy(tempV2,ceParams);
  gsl_vector_scale(tempV2,inParams->parMu);   

  gsl_vector_sub(tempV,inParams->parB); //d-b
  gsl_vector_sub(tempV,tempV2); // d-b-mu*x
  gsl_blas_ddot(tempV,tempV,&term2);  
  term2=term2*0.5*(inParams->parLambda);
  //std::cout<<term1<< " "<<term2<<std::endl;

  gsl_vector_free(Ax);
  gsl_vector_free(tempV);
  gsl_vector_free(tempV2);

  return term1+term2;
}


void my_dfSB(const gsl_vector * ceParams, void *params , gsl_vector * df)
{ 
  minPar *inParams = (minPar *) params;
  
  gsl_blas_dgemv(CblasNoTrans,1.0,inParams->parAtA,ceParams,0.0,df); // N.B transpose here 
  gsl_vector_sub(df,inParams->parftA);

 //so as not to change the paramaeters, might be superfluous because of const in argument
  gsl_vector * tempV = gsl_vector_alloc(inParams->parColumns);
  gsl_vector * tempV2= gsl_vector_alloc(inParams->parColumns);
  gsl_vector_memcpy(tempV,inParams->parD);
  gsl_vector_memcpy(tempV2,ceParams);
  gsl_vector_scale(tempV2,inParams->parMu);

  gsl_vector_sub(tempV,tempV2); //d-mu*x
  gsl_vector_sub(tempV,inParams->parB); //d-mu*x-b
  double scaleTemp=(inParams->parLambda)*(inParams->parMu);
  gsl_vector_scale(tempV,scaleTemp); // lmbda*mu*(d - mu * x - b)

  // AtA*x-Atf-lmbda*mu*(d - mu * x - b)
  gsl_vector_sub(df,tempV); 

  gsl_vector_free(tempV);
  gsl_vector_free(tempV2);
}



void my_fdfSB(const gsl_vector *ceParams, void *params, double * f, gsl_vector * df)
{
  *f =my_fSB(ceParams,params);
  my_dfSB(ceParams,params,df);
}


double my_f(const gsl_vector * ceParams, void *params)
{
  minPar *inParams = (minPar *) params;
  gsl_vector * Ax = gsl_vector_alloc(inParams->parRows);
  gsl_blas_dgemv(CblasNoTrans,1.0,inParams->parA,ceParams,0.0,Ax);
  gsl_vector_sub(Ax,inParams->parEnergy); //Ax-y
  double term1;
  gsl_blas_ddot(Ax,Ax,&term1);
  return term1;

}

void my_df(const gsl_vector * ceParams, void *params , gsl_vector * df)
{ 
  minPar *inParams = (minPar *) params;
  gsl_blas_dgemv(CblasNoTrans,1.0,inParams->parAtA,ceParams,0.0,df); // N.B transpose here 
  gsl_vector_sub(df,inParams->parftA);  
  // return ret;
}
void my_fdf(const gsl_vector *ceParams, void *params, double * f, gsl_vector * df)
{
  *f =my_f(ceParams,params);
  my_df(ceParams,params,df); 
}



std::vector<double> energyFromParams(std::vector<double> inParams, std::vector<double> inA)
{  
  std::vector<double> ret;
  int columns = inParams.size();
  int rows = inA.size()/columns;
  ret.resize(rows);
  // if(columspos.size() !=0)
  //   {
  //     std::cout<<"size mismatch in energyFrom params"<<std::endl;
  //   }
  gsl_matrix * A = gsl_matrix_alloc(rows, columns);
  transformToGSLMatrix(inA,A,rows,columns);
  gsl_vector * params = gsl_vector_alloc(columns);
  for(int i=0; i<columns; i++)
    {
      gsl_vector_set(params,i,inParams[i]);
    }
  gsl_vector * Ai = gsl_vector_alloc(columns);
  double tempResult;
  for(int k=0; k<rows; k++)
    {
      gsl_matrix_get_row(Ai,A,k);
      gsl_blas_ddot(Ai,params,&tempResult);
      //std::cout<<"result: "<<tempResult<<std::endl;
      ret[k]=tempResult;
    }
  //std::cout<<"rows"<<rows<<" size of ret "<<ret.size()<<std::endl;
  gsl_matrix_free(A);
  gsl_vector_free(params);
  gsl_vector_free(Ai);
  return ret;
}


void printCVCorr(std::string confDirectory,int numberOfConfigs,std::vector<std::string> subElements,double cutOff)
{
  PairList pl = PairList();
  std::vector<LatticeList> dftPos = readConfig(confDirectory,numberOfConfigs,46,4,subElements);
  LatticeList lista = LatticeList(1,1,1);
  pl.initializePairs(lista,subElements,cutOff);  

  std::vector<double> X = getAMatrixWith1(pl,dftPos,numberOfConfigs);
  //std::vector<double> X = getAMatrix(pl,dftPos,numberOfConfigs);
  for(int i= 0; i<numberOfConfigs; i++)
    {
      for(int j=0; j<pl.getNbrOfPairs()+1; j++)
	{
	  std::cout<<X[i*(pl.getNbrOfPairs()+1)+j]<<" ";
	}
      std::cout<<"\n";
    }
  std::vector<double> corrs = getCVCorrection(X,numberOfConfigs,pl.getNbrOfPairs()+1);

  for(int i = 0; i<numberOfConfigs; i++)
    {
      std::cout<<corrs[i]<<std::endl;
    }
}


void printCVCorr2(std::string confDirectory,std::string parameterFile,int numberOfConfigs,std::vector<std::string> subElements,double cutOff)
{
  PairList pl = PairList();
  std::vector<LatticeList> dftPos = readConfig(confDirectory,numberOfConfigs,46,4,subElements);
  LatticeList lista = LatticeList(1,1,1);
  pl.initializePairs(lista,subElements,cutOff);  
  ParameterList paramList = ParameterList(parameterFile,0.0000001);
  std::vector<NeighbourList> allNbrList; // this really should be a class    
  std::vector< std::vector<double> > alphas;
  std::vector<double> alpha;
  alpha.resize(lista.getNbrOfSites());
  for(int i=0; i<lista.getNbrOfSites(); i++)
    {
      NeighbourList nl = NeighbourList(i,lista,paramList);      
      allNbrList.push_back(nl);
    }
  
  Neighbour tempNbr= Neighbour();
  for(int i=0; i<pl.getNbrOfPairs(); i++)
    {
      tempNbr.setSite1(pl.getPair(i).getSite1());
      tempNbr.setSite2(pl.getPair(i).getSite2());
      tempNbr.setDistance(pl.getPair(i).getDistance());
      for(int j=0; j<lista.getNbrOfSites(); j++)
	{
	  alpha[j]=0;
	  if(allNbrList[j].isThisMatchingNeighbour(dftPos[0],tempNbr)==1)
	    {
	      alpha[j]=1;
	      std::cout<<j<<" = "<< 1<<std::endl;
	    }
	  else
	    std::cout<<j<<" = "<< 0<<std::endl;
	}
      alphas.push_back(alpha);
    }
}

std::vector<double> getAwithATAT(std::vector<LatticeList> dftPos,int numberOfConfigs,std::vector<std::string> subElements,double cutOff,std::vector<double> & dist,bool doAverage)
{
  const double PI = 3.1415926535897932384626;
  bool doTriplet=false;

  double tripCO=4.6;
  int Mi=subElements.size(); //same notation as Wal09 
	
  PairList pl = PairList();
  //std::vector<LatticeList> dftPos = readConfig(confDirectory,numberOfConfigs);
  //LatticeList lista = LatticeList(1,1,1);
  pl.initializePairs(dftPos[0],subElements,cutOff); 
  // pl.printList();

  ParameterList paramList = ParameterList(pl);
  std::vector<NeighbourList> allNbrList; // this really should be a class    
  for(int i=0; i<dftPos[0].getNbrOfSites(); i++)
    {
      NeighbourList nl = NeighbourList(i,dftPos[0],paramList);      
      allNbrList.push_back(nl);
    }
  bool addDistance=true;
  std::vector<double> distances; //keep track on all individual distances..
  std::vector<double> uniq_dist; //keep only the unique distances..
  for(int i=0; i<pl.getNbrOfPairs(); i++)
    {
      addDistance=true;
      for(int j=0; j<distances.size(); j++)
	{
	  if(distances[j]==pl.getPair(i).getDistance())
	    {
	      addDistance=false;
	    }
	}
      if(addDistance)
	{
	  uniq_dist.push_back(pl.getPair(i).getDistance());
	  for(int m=2; m<= subElements.size(); m++)
	    {
	      for(int t=0; t<m-1; t++)
		{
		  distances.push_back(pl.getPair(i).getDistance());
		}
	    }
	}
    }
  uniq_dist.push_back(0);
  for(int m=2; m<= subElements.size(); m++)
    {
       distances.push_back(0);		
    }

  
  std::sort (distances.begin(),distances.end());
  std::sort (uniq_dist.begin(),uniq_dist.end());



  std::vector<double> X;
  double singletCount;
  for(int confIter=0; confIter<numberOfConfigs; confIter++)
    {   
 
      TripletList t3 = TripletList();
      t3.initializeTriplets(dftPos[confIter],subElements,tripCO);


      singletCount=0;	
      for(int i=0; i<uniq_dist.size(); i++)
	{
	  for(int m=2; m<=subElements.size(); m++)
	    {
	      for(int t=0; t<m-1; t++)
		{

		  //std::cout<<i<< " "<<(m-2)%2<< " "<<t%2<<std::endl;
		  //std::cout<< i<< " "<<m<< " "<<t<<" "<<t%2<<std::endl;
		  
		  singletCount=0;	

		  double totalNbrForDistance=0;
		  double tempAverage=0;
		  double tempTotal=0; //divide this always by two becayse you always get twice in isThisMatchingNEigbhour


		  // for(int jj=0; jj<tempNbrs.size(); jj++)// set the right distances for the neighbours
		  //   {
		  //     tempNbrs[jj].setDistance(distances[i]);
		  //   }

		  for(int j=0; j<dftPos[0].getNbrOfSites(); j++) // loop over all sites on config at confIter
		    {	      
		      if(uniq_dist[i]!=0) // if not singlets
			{	
			  double tempVal=0;
			  int tempT=0;
			  for(int ii=0; ii<subElements.size(); ii++)
			    {
			      for(int jj=ii; jj<subElements.size(); jj++)
				{
				  Neighbour tempNbr = Neighbour(subElements[ii],subElements[jj],1,uniq_dist[i]);
				  tempTotal=(double)allNbrList[j].isThisMatchingNeighbour(dftPos[confIter],tempNbr)/2.0;
				  totalNbrForDistance +=tempTotal;
				  // std::cout<<totalNbrForDistance<<std::endl;

				  tempT=(m/2); //round down aye
				  if(((m-2)%2==0))
				    {
				      tempVal=-cos(2*PI*ii*tempT/(subElements.size()));
				    }
				  else
				    {
				      tempVal=-sin(2*PI*ii*tempT/(subElements.size()));
				    }
				  tempT=((t+2)/2); //round down aye
				  
				  if((t%2==0))
				    {
				      tempVal*=-cos(2*PI*jj*tempT/(subElements.size()));
				    }
				  else
				    {
				      tempVal*=-sin(2*PI*jj*tempT/(subElements.size()));
				    }	
				  tempAverage +=tempTotal*tempVal;
				}
			    }
			}
		      else
			{
			  if(t==0)
			    {
			      for(int jj=0; jj<subElements.size(); jj++) // loop through singlets
				{
				  if(dftPos[confIter].getSite(j)==subElements[jj])
				    {		      
				      int tempT=m/2;
				      if(((m-2)%2==0))
					{
					  singletCount +=-cos(2.0*PI*jj*tempT/(subElements.size()));    
					}
				      else
					{
					  singletCount +=-sin(2.0*PI*jj*tempT/(subElements.size()));
					}
				    }
				}
			    }			
			}
		    }
		  // std::cout<<"pushing back, i="<<i<< " m="<<m<< " t="<<t<<" X.size()="<<X.size()<<std::endl;
		  if(uniq_dist[i] != 0) //for pairs
		    {
		      // std::cout<<uniq_dist[i]<< " "<< tempAverage<< " "<<totalNbrForDistance<< std::endl;
		      if(doAverage)
			{			  
			  X.push_back((double)(tempAverage)/(double)totalNbrForDistance);
			}
		      else
			{
			  X.push_back((double)(tempAverage));///totalNbrForDistance);
			}
		    }
		  else // for singlets
		    {	      
		      if(t==0)
			{
			  if(doAverage)
			    {	
			      X.push_back(singletCount/((double)dftPos[confIter].getNbrOfSites()));
			    }
			  else
			    { 
			      X.push_back(singletCount);		  
			    }
			}
		    }
		}
	    }
	       
	}
      if(doTriplet)
	{

	  t3=countTriplets(dftPos[confIter],t3);      
	  //t3.printList();
	  std::vector<double> tripletCluster = t3.getClusterVector(subElements,tripCO);
       
	  for(int kl=0; kl<tripletCluster.size(); kl++)
	    {
	      //	  std::cout<<tripletCluster[kl]<< " ";
	      X.push_back(tripletCluster[kl]);
	    }
       
	}

    }

  
//    std::cout<<"sizeX "<<X.size()<< " distances*nbrofConfigs= "<<numberOfConfigs*distances.size()<<std::endl;

// for( int j=0; j<numberOfConfigs; j++)
//   {
//     for(int i=0; i<distances.size(); i++)
//       {
// 	std::cout<<X[j*distances.size()+i]<< " ";	  
//        }
//     std::cout<< std::endl;
//   }

  if(doTriplet)
    {
      TripletList t3 = TripletList(dftPos[0],subElements,tripCO);
      std::vector<double> tripletCluster = t3.getClusterVector(subElements,tripCO);
      distances.resize(distances.size()+tripletCluster.size());
    }
dist=distances;
//std::cout<<X.size()/(double)distances.size()<<std::endl;


if(doAverage)
  {
    //  std::cout<<"rows : "<< numberOfConfigs<< " columns "<<distances.size()<< " size X "<< X.size()<<std::endl;
    std::vector<double> cvCorr = getCVCorrection(X,numberOfConfigs,distances.size());
    return cvCorr;
  }
return X;
}

    







void getFileNames(int option, std::string &energyTrainFile, std::string &energyValidFile,std::string &confFileName,std::string &confDirectory,std::vector<std::string> &subElements)
{
  //option 1 BaGaGe
  //option 2 BaGaSi
  //option 3 BaAlSi
  //option 4 BaAlGe
  //option 5 BaAlGaSi4
  //option 6 BaAlGaSi8
  //option 7 BaAlGaSi12
  //option 8 BaAlGaSiAll
  //option 9 BaAlGaSi -{4,8,12}
  // option 10 QUATERNARY
  // option 11 BaAlGaGe

  if(option==1)
    {
      energyTrainFile="dataFiles/energyTrainBaGaGe.data";
      energyValidFile="dataFiles/energyValidBaGaGe.data";

      confFileName = "configs/confBaGaGe/config_0";
      confDirectory="configs/confBaGaGe/config_";

  // confFileName = "configs/confs4/config_0";
  // confDirectory="configs/confs4/config_";


      subElements.push_back("Ga");
      subElements.push_back("Ge");
    }
  else if(option==2)
    {
      energyTrainFile="dataFiles/energyTrainBaGaSi.data";
      energyValidFile="dataFiles/energyValidBaGaSi.data";
      confFileName = "configs/confBaGaSi/config_0";
      confDirectory="configs/confBaGaSi/config_";
      subElements.push_back("Ga");
      subElements.push_back("Si");
    }
  else if(option==3)
    {
      energyTrainFile="dataFiles/energyTrainBaAlSi.data";
      energyValidFile="dataFiles/energyValidBaAlSi.data";
      confFileName = "configs/confBaAlSi/config_0";
      confDirectory="configs/confBaAlSi/config_";
      subElements.push_back("Al");
      subElements.push_back("Si");
    }
  else if(option==4)
    {
      energyTrainFile="dataFiles/energyTrainBaAlGe.data";
      energyValidFile="dataFiles/energyValidBaAlGe.data";
      confFileName = "configs/confBaAlGe/config_0";
      confDirectory="configs/confBaAlGe/config_";
      subElements.push_back("Al");
      subElements.push_back("Ge");
    }
  else if(option==5)
    {
      energyTrainFile="dataFiles/energyTrainBaAlGaSi4.data";
      energyValidFile="dataFiles/energyValidBaAlGaSi4.data";
      confFileName = "configs/confBaAlGaSi4/config_0";
      confDirectory="configs/confBaAlGaSi4/config_";
      subElements.push_back("Al");
      subElements.push_back("Ga");
      subElements.push_back("Si");
    }

  else if(option==6)
    {
      energyTrainFile="dataFiles/energyTrainBaAlGaSi8.data";
      energyValidFile="dataFiles/energyValidBaAlGaSi8.data";
      confFileName = "configs/confBaAlGaSi8/config_0";
      confDirectory="configs/confBaAlGaSi8/config_";
      subElements.push_back("Al");
      subElements.push_back("Ga");
      subElements.push_back("Si");
    }
 else if(option==7)
    {
      energyTrainFile="dataFiles/energyTrainBaAlGaSi12.data";
      energyValidFile="dataFiles/energyValidBaAlGaSi12.data";
      confFileName = "configs/confBaAlGaSi12/config_0";
      confDirectory="configs/confBaAlGaSi12/config_";
      subElements.push_back("Al");
      subElements.push_back("Ga");
      subElements.push_back("Si");
    }
 else if(option==8)
    {
      energyTrainFile="dataFiles/energyTrainBaAlGaSiAll.data";
      energyValidFile="dataFiles/energyValidBaAlGaSiAll.data";
      confFileName = "configs/confBaAlGaSiAll/config_0";
      confDirectory="configs/confBaAlGaSiAll/config_";
      subElements.push_back("Al");
      subElements.push_back("Ga");
      subElements.push_back("Si");
    }
 else if(option==9)
    {
      energyTrainFile="dataFiles/energyTrainBaAlGaSi.data";
      energyValidFile="dataFiles/energyValidBaAlGaSi.data";
      confFileName = "configs/confBaAlGaSi/config_0";
      confDirectory="configs/confBaAlGaSi/config_";
      subElements.push_back("Al");
      subElements.push_back("Ga");
      subElements.push_back("Si");
    }
   else if(option==10)
    {
      energyTrainFile="dataFiles/energyTrainQUAT.data";
      energyValidFile="dataFiles/energyValidQUAT.data";
      confFileName = "configs/confQUAT/config_0";
      confDirectory="configs/confQUAT/config_";
      subElements.push_back("Al");
      subElements.push_back("Ga");
      subElements.push_back("Si");
      subElements.push_back("Ge");
    }
   else if(option==11)
     {
       energyTrainFile="dataFiles/energyTrainBaAlGaGe.data";
       energyValidFile="dataFiles/energyValidBaAlGaGe.data";
       confFileName="configs/confBaAlGaGe/config_0";
       confDirectory="configs/confBaAlGaGe/config_";
       subElements.push_back("Al");
       subElements.push_back("Ga");
       subElements.push_back("Ge");
     }
   else if(option==12)
     {
       energyTrainFile="dataFiles/energyTrainBaGaGepwr.data";
       energyValidFile="dataFiles/energyValidBaGaGepwr.data";
       confFileName="configs/confpowerBGG/config_0";
       confDirectory="configs/confpowerBGG/config_";
       subElements.push_back("Ga");
       subElements.push_back("Ge");
     }



 else
   {
     std::cout<<"No config matching option: "<<option<<std::endl;
   }
}

std::vector<NeighbourList> getNLVector(LatticeList ll, ParameterList pl)
{
  std::vector<NeighbourList> ret;
  //std::cout<<" int get NL "<<std::endl;
  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      NeighbourList nl = NeighbourList(i,ll,pl);
      //nl.printList();
      ret.push_back(nl);
    }
  //std::cout<<"returning with size "<<ret.size()<<std::endl;
  return ret;
}




void doClusterStuff(std::string configFolder,int numberOfConfigs, int numberOfProperties,int numberOfAtoms,std::string outputFolder,double cutoff,int tuplets,std::vector<std::string> subElements)
{
  std::cout<<"Do cluster stuff.."<<std::endl;
  
  //read configs
  //generate X matrix with atat
  // write out CV from soon to be function
  // print parameters for nbr of configs with best cv score
  //profit?

  //std::vector<double> getAwithATAT(std::vector<LatticeList> dftPos,int numberOfConfigs,std::vector<std::string> subElements,double cutOff,std::vector<double> & dist,bool doAverage)
  std::cout<<"Reading configurations...."<<std::endl;
  std::vector<LatticeList> dftPos = readConfig(configFolder,numberOfConfigs,numberOfAtoms,numberOfProperties,subElements);
 std::vector<double> dist;
 
// std::vector<double> X = getAwithATAT(dftPos,numberOfConfigs,subElements,cutoff,dist,false);
 std::cout<<"Print CV data..."<<std::endl;
 printCVdata(dftPos,cutoff,numberOfConfigs,subElements);
 
}


void printCVdata(std::vector<LatticeList> dftPos,double cutoff,int numberOfConfigs,std::vector<std::string> subElements)
{

  int nbrOfAverageSteps=10;
  std::vector<double> params;
  std::vector<double> distances;
  std::vector<double> X2=getAwithATAT(dftPos,numberOfConfigs,subElements,cutoff,distances,false);
  int columns=X2.size()/numberOfConfigs;
  std::vector<double> cv_corr = getAwithATAT(dftPos,numberOfConfigs,subElements,cutoff,distances,true);

  std::vector<double> property;
  property.resize(numberOfConfigs);
  std::vector<double> X2_mini;
  std::vector<double> cv_corr_mini;
  std::vector<double> tempEnergy;

  std::cout<<" Printing CV and standard deviation of CV  score for different properties "<<std::endl;
  std::cout<<" Number of configs in training,"<< nbrOfAverageSteps<<" times averaged "<<std::endl;
  std::cout<<"================================================================="<<std::endl;
  for(int j=0; j<dftPos[0].getNumberOfProperties(); j++)
    {
      for(int jj=0; jj<numberOfConfigs; jj++)
	{
	  property[jj]=(dftPos[jj].getProperty(j));
	}

      

      std::cout<<"#Property "<<j<< " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
      
      for(int i=columns+1; i<numberOfConfigs; i++)
	{
	  double tempCV=0;
	  double cv=0;
	  double cvSquare=0;
	  for(int ii =0; ii<nbrOfAverageSteps; ii++) //some averaging
	    {
	      tempCV=0;
	      X2_mini=getAMatrixFromExistingOne(X2,i,X2.size()/numberOfConfigs);
	      cv_corr_mini=getAMatrixFromExistingOne(cv_corr,i,X2.size()/numberOfConfigs);
	      params = standardParameters(X2_mini,property,columns);
	      tempEnergy = energyFromParams(params,X2_mini);
	      
	      for(int ii2=0; ii2<params.size(); ii2++)
		{
		  tempCV += pow(((tempEnergy[ii2]-property[ii2])/(1.0-cv_corr[ii2])),2.0);
		}
	      cv += tempCV/((double)cv_corr_mini.size());
	      cvSquare += pow(tempCV/((double)cv_corr_mini.size()),2.0);
	      shuffleXMatrix(X2,cv_corr,property);
	    }

	  std::cout<<i<< " "<<sqrt(cv/(double)nbrOfAverageSteps)<< " "<< sqrt(cvSquare/(double)nbrOfAverageSteps-pow(cv/((double)nbrOfAverageSteps),2.0))<<std::endl;
	}
    }
  
  std::cout<<"================================================================="<<std::endl;

}


void shuffleXMatrix(std::vector<double> &X, std::vector<double> &cvCorr, std::vector<double> &property) //
{
  //X=[property.size()==numberOfConfigs][columns == X.size/property.size()]
  
  int columns=X.size()/property.size();
  
  double tempRow;
  double tempProperty;
  for(int i=0; i< columns*2; i++)
    {
      int temp =rand()%(property.size());
      int temp2 = rand()%(property.size());
      while(temp == temp2)
	temp2=rand()%(property.size());
      
      tempProperty=property[temp2];
      property[temp2]=property[temp];
      property[temp]=tempProperty;
      
      tempProperty=cvCorr[temp2];
      cvCorr[temp2]=cvCorr[temp];
      cvCorr[temp]=tempProperty;


      
      for(int j=0; j<columns; j++)
	{
	  tempRow=X[temp2*columns+j];
	  X[temp2*columns+j]=X[temp*columns+j];
	  X[temp*columns+j]=tempRow;
	  

	}


 
      
    }




}


std::vector<double> standardParameters(std::vector<double> X,std::vector<double> property,int columns)
{

  const double alpha=0.1;
  const double mu=0.65;
  //const double mu=0.65;
  const double lambda=100;
  const bool doSB=true;
  const int sbIters=1000;
  const double sbTol=1e-6;
  const int bfgsIters=5000;
  const double bfgsTol=1e-5;
  const bool verbal=false;


  return doMinimize(X,columns,property,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);

}
