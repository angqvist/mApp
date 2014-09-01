#include <iostream>
#include <string>
#include <vector>
#include "LatticeList.hpp"
#include "Pair.hpp"
#include "Triplet.hpp"
#include "TripletList.hpp"
#include "ParameterList.hpp"
#include "Neighbour.hpp"
#include "NeighbourList.hpp"
#include "PairList.hpp"
#include <typeinfo>
#include <stdlib.h>
#include <math.h> 
#include "helperFunctions.hpp"
#include <fstream>
#include "MC.hpp"
#include "MC_methods.hpp"
#include <sstream>

//#include "optimize.hpp"


  // LatticeList madList = LatticeList(1,1,1);
  // int nbrMadness=40;
  // vector<LatticeList> llMadness;
  // llMadness.resize(nbrMadness);
  // for(int i=0; i<nbrMadness; i++)
  //   {
  //     madList.setRandomSites(16,"Ga","Ge");
  //     llMadness[i]=madList;
  //   }


using namespace std;
int main()
{

  // LatticeList madList = LatticeList(1,1,1);
  // int nbrMadness=300;
  // vector<LatticeList> llMadness;
  // llMadness.resize(nbrMadness);
  // for(int i=0; i<nbrMadness; i++)
  //   {
  //     madList.setRandomSites(16,"Ga","Ge");
  //     llMadness[i]=madList;
  //   }


  
    //std::string test="confs2/" + "config_" + ss.str();

  // for(int i=0; i<300; i++)
  //   {  
  // std::ostringstream ss;
  // ss << "confs2/config_"<<i;

  // LatticeList l2 = LatticeList(1,1,1,46,ss.str());
  // //l2.setRandomSites(46,"Ga","Ge");
  // LatticeList l1 = LatticeList();
  // //l2.printList();
  // cout<<l2.getEnergy()<<" ";s
  // //cout<<l2.getNbrOfSites()<<endl;








  // cout<<MC_totalEnergy(l2,testNbrList)<<endl;




  //   }

  // return 0;

	LatticeList l3 = LatticeList(1, 1, 1, 9, "configs/confSimple/config_0");

	l3.printList();

	TripletList t3 = TripletList();
	std::vector<string> subE;
	subE.push_back("Al");
	subE.push_back("Ga");
	subE.push_back("Ge");

	t3.initializeTriplets(l3, subE, 40);
	
	//t3.resetCounts();
	t3= countTriplets(l3, t3);
	t3.printList();
	int b;
	cin >> b;

	return 0;




  // LatticeList l2 = LatticeList(1,1,1,46,"configs/confBaAlGe/config_0");
  // cout<<"avg latt "<<l2.getVolume()<<endl;
  // vector<string> subelements;
  // subelements.push_back("Al");
  // subelements.push_back("Ge");
  // PairList pl = PairList();    
  // pl.initializePairs(l2,subelements,10); 
  // PairList pl2 = countPairs(l2,pl);
  // pl2.printList();

  int configStart=95;
  int configMax=101;
  int configStep=1;
  int nbrOfConfigs=101;
  int printParamAtIndex=100;
  std::string posFile="testPosT";
  std::string energyFile="torsdag2";
  double paramCutOff=0.02;
  double doubleCO=10;
  double tripletCO=2.6;
  std::string energyOutFile;//="EnergyBaAlGaSiAll.data";
  std::string paramOutFile;//="paramdataBG300.data";
  int nbrOfConfigsInParams=99;
  bool doNrgNorms=false;
  bool doParams=false;
  
  bool doBandGap=false;
  bool doVolume=false;
  bool doLattice=true; 
  bool printTrainingAndValidationEnergy=true;
  bool doATAT=true;
  bool doCV=false;

  std::string energyTrainFile;//="dataFiles/energyTrainBaAlGaSiAll.data";
  std::string energyValidFile;//="dataFiles/energyValidBaAlGaSiAll.data";
  std::string confFileName;// = "configs/confBaAlGaSiAll/config_0"; //just for initiating
  std::string confDirectory;// ="configs/confBaAlGaSiAll/config_"; //the config_ makes it easier for function readConfigs
  vector<string> subElements;

  //option 1 BaGaGe*
  //option 2 BaGaSi
  //option 3 BaAlSi
  //option 4 BaAlGe
  //option 5 BaAlGaSi4
  //option 6 BaAlGaSi8
  //option 7 BaAlGaSi12  
  //option 8 BaAlGaSiAll*
  //option 9 BaAlGaSi -{4,8,12}
  //option 10 QUATERNARY 1098 max
  // option 11 BaAlGaGe* 299



  getFileNames(4,energyTrainFile,energyValidFile,confFileName,confDirectory,subElements);
  std::string ECIParamOutFile ="ECIParams/basd.param";
  //  confFileName = "configs/confSimple/config_0";
  // confDirectory = "configs/confSimple/config_";

energyTrainFile="dataFiles/ET_55.data";
  energyValidFile="dataFiles/EV_55.data";


  bool printECIParams=false;
   bool doTriplet=false;
   std::vector<LatticeList> dftP;
   //LatticeList llP = LatticeList(1,1,1,9,confFileName);
   // std::vector<int> test = neighbourCount(llP,subElements,0,1.01);
   // std::cout<<subElements[0]<< " "<<subElements[1]<<std::endl;
   // for(int i=0; i<9; i++)
   //   {
   //     std::cout<<llP.getSite(i)<<" ";
   //     for(int j=0; j<subElements.size(); j++)
   // 	 {
   // 	   cout<<test[i*subElements.size()+j]<< " ";
   // 	 }
   //     cout<<endl;
   //   }
   // return 0;
   // //  // std::vector<string> subE;
   // //  // subE.push_back("Al");
   // //  // subE.push_back("Ga");
   // //  // subE.push_back("Ge");
   
   // //  //   dftP.push_back(llP);
   // // //llP.printList();
   // std::vector<double> dists;
  //LatticeList lll = LatticeList(1,1,1,46,confFileName);
   //  std::vector<double> tempp = getAwithATAT(dftP,1,subE,10,dists,false);
   //  return 0;
  //cout<<"test"<<endl;
  //lll.printList();
  //return 0;
  //printCVCorr(confDirectory,300,subElements,14);
  //printCVCorr2(confDirectory,"paramApr",300,subElements,14.0);
  // cout<<llP.getDistance(0,6)<<std::endl;
  //getAwithATAT(dftP,1,subElements,14.0,dists);
  
  // return 0;
       printData(configStart,configMax,configStep,nbrOfConfigs,posFile,energyFile,paramCutOff,energyOutFile,paramOutFile,doNrgNorms,nbrOfConfigsInParams,doTriplet,doParams,doubleCO,tripletCO,doBandGap,doVolume,doLattice,printTrainingAndValidationEnergy,energyTrainFile,energyValidFile,confFileName,confDirectory,ECIParamOutFile,subElements,printECIParams,doATAT,printParamAtIndex,doCV);
     return 0;
  //===================================================================


  //==================================================================================
  // do the calculation for the mixing energy dependence on concentration
   //   ParameterList paraEm = ParameterList("ECIParams/bagasi_E.param",0.0,subElements);
   // ParameterList paraBG = ParameterList("ECIParams/BAGS_BG.params",0.0,subElements);
   // ParameterList paraVOL = ParameterList("ECIParams/BAGS_VOL.params",0.0,subElements);
   // ParameterList paraLAT = ParameterList("ECIParams/BAGS_LAT.params",0.0,subElements);

   // std::vector<NeighbourList> nbrE = getNLVector(

   //ParameterList para = ParameterList("ECIParams/BaAlGaSiAll_params",0.001);
   //para2.printList();
   //return 0;
    // std::cout<<"param size "<<para2.getNbrOfParams()<<std::endl;
   // std::cout<<"nbr sites " <<l2.getNbrOfSites()<<std::endl;

   //  LatticeList minconf = LatticeList(1,1,1,46,"configs/finalconfs/bagasi_conf");
  // std::vector<NeighbourList> nlVectorm = getNLVector(minconf,paraEm);

   // std::vector<NeighbourList> nlVectorm = getNLVector(minconf,paraEm);
   //   std::cout<<MC_totalEnergy(minconf,nlVectorm)<<" "<< minconf.getEnergy()<<std::endl;
   //    return 0;
   //std::cout<<"Mi= "<<subElements.size()<<std::endl;

   //vector<LatticeList> many = readConfig(confDirectory,500);
   // for(int i=0; i<500; i++)
   //   {
   //     std::cout<<"total energy "<<MC_totalEnergy(many[i],nlVector)<< " target "<<many[i].getEnergy()<< " diff "<< MC_totalEnergy(many[i],nlVector)-many[i].getEnergy()<< std::endl;
   //   }
   //nlVector[1].printList();
   //return 0;


   //ParameterList paramList = ParameterList("ECIParams/BaAlGaSiAll_params",0.05);
   // paramList.printList();
   //vector<NeighbourList> testNbrList; // this really should be a class      


 
  // for(int i=0; i<l2.getNbrOfSites(); i++)
  //   {
  //     NeighbourList nl = NeighbourList(i,l2,paramList);
  //     testNbrList.push_back(nl);
  //   }
  //testNbrList[500].printList();
  //return 0;
  // PairList tempPL=PairList();
  // tempPL.initializePairs(l2,subElements,7.7);

  // for(int c=0; c<=16; c++)
  //   {
  //     double tempis=0;
	
  //     for(int average=0; average<1000;average++)
  // 	{
  // 	  l2.setRandomSites(c,16-c,"Al","Ga","Si");
  // 	  tempis+=MC_totalEnergy(l2,testNbrList);
  // 	}

  //     //tempPL.initialise
  //     //cout<<"Nbr of pairs " <<tempPL.getNbrOfPairs()<<std::endl;
  //     cout<<(double)c/16.0<<" "<<tempis/1000.0<<endl;
  //   }
  // return 0;
  //end mixing energy
  // ==================================================================================================


  //Monte carlo part













  MC testMC = MC();
  double hmm;
  double stdE=0;
  double i16;
  double k24;
  double c6;
  double stdk24;
  double stdi16;
  double stdc6;
  double orderP=0;
  std::vector<int> temperature;
  temperature.push_back(3000);
  temperature.push_back(2000);
  temperature.push_back(1000);
  temperature.push_back(800);
  temperature.push_back(600);
  temperature.push_back(400);
  


  int size = 1;
  int c=16;
  int mcSteps=40000;
  int mcThrowAway=1500;



   ParameterList paraE = ParameterList("ECIParams/binParam/BGG_E.param",0.0,subElements);
   ParameterList paraBG = ParameterList("ECIParams/binParam/BGG_BG.param",0.0,subElements);
   ParameterList paraVOL = ParameterList("ECIParams/binParam/BGG_VOL.param",0.0,subElements);
   ParameterList paraLAT = ParameterList("ECIParams/binParam/BGG_LAT.param",0.0,subElements);




  LatticeList l2 = LatticeList(size,size,size,46,confFileName);
     std::string Wtype ="Ga";

  l2.setRandomSites(16,"Ga","Ge");
  
   //  std::vector<NeighbourList> nlVector = getNLVector(l2,para2);

   std::vector<NeighbourList> nbrE = getNLVector(l2,paraE);
   std::vector<NeighbourList> nbrBG = getNLVector(l2,paraBG);
   std::vector<NeighbourList> nbrVOL = getNLVector(l2,paraVOL);
   std::vector<NeighbourList> nbrLAT = getNLVector(l2,paraLAT);
   std::vector<std::vector<NeighbourList> > nlVectors;
   nlVectors.push_back(nbrE);
   nlVectors.push_back(nbrBG);
   nlVectors.push_back(nbrVOL);
   nlVectors.push_back(nbrLAT);
   int conc;
   
   std::vector<double> data;
   //first ensemble average and std adn then wyckoffs site
   data.resize(nlVectors.size()*2+6);

  

  // for(int c=0; c<=16; c++)
  //   {

  cout<<"# "<<mcThrowAway<< " mc step equilibration, concentration "<<c/16.0<< " , size "<< size<< ", "<<mcSteps<<" mcSteps equilibrate"<<std::endl;
  cout<<"T,conc,E, accProb,stdE,BG,stdBG,VOL,stdVOL,LAT,stdLAT,c6,k24,i16,stdc6,stdk24,stdi16)"<<std::endl;
   //l2.setRandomSites(c,16-c,"Al","Ga","Si");

  int i;
  for(int i=1200; i>=0; i -= 60)
    //for(int k=0; k<temperature.size(); k++)
    {
      //i=temperature[k];
      testMC.setTemperature(i);
      hmm= testMC.step(mcThrowAway*46*size*size*size,l2,nlVectors[0]); // equilibrate
      testMC.resetSwapCounts();
      //averageStep(mcSteps,averageStep...)
      hmm=testMC.averageStep(size*size*size*46,mcSteps,l2,nlVectors,data,Wtype);
      //hmm=testMC.averageStep(size*size*size*46,mcSteps,l2,nlVectors,stdE,orderP,c6,k24,i16,stdc6,stdk24,stdi16);
  //data:  E-0, std E-1,BG-2, stdBG-3,VOL-4,stdVOL-5,Lat-6,stdLat-7, c6-8, k24-9, i16-10, stdC6-11, stdk24-12, stdi16-13

      cout<<i<<" "<<c<<" "<<data[0]<<" "<<testMC.getSwapRatio()<<" "<<data[1]<< " "<<data[2] << " "<<data[3] << " "<<data[4] << " "<<data[5] << " "<<data[6] <<" "<<data[7] << " "<< data[8] << " " <<data[9] << " " <<data[10] << " " << data[11]<< " "  << data[12] << " "  <<data[13]<<std::endl;
      // cout<<"fin "<<<<endl;
    }
      // }
  return 0;

}
