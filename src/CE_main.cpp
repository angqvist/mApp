#include <iostream>
#include <fstream>
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
using namespace std;


void readInputCE(char* inputFile)
{
  string configFolder;
  int numberOfConfigs;
  int numberOfAtoms;
  int numberOfElements;
  
  std::vector<std::string> subElements;
  int numberOfProperties;
  string outputFolder;
  double cutoff;
  int tuplets;

  ifstream the_file ( inputFile );
  if ( !the_file )
    {
      cout<<"Could not open file "<< inputFile<<std::endl;
      return;
    }
  else
    { 
      the_file>>configFolder;	
      the_file>>numberOfConfigs;
      the_file>>numberOfAtoms;
      the_file>>numberOfElements;
      for(int i=0; i<numberOfElements; i++)
	{
	  string tempElement;
	  the_file >> tempElement;
	  subElements.push_back(tempElement);
	}
      the_file >> numberOfProperties;	  
      the_file>>outputFolder;	
      the_file>>cutoff;
      the_file>> tuplets;   	
    }

  doClusterStuff(configFolder,numberOfConfigs,numberOfProperties,numberOfAtoms,outputFolder,cutoff,tuplets,subElements);

 // int readingLineNumber=0;
 //      while(the_file.get( x ) )
 // 	{
 // 	  if(readingLineNumber==0)
 // 	    {
 // 	      configFolder<<x;
 // 	    }
 // 	  else if(readingLineNumber==1)
 // 	    {
 // 	      numberOfProperties<<x;
 // 	    }
 // 	  else if(readingLineNumber==2)
 // 	    {
 // 	      outputFolder<<x;
 // 	    }
 // 	  else if(readingLineNumber==3)
 // 	    {
 // 	      cutoff<<x;
 // 	    }
 // 	  else if(readingLineNumber==4)
 // 	    {
 // 	      tuplets<<x;
 // 	    }

 // 	  readingLineNumber++;
	
 // 	}
  
  cout<<configFolder<< " "<<numberOfConfigs<<" "<<numberOfProperties<< " "<<outputFolder<< " "<<cutoff<< " "<<tuplets<<std::endl;




}




int main(int argc, char *argv[])
{
  if(argc == 1)
    {
      cout<<"you need help. Help only available in pro version (339.95$)"<<std::endl;
    }

  cout<<"number of arguments "<<argc<<endl;
  for(int i=1; i<argc; i++)
    {
      readInputCE(argv[i]);
    }


  // if ( argc != 2 ) // argc should be 2 for correct execution
  //   // We print argv[0] assuming it is the program name
  //   cout<<"usage: "<< argv[0] <<" <filename>\n";
  // else {
  //   // We assume argv[1] is a filename to open
  //   ifstream the_file ( argv[1] );
  //   // Always check to see if file opening succeeded
  //   if ( !the_file.is_open() )
  //     cout<<"Could not open file\n";
  //   else {
  //     char x;
  //     // the_file.get ( x ) returns false if the end of the file
  //     //  is reached or an error occurs
  //     while ( the_file.get ( x ) )
  //       cout<< x;
  //   }
  //   // the_file is closed implicitly here
  // }

    return 0;

}
