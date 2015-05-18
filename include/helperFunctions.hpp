#pragma once
#include<vector>
#include<string>
#include <gsl/gsl_matrix.h>


  // // struct minPar
  // {
  //   gsl_matrix * parA;
  //   gsl_matrix * parAtA;
  //   gsl_vector * parEnergy;
  //   gsl_vector * parftA;
  //   int parRows;
  // };
  



//remove
void shuffleLists(std::vector<double> &,std::vector<class LatticeList> &);
void shuffleLists(std::vector<double> &,std::vector<double> &,std::vector<class LatticeList> &);


std::vector<class LatticeList> readConfig(std::string,int,int,int,std::vector<std::string>);

//add all counting into respective classes
void  countPairs(class LatticeList &,class PairList &);
void  countPairs(class LatticeList &,class PairList &,double &);

void  countTriplets(class LatticeList &,class TripletList &,double &);

std::vector<int> getPairCounts(class PairList &);
// from the pairlist, get the rows of pairs you got, int is the number of configs you want to take from ll vector
std::vector<double> getAMatrixFromExistingOne(std::vector<double> &,int &,int &);
void sortOrder(std::vector<double> &,std::vector<int> &);

void transformToGSLMatrix(std::vector<double> &,gsl_matrix,int &,int &); 
std::vector<double> getCVCorrection(std::vector<double> &,int &,int &);

//consider adding the minimizer into separate file
std::vector<double> doMinimize(std::vector<double> &,int,std::vector<double> &,double,double,double,bool,int,double,int,bool,double);
double my_f(const gsl_vector *, void*);
void my_df(const gsl_vector * , void*,gsl_vector *);
void my_fdf(const gsl_vector *, void*, double *, gsl_vector *);

double my_fSB(const gsl_vector *, void*);
void my_dfSB(const gsl_vector * , void*,gsl_vector *);
void my_fdfSB(const gsl_vector *, void*, double *, gsl_vector *);

std::vector<double> energyFromParams(std::vector<double> &, std::vector<double> &);
void shrink(const gsl_vector *, void*);
void updateB(const gsl_vector *, void*);




std::vector<class NeighbourList> getNLVector(class LatticeList,class ParameterList);



void shuffleXMatrix(std::vector<double> &,std::vector<double> &,std::vector<double> &);



std::vector<double> standardParameters(std::vector<double> &,std::vector<double> &,int &,bool &);

std::vector<double> standardParameters(std::vector<double> &,std::vector<double> &,int &,bool &, double,  double);


std::vector<double> getSingleClusterVector(std::string,std::vector<double>,std::vector<std::string>,int,int,bool);
std::vector<double> getSingleClusterVector(class LatticeList & , std::vector<double> &, std::vector<std::string> &,bool &,bool &);
std::vector<double> getSingleClusterVector(class LatticeList & , class ClusterList &, std::vector<double> &, std::vector<std::string> &,bool &, int);


void getClusterVectors(std::vector<class LatticeList> &, std::vector<double> &,std::vector<double> &,std::vector<std::vector<double> > &,std::vector<double> &, std::vector<std::string> &, int &,  int &,bool &,bool &);

void shuffleFittingObject(std::vector<double> &,std::vector<double> &,std::vector<std::vector<double> > &,int);
void pushToBack(std::vector<double> &,std::vector<double> &,std::vector<std::vector<double> > &);

void printInterfaceParameters(std::vector<double> &,std::vector<double> &,class LatticeList &,std::string &,std::vector<std::string> &);


void printVector(std::vector<double> );
