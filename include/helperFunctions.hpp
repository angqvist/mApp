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
  

void printData(int,int,int,int,std::string,std::string,double,std::string,std::string,bool,int,bool,bool,double,double,bool,bool,bool,bool,std::string,std::string,std::string,std::string,std::string,
	       std::vector<std::string>,bool,bool,int,bool);

std::vector<double> readEnergies(std::string,int);
void shuffleLists(std::vector<double> &,std::vector<class LatticeList> &);
std::vector<class LatticeList> readConfig(std::string,int,int,int);
class PairList countPairs(class LatticeList,class PairList);
class TripletList countTriplets(class LatticeList,class TripletList);

std::vector<int> getPairCounts(class PairList);
// from the pairlist, get the rows of pairs you got, int is the number of configs you want to take from ll vector
std::vector<double> getAMatrix(class PairList,std::vector<class LatticeList>,int);
std::vector<double> getAMatrixFromExistingOne(std::vector<double>,int,int);
std::vector<double> getAMatrixWith1(class PairList,std::vector<class LatticeList>,int);

std::vector<double> getAMatrixWTriplets(class PairList,class TripletList, std::vector<class LatticeList>,int);
void sortOrder(std::vector<double> &,std::vector<int> &);

std::vector<double> getA2Matrix(std::vector<double>,int,int);
void transformToGSLMatrix(std::vector<double>,gsl_matrix,int,int); 
std::vector<double> getCVCorrection(std::vector<double>,int,int);
std::vector<double> doMinimize(std::vector<double>,int,std::vector<double>,double,double,double,bool,int,double,int,bool,double);
double my_f(const gsl_vector *, void*);
void my_df(const gsl_vector * , void*,gsl_vector *);
void my_fdf(const gsl_vector *, void*, double *, gsl_vector *);

double my_fSB(const gsl_vector *, void*);
void my_dfSB(const gsl_vector * , void*,gsl_vector *);
void my_fdfSB(const gsl_vector *, void*, double *, gsl_vector *);

std::vector<double> energyFromParams(std::vector<double>, std::vector<double>);
void shrink(const gsl_vector *, void*);
void updateB(const gsl_vector *, void*);

void printCVCorr(std::string,int ,std::vector<std::string> ,double );

void printCVCorr2(std::string,std::string,int ,std::vector<std::string> ,double );
std::vector<double> getAwithATAT(std::vector<class LatticeList>,int ,std::vector<std::string> ,double,std::vector<double> &,bool );

void getFileNames(int , std::string &, std::string &,std::string &,std::string &,std::vector<std::string> &);

std::vector<class NeighbourList> getNLVector(class LatticeList,class ParameterList);

void doClusterStuff(std::string,int,int,int,std::string,double,int);
