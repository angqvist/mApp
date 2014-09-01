#include <tuple>
#include <vector>
#include "optimize.hpp"
#include "dlib/optimization.h"
#include "helperFunctions.hpp"
#include <iostream>
using namespace dlib;

double minObjFunc(std::vector<double> A,std::vector<double> x,std::vector<double> f,double mu, double lambda,std::vector<double> d,std::vector<double> b,std::vector<int> AtA,std::vector<double> ftA, int rows, int columns)
{
  // Ax-f
  std::vector<double>  AxSubf = vectorSubtractVector(matrixDotVector(A,x,rows,columns),f);
  // ||ax-f||^2
  double term1 =0.5*vectorDotVector(AxSubf,AxSubf);
  std::vector<double> temp2 = vectorSubtractVector(d,b);
  temp2=vectorSubtractVector(temp2,scalarOnVector(mu,x));
  
  double term2= 0.5*lambda*vectorDotVector(temp2,temp2);
  
  return (term1+term2);
  
}


std::tuple<std::vector<double>,int,double,double> doOptimize(std::vector<double> A,std::vector<double> x,std::vector<double> f,bool doSplitB,double tolerance,double mu,double lambda,int rows,int columns)
{
  std::vector<double> AtA;
  std::vector<double> ftA;

  AtA=matrixTmatrix(A,rows,columns);
  ftA=vectorDotMatrix(f,A,rows,columns);
  
  std::vector<double> data;
  data.push_back(2.2);

  return std::make_tuple(data,2,3.3,2.2);
}
