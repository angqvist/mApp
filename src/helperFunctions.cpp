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
#include <omp.h>
#include "clust.hpp"
#include "QuatupletList.hpp"
#include <iomanip>      // std::setprecision
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
     LatticeList ll = LatticeList(1,1,1,nbrOfAtoms,nbrOfProperties,ss.str(),subElements);
     // LatticeList ll = LatticeList(1,1,1,nbrOfAtoms,ss.str());
     
     ret.push_back(ll);
   }
 return ret; 
}



void countPairs(LatticeList &ll, PairList &pl)
{
  double dr;
  ll.calculate_lookup_table();
  pl.resetCounts();
  Pair tempPair;
  bool addPairWhileUpdating=false;
  std::string site1;
  std::string site2;
  std::string tempSite;
  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      for(int j=i+1; j<ll.getNbrOfSites(); j++)
	{
	 
	  site1=ll.getSite(i);
	  site2=ll.getSite(j);
	  if(site1.compare(site2)>0) //use tuple_remodulator instead you archaic buffoon
	    {
	      tempSite=site1;
	      site1=site2;
	      site2=tempSite;
	    }
	  // tempPair.printPair();
	  
	  dr=ll.fast_distance(i,j); // add cutoff here?
	  tempPair=Pair(dr,site1,site2);
	  pl.updatePair(tempPair,addPairWhileUpdating);
	}
    }
  // pl.divideCountByTwo();
}




void countPairs(LatticeList &ll, PairList &pl, double &cutoff)
{
  double dr;
  ll.calculate_lookup_table();
  pl.resetCounts();
  Pair tempPair;
  bool addPairWhileUpdating=false;

  std::string site1;
  std::string site2;
  std::string tempSite;
  std::vector<double> dists;
  dists.resize(1);
  for(int i=0; i<ll.get_original_atoms_count(); i++)
    {
      for(int j=0; j<ll.getNbrOfSites(); j++)
	{
	  if(i==j)
	    {
	      continue;
	    }
	 
	  dr=ll.fast_distance(i,j);
	  if(dr>cutoff)
	    {
	      continue;
	    }
	 
	  dists[0]=dr;
	  site1=ll.getSite(i);
	  site2=ll.getSite(j);
	  if(site1.compare(site2)>0) //use tuple_remodulator instead you archaic buffoon
	    {
	      tempSite=site1;
	      site1=site2;
	      site2=tempSite;
	    }
	  //dr=dists[d];
	  //dr=ll.fast_distance(i,j); // add cutoff here?
	  // tempPair.printPair();
	  tempPair=Pair(dists[0],site1,site2);
	  if(pl.updatePair(tempPair,addPairWhileUpdating))
	    {
	      std::cout<<"found a new pair while counting..."<<std::endl;
	      tempPair.printPair();
	    }
	      
	}
    }
  // std::cout<<"counting pairs1"<<std::endl;
  // pl.printList();
  // std::cout<<"counting pairs2"<<std::endl;
  //pl.divideCountByTwo();
}





void countTriplets(LatticeList &ll, TripletList &tl, double &cutoff)
{
  bool lexical_sort=false;
  bool addNewTripletWhenCounting=false;
  double dr1;
  double dr2;
  double dr3;
  std::vector<double> orderDr;
  std::vector<int> orderIndex;
  orderDr.resize(3);
  orderIndex.resize(3);
  tl.resetCounts(); 
  ll.calculate_lookup_table();
  Triplet tempTriplet;
  std::vector<std::string> elements;
  elements.resize(3);
  int multiplicity;
  //#pragma omp parallel for private(dr1,dr2,dr3,orderDr,orderIndex,tempTriplet,tl) shared(ll)
  for(int i=0; i<ll.get_original_atoms_count(); i++)
    {
      for(int j=i+1; j<ll.getNbrOfSites(); j++)
	{
	  
	  if(ll.fast_distance(i,j)>cutoff)
	    {
	      continue;
	    }
	  for(int k=j+1; k<ll.getNbrOfSites(); k++)
	    {
	      

	      
	      dr1=ll.fast_distance(i,j);
	      dr2=ll.fast_distance(i,k);
	      dr3=ll.fast_distance(j,k);

	      
	      if(dr2>cutoff || dr3 > cutoff)
		{
		  continue;
		}
	      
	      multiplicity=1;
	      
	      if(j<ll.get_original_atoms_count())
		{
		  multiplicity++;
		}
	      
	      if(k<ll.get_original_atoms_count())
		{
		  multiplicity++;
		}


	      
	      orderDr[0]=dr1;
	      orderDr[1]=dr2;
	      orderDr[2]=dr3;
	      elements[0]=ll.getSite(i);
	      elements[1]=ll.getSite(j);
	      elements[2]=ll.getSite(k);
	      //false for sorting alphabetically
	      tuple_remodulator(orderDr,elements,lexical_sort);

	      tempTriplet=Triplet(orderDr[0],orderDr[1],orderDr[2],elements[0],elements[1],elements[2]);

	      
	      if(tl.updateTriplet(tempTriplet,addNewTripletWhenCounting ,multiplicity))
		{
		  tempTriplet.printTriplet();
		  
		   std::cout<<" not found in tripletList whikle counting"<<std::endl;
		   // tl.printList();
		}
	    }		
	}
    }
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
  else
    {
      while(swapped)
	{

	  swapped=false;
	  for(int i=0; i<2; i++)
	    {
	      if(orderDr[i]>orderDr[i+1])
		{
		  if(i==0)
		    {
		      //swap atoms 2 and 3
		      tempIndex=orderIndex[1];
		      orderIndex[1]=orderIndex[2];
		      orderIndex[2]=tempIndex;
		    }
		  if(i==1)
		    {
		      //swap atoms 1 and 2
		      tempIndex=orderIndex[0];
		      orderIndex[0]=orderIndex[1];
		      orderIndex[1]=tempIndex;

		    }
  
	
		  tempDr=orderDr[i];
		  orderDr[i]=orderDr[i+1];
		  orderDr[i+1]=tempDr;

		  swapped=true;
		}
	    }
	}
    }
}
std::vector<int> getPairCounts(PairList &pl)
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






std::vector<double> getAMatrixFromExistingOne(std::vector<double> &X2,int &nbrOfConfigs,int &columns)
{
  std::vector<double> ret;
  ret.resize(nbrOfConfigs*columns);
  for(int i=0; i< ret.size(); i++)
    {
      ret[i]=X2[i];
    }
  return ret;
}









void transformToGSLMatrix(std::vector<double> &mat,gsl_matrix * ret,int &rows, int &columns)
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

std::vector<double> getCVCorrection(std::vector<double> &X,int &rows, int &columns)
{
  for(int i=0; i<rows; i++)
    {
      for(int j=i+1; j<rows; j++)
	{
	  int equal=0;
	  for(int k=0; k<columns; k++)
	    {
	      if(fabs(X[i*columns+k]-X[j*columns+k])<1e-9)
		{
		  equal++;
		}
	    }
	  if(equal==columns)
	    {
	      std::cout<<"Error: two rows equal: rows"<<i<<", "<<j<<std::endl;
	    }
	}
    }
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
  std::cout<<"LU DECOMP"<<std::endl;
  gsl_linalg_LU_decomp(C,perm,&s);
  // for(int i=0; i<columns; i++)
  //   {
  //     std::cout<<std::endl;
  //     for(int j=0; j<columns; j++)
  // 	{
  // 	  std::cout<<gsl_matrix_get(C,i,j)<< " ";
  // 	}
  //   }
  
  std::cout<<"LU DECOMP DONE1"<<std::endl;
  gsl_linalg_LU_invert(C,perm,inverse);
  std::cout<<"LU DECOMP DONE2"<<std::endl;
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



std::vector<double> doMinimize(std::vector<double> &inA, int columns,std::vector<double> &nrgy,double alpha,double lambda,double mu,bool doSB,int sbIterMax,double sbTol,int bfgsIters,bool verbal,double bfgsTol)
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
      gsl_multimin_fdfminimizer_set(s,&my_func,ceParams,0.1,bfgsTol); //1e-7
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
	      std::cout<<"SBiter: "<<sbIter<<"fSB: "<< " initial value: :"<<initialValue<<" end value "<< endVal<<" Difference: "<< initialValue-endVal<<" New norm "<<newNorm<< " normDiff "<<fabs(newNorm-oldNorm)<<" BFGS iters "<<iter<<std::endl;
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
      //returnValues.push_back(gsl_vector_get(res,i));
      returnValues.push_back(gsl_vector_get(ceParams,i));
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
  gsl_vector * tmpY = gsl_vector_alloc(inParams->parColumns); //vector that replaces d
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
  double term1;// = gsl_blas_dnrm2(Ax);
  //std::cout<<"my_fsb term1 "<<term1;

   gsl_blas_ddot(Ax,Ax,&term1); //||ax-y||^2
  // std::cout<<term1<<" "<<sqrt(term1)<<std::endl;
  // term1=sqrt(term1);
  // if(term1<1e-5)
  //   {
  //     std::cout<<"printing params"<<std::endl;
  //     for(int i=0; i<inParams->parColumns; i++)
  // 	{
  // 	  std::cout<<gsl_vector_get(ceParams,i)<<std::endl;
  // 	}

  //   }
  
  term1=0.5*(term1);
  double term2;
  gsl_vector * tempV = gsl_vector_alloc(inParams->parColumns);
  gsl_vector * tempV2= gsl_vector_alloc(inParams->parColumns);
  gsl_vector_memcpy(tempV,inParams->parD); //tempV=d
  gsl_vector_memcpy(tempV2,ceParams);
  gsl_vector_scale(tempV2,inParams->parMu); //tempV2=mu*x

  gsl_vector_sub(tempV,inParams->parB); //tempV=tempV-b = d-b
  gsl_vector_sub(tempV,tempV2); // tempV=tempV-my*x =d - b - mu*x
  // term2=gsl_blas_dnrm2(tempV);
   gsl_blas_ddot(tempV,tempV,&term2);  // ||d-b-mu*x||^2
  // std::cout<<term2<<" "<<sqrt(term2)<<std::endl;
  // term2=sqrt(term2);
   term2=(term2)*0.5*(inParams->parLambda);
  //std::cout<<term1<< " "<<term2<<std::endl;

  gsl_vector_free(Ax);
  gsl_vector_free(tempV);
  gsl_vector_free(tempV2);
  // std::cout<<"term 1 "<<(term1)<< " term2 " <<(term2)<<std::endl;

  return (term1 + term2);
}


void my_dfSB(const gsl_vector * ceParams, void *params , gsl_vector * df)
{ 
  minPar *inParams = (minPar *) params;
  
  gsl_blas_dgemv(CblasNoTrans,1.0,inParams->parAtA,ceParams,0.0,df); // N.B no transpose here 
  gsl_vector_sub(df,inParams->parftA);
  

  // gsl_vector * Ax = gsl_vector_alloc(inParams->parRows);
  // gsl_blas_dgemv(CblasNoTrans,1.0,inParams->parA,ceParams,0.0,Ax);
  // gsl_vector_sub(Ax,inParams->parEnergy); //Ax-y
  // double term1 = gsl_blas_dnrm2(Ax);  
  // std::cout<<"dfsb term1 "<<term1<<std::endl;
  
  // gsl_vector_scale(df,1.0/term1);

  //second part is -lambda*mu*(d-mu*x-b)/||d-b-mu*x||2

  //stuff 

  // double term2;
  // gsl_vector * tempV = gsl_vector_alloc(inParams->parColumns);
  // gsl_vector * tempV2= gsl_vector_alloc(inParams->parColumns);
  // gsl_vector_memcpy(tempV,inParams->parD); //tempV=d
  // gsl_vector_memcpy(tempV2,ceParams);
  // gsl_vector_scale(tempV2,inParams->parMu); //tempV2=mu*x

  // gsl_vector_sub(tempV,inParams->parB); //tempV=tempV-b = d-b
  // gsl_vector_sub(tempV,tempV2); // tempV=tempV-my*x =d - b - mu*x
  // term2 = gsl_blas_dnrm2(tempV);
  // //  gsl_blas_ddot(tempV,tempV,&term2);  // ||d-b-mu*x||^2
  // term2=1.0/(term2); //this is the thing to divide everytghin in the second part with
  // std::cout<<"dfsb term 2 "<<term2<<std::endl;
  
  
  // gsl_vector_memcpy(tempV,inParams->parD); //tempV=d
  // gsl_vector_memcpy(tempV2,ceParams);
  // gsl_vector_scale(tempV2,inParams->parMu); //tempV2=mu*x
  // gsl_vector_sub(tempV,inParams->parB); //tempV=tempV-b = d-b
  // gsl_vector_sub(tempV,tempV2); // tempV=tempV-my*x =d - b - mu*x
  // gsl_vector_scale(tempV,-term2*inParams->parLambda*(inParams->parMu));

  // gsl_vector_add(df,tempV);

  /// stuff

 //so as not to change the paramaeters, might be superfluous because of const in argument


  //old /working

  gsl_vector * tempV = gsl_vector_alloc(inParams->parColumns);
  gsl_vector * tempV2= gsl_vector_alloc(inParams->parColumns);
  gsl_vector_memcpy(tempV,inParams->parD);
  gsl_vector_memcpy(tempV2,ceParams);
  gsl_vector_scale(tempV2,inParams->parMu);



  gsl_vector_sub(tempV,tempV2); //d-mu*x
  gsl_vector_sub(tempV,inParams->parB); //d-mu*x-b
  //  double term2 = 1.0/gsl_blas_dnrm2(tempV);
  double scaleTemp=(inParams->parLambda)*(inParams->parMu);
  gsl_vector_scale(tempV,scaleTemp); // lmbda*mu*(d - mu * x - b)

  // AtA*x-Atf-lmbda*mu*(d - mu * x - b) 
  gsl_vector_sub(df,tempV); 

  //old
  // gsl_vector_free(Ax);

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
  gsl_blas_dgemv(CblasTrans,1.0,inParams->parA,ceParams,0.0,Ax);
  gsl_vector_sub(Ax,inParams->parEnergy); //Ax-y
  double term1;
  gsl_blas_ddot(Ax,Ax,&term1);
  return term1;

}

void my_df(const gsl_vector * ceParams, void *params , gsl_vector * df)
{ 
  minPar *inParams = (minPar *) params;
  gsl_blas_dgemv(CblasTrans,1.0,inParams->parAtA,ceParams,0.0,df); // N.B transpose here 
  gsl_vector_sub(df,inParams->parftA);  
  // return ret;
}
void my_fdf(const gsl_vector *ceParams, void *params, double * f, gsl_vector * df)
{
  *f =my_f(ceParams,params);
  my_df(ceParams,params,df); 
}



std::vector<double> energyFromParams(std::vector<double> &inParams, std::vector<double> &inA)
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




std::vector<NeighbourList> getNLVector(LatticeList ll, ParameterList pl)
{
  std::vector<NeighbourList> ret;
  //std::cout<<" int get NL "<<std::endl;
  for(int i=0; i<ll.get_original_atoms_count(); i++)
    {
      NeighbourList nl = NeighbourList(i,ll,pl);
      //nl.printList();
      ret.push_back(nl);
    }
  return ret;
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


std::vector<double> standardParameters(std::vector<double> &X,std::vector<double> &property,int &columns,bool &newVerbal)
{
  const double alpha=0.1;
  const double mu=0.0065;
  //const double mu=0.65;
  const double lambda=100;
  const bool doSB=true;
  const int sbIters=1000;
  const double sbTol=1e-6;
  const int bfgsIters=5000;
  const double bfgsTol=1e-5;
  const bool verbal=newVerbal;
  return doMinimize(X,columns,property,alpha,lambda,mu,doSB,sbIters,sbTol,bfgsIters,verbal,bfgsTol);

}


std::vector<double> getSingleClusterVector(std::string fileName,std::vector<double> cutoffs,std::vector<std::string> subElements, int properties, int atoms,bool average)
{
  
  const double PI = 3.1415926535897932384626;

  LatticeList ll = LatticeList(1,1,1,atoms,properties,fileName,subElements);

  // first singlets
  std::vector<double> singletVector;
  double tempAverage;
  int tempT;
  int s1;
  double tempVal;
  for(int m=2; m<=subElements.size(); m++)
    {
      tempAverage=0;
      for(int i=0; i<ll.getNbrOfSites(); i++)
	{
	  for(int ii=0; ii<subElements.size(); ii++)
	    {
	      if(ll.getSite(i)==subElements[ii])
		{
		  s1=ii;
		}
	    }
	  
	  tempT=(m/2); //round down aye
	  if(((m-2)%2==0))
	    {
	      tempVal=-cos(2.0*PI*s1*tempT/(subElements.size()));
	    }
	  else
	    {
	      tempVal=-sin(2.0*PI*s1*tempT/(subElements.size()));
	    }
	  tempAverage +=tempVal;
	}
      if(average)
	{
	  singletVector.push_back(tempAverage/(double)ll.getNbrOfSites());
	}
      else
	{
	   singletVector.push_back(tempAverage);
	}
	  
    }

  if(cutoffs.size() >=1)
    {
      PairList pl = PairList();
      pl.initializePairs(ll,subElements,cutoffs[0]);
      countPairs(ll,pl,cutoffs[0]);
      std::vector<double> pairVector = pl.getClusterVector(subElements,cutoffs[0],average);
      if(pairVector.size()>0)
	{
	  singletVector.insert(singletVector.end(),pairVector.begin(), pairVector.end());
	}
    }

  if(cutoffs.size() >=2)
    {
      TripletList tl = TripletList(ll,subElements,cutoffs[1]);
      countTriplets(ll,tl,cutoffs[1]);
      std::vector<double> tripletVector = tl.getClusterVector(subElements,cutoffs[1],average);
      if(tripletVector.size()>0)
	{
	  singletVector.insert(singletVector.end(),tripletVector.begin(), tripletVector.end());
	}
    }
  return singletVector;
}


//version 2 OVERLOADED
std::vector<double> getSingleClusterVector(LatticeList &ll, std::vector<double> &cutoffs,std::vector<std::string> &subElements,bool &average,bool &ATAT)
{
  std::cout<<"calculating lookup table"<<std::endl;
  ll.calculate_lookup_table();
  std::cout<<"initializing pairs"<<std::endl;
  PairList pl = PairList();
  std::cout<<"initializing pairs2"<<std::endl;

  int nbrOfElements=subElements.size();

  pl.initializePairs(ll,subElements,cutoffs[0]);
  //pl.printList();
  // std::cout<<"number of pairs "<<pl.getNbrOfPairs<<std::endl;
  // for(int i=0; i<pl.getNbrOfPairs(); i++)
  //   {
  //     pl.getPair(i).printPair();
  //   }
  TripletList tl;
  if(cutoffs.size() >=2)
    {
      tl = TripletList(ll,subElements,cutoffs[1]);
      tl.sortTripletList();
    }
  QuatupletList ql;// = QuatupletList();

  if(cutoffs.size() >=3)
    {
      ql.initializeQuatuplets(ll,subElements,cutoffs[2]);
    }

  std::cout<<"done initializing tuplelists"<<std::endl;
  
  const double PI = 3.1415926535897932384626;

  //LatticeList ll= LatticeList(1,1,1,atoms,properties,fileName,subElements);

  // first singlets
  std::vector<double> singletVector;
  singletVector.push_back(ll.get_original_atoms_count()); //zerolets
  double tempAverage;
  int tempT;
  int s1;
  double tempVal;
  std::cout<<"do singlets"<<std::endl;
  if(ATAT)
    {
      std::vector<double> tomt;
      bool reverseSort=true;
      std::vector<std::vector<int> > clust_singlet = symmetric_cluster_function(tomt,nbrOfElements,reverseSort);
      //====
      // std::cout<<"singlets"<<std::endl;
      for(int i=0; i<clust_singlet.size(); i++)
      	{
	  // std::cout<<clust_singlet[i][0]<<std::endl;
      	  tempAverage=0;
      	  for(int j=0; j<ll.get_original_atoms_count(); j++)
      	    {
      	      int s1;
      	      for(int jj=0; jj<subElements.size(); jj++)
      		{
      		  if(ll.getSite(j)==subElements[jj])
      		    {
      		      s1=jj;
      		    }
		}
	      //   std::cout<<clust_singlet.size()<<" "<<clust_singlet[i][0]<<" "<< clusterFunction(subElements.size(),s1,clust_singlet[i][0])<<  " clust singlet "<<std::endl;
	      tempAverage += clusterFunction(nbrOfElements,s1,clust_singlet[i][0]);
	    }
      	    

      	  if(average)
      	    {
      	      singletVector.push_back(tempAverage/((double)ll.get_original_atoms_count()));
      	    }
      	  else
      	    { 
	      singletVector.push_back(tempAverage);
      	    }
      	}
    }
  
  else
    {
      for(int m=0; m<subElements.size(); m++)
	{
	  tempAverage=0;
	  for(int i=0; i<ll.get_original_atoms_count(); i++)
	    {
	      if(ll.getSite(i)==subElements[m])
		{
		  tempAverage++;
		}
	  
	    }
	  if(average)
	    {
	      singletVector.push_back(tempAverage/(double)ll.get_original_atoms_count());
	    }
	  else
	    {
	      singletVector.push_back(tempAverage);
	    }
	  
	}
    }
  //std::cout<<"singlets"<<std::endl;

  std::cout<<"starting with pairs" <<std::endl;
  if(cutoffs.size() >=1)
    {
      // PairList pl = PairList();
      //pl.initializePairs(ll,subElements,cutoffs[0]);
      std::cout<<"counting pairs"<<std::endl;
      countPairs(ll,pl,cutoffs[0]);
      std::cout<<"done counting pairs"<<std::endl;
      if(pl.getNbrOfPairs()>0)
	{

	  if(ATAT)
	    {
	      pl.printList();
      
	      //pl.getPair(2).printPair();      
	      std::vector<double> pairVector = pl.getClusterVector(subElements,cutoffs[0],average);
	      std::cout<<"pairvector size "<<pairVector.size()<<std::endl;
	      //std::cout<<"Number of pairs "<<pl.getNbrOfPairs()<< " pairvector.size: "<<pairVector.size()<<std::endl;
	      if(pairVector.size()>0)
		{
		  singletVector.insert(singletVector.end(),pairVector.begin(), pairVector.end());
		}      
	    }
	  else
	    {

	      for(int i=0; i<pl.getNbrOfPairs(); i++)
		{
		  if(pl.getPair(i).getDistance()<cutoffs[0])
		    {
		      singletVector.push_back(pl.getPair(i).getCount());
		    }

		}

	    }
	}
    }

  //  std::cout<<"starting with triplets"<<std::endl;
  if(cutoffs.size() >=2)
    {
      // TripletList tl = TripletList(ll,subElements,cutoffs[1]);
      countTriplets(ll,tl,cutoffs[1]);
      

      //tl.printList();

      if(tl.getNbrOfTriplets()>0)
	{
	  //ATAT below
	  if(!ATAT)
	    {
	      //	  tl.getTriplet(3).printTriplet();
	      tl.printList();
	      // std::cout<<"-------------------------"<<std::endl;
      

	      std::vector<double> tripletVector = tl.getClusterVector(subElements,cutoffs[1],average);
	      std::cout<<"Number of triplets: "<<tl.getNbrOfTriplets()<< " tripvector "<<tripletVector.size()<<std::endl;
	      if(tripletVector.size()>0)
		{
		  singletVector.insert(singletVector.end(),tripletVector.begin(), tripletVector.end());
		}
	    }
	  else
	    {            
	      for(int i=0; i<tl.getNbrOfTriplets(); i++)
		{
		  if(tl.getTriplet(i).getDistance2()<cutoffs[1])
		    {
		      singletVector.push_back(tl.getTriplet(i).getCount());
		    }
		}
	    }
	}
    }

  if(cutoffs.size() >=3)
    {
      ql.count_quatuplets(ll,cutoffs[2]);


      if(!ATAT)
	{
	  std::vector<double> quatVector = ql.getClusterVector(subElements,cutoffs[2],average);
	  ql.printList();
	  std::cout<<"NUmber of quatuplets "<<ql.getNbrOfQuatuplets()<< " quatvector.size() "<<quatVector.size()<<std::endl;
	  
	  if(quatVector.size()>0)
	    {
	      singletVector.insert(singletVector.end(),quatVector.begin(),quatVector.end());
	    }
	}
      else
	{
	  for(int i=0; i<ql.getNbrOfQuatuplets(); i++)
	    {
	      //ql.getQuatuplets.getDis no distance should be over cutoff	      
	      singletVector.push_back(ql.getQuatuplet(i).getCount());
	    }
	}
    }

  

  // if(false)
  //   {
  //     std::cout<<std::setprecision(6)<<ll.getConcentration(subElements[0])<<" "<<singletVector.size()<<" ";
  //     for(int i=0; i<singletVector.size(); i++)
  // 	{
  // 	  std::cout<<singletVector[i]<< " ";
  // 	}
  //     std::cout<<std::endl;
  //   }
  std::cout<<"normalizing"<<std::endl;
  
  if(!average)
    {
      for(int i=0; i<singletVector.size(); i++)
  	{
  	  singletVector[i] = singletVector[i]/((double)ll.get_original_atoms_count());
  	}
      
    }
  
  return singletVector;
}




//getSingleClusterVector(std::string fileName,std::vector<double> cutoffs,std::vector<std::string> subElements, int properties, int atoms,bool average)


void getClusterVectors(std::vector<LatticeList> &lattice_list, std::vector<double> &X, std::vector<double> &X_avg,std::vector<std::vector<double> > &properties,std::vector<double> &cutoffs, std::vector<std::string> &subelements,int &nbrOfproperties,int &atoms,bool &verbal,bool &ATAT)
{
  std::vector<double> temp;
  //std::vector<double> X_avg;
  if(verbal)
    {
      std::cout<<"Starting to parse "<<lattice_list.size()<<" structures "<<std::endl;
    }
  
  // LatticeList ll= LatticeList(1,1,1,atoms,nbrOfproperties,filenames[0],subelements);


  // //ll.printList();
  // PairList pl = PairList();
  // pl.initializePairs(ll,subelements,cutoffs[0]);
  
  // pl.printList();
  // //std::cout<<"done initialize pairs"<<std::endl;

  // TripletList tl;
  // // std::cout<<"starting initialize triplets"<<std::endl;

  // if(cutoffs.size() >=2)
  //   {
  //     //tl = TripletList(ll,subelements,cutoffs[1]);
  //     tl = TripletList(ll,subelements,cutoffs[1]);
  //     tl.sortTripletList();
  //     //tl.printList();
  //   }

  // //std::cout<<"doing quatupletlIst"<<std::endl;
  // QuatupletList ql;// = QuatupletList();

  // if(cutoffs.size() >=3)
  //   {
  //     //ql=QuatupletList
  //     ql.initializeQuatuplets(ll,subelements,cutoffs[2]);
  //   }
  // LatticeList ll2;
  std::cout<<"starting to count cluster vectors"<<std::endl;
  bool noAverage=false;
  bool doAverage=true;
  for(int i=0; i<lattice_list.size(); i++)
    {
      std::cout<<i<<std::endl;
      //ll= LatticeList(1,1,1,atoms,nbrOfproperties,filenames[i],subelements,ll.getLookupTable());
      //   std::cout<<"Get clustervector: "<<std::endl;
      temp=getSingleClusterVector(lattice_list[i],cutoffs,subelements,noAverage,ATAT);
      X.insert(X.end(),temp.begin(),temp.end());
      temp=getSingleClusterVector(lattice_list[i],cutoffs,subelements,doAverage,ATAT);
      X_avg.insert(X_avg.end(),temp.begin(),temp.end());
      for(int j=0; j<nbrOfproperties; j++)
	{
	  properties[j].push_back(lattice_list[i].getProperty(j));
	}
      lattice_list[i].clear_lookup_table();
    }
  //  std::vector<double> cvCorr = getCVCorrection(X,numberOfConfigs,distances.size());

  //  cv_correction=getCVCorrection(X_avg,filenames.size(),X_avg.size()/filenames.size()); //X_avg,rows,columns
}


      
void shuffleFittingObject(std::vector<double> &X,std::vector<double> &X_avg,std::vector<std::vector<double> > &properties,int n)
{
 
  int columns=X.size()/properties[0].size();
  int rows = properties[0].size();
  double tempRow;
  double tempProperty;
  for(int i=0; i< columns*n; i++)
    {
      int temp =rand()%(rows);
      int temp2 = rand()%(rows);
      while(temp == temp2)
	temp2=rand()%(rows);
      
      for(int j=0; j<properties.size(); j++)
	{
	  tempProperty=properties[j][temp2];
	  properties[j][temp2]=properties[j][temp];
	  properties[j][temp]=tempProperty;
	}
     
      for(int j=0; j<columns; j++)
	{
	  tempRow=X[temp2*columns+j];
	  X[temp2*columns+j]=X[temp*columns+j];
	  X[temp*columns+j]=tempRow;

	  tempRow=X_avg[temp2*columns+j];
	  X_avg[temp2*columns+j]=X_avg[temp*columns+j];
	  X_avg[temp*columns+j]=tempRow;	  
	}      

    }
}


void pushToBack(std::vector<double> &X,std::vector<double> &X_avg,std::vector<std::vector<double> > &properties)
{

  
  int columns=X.size()/properties[0].size();
  int rows = properties[0].size();
  double tempRow;
  double tempProperty;

  for(int i=0; i<properties.size(); i++)
    {
      properties[i].push_back(properties[i][0]);
      properties[i].erase(properties[i].begin());
    }

  for(int i=0; i<columns; i++)
    {
      X.push_back(X[0]);
      X.erase(X.begin());      
      X_avg.push_back(X_avg[0]);
      X_avg.erase(X_avg.begin());
    }
}
      
      
void printInterfaceParameters(std::vector<double> &parameters,std::vector<double> &cutoffs,LatticeList &ll,std::string &fileName,std::vector<std::string> &subelements)
{
  //format is:
  //tuple_order dist0...distn energy
  // tuple order 0 is zerolet, 1 is singlet, 2 is pair, 3 is triplets, 4 is quatuplets
  std::cout<<"parameter length: "<<parameters.size()<<std::endl;
  PairList pl = PairList();
  bool sortAlphabetically=true;
  int numberOfElements= subelements.size();
  pl.initializePairs(ll,subelements,cutoffs[0]);
  TripletList tl;
  if(cutoffs.size() >=2)
    {
      tl = TripletList(ll,subelements,cutoffs[1]);
      tl.sortTripletList();
    }

  QuatupletList ql;
  if(cutoffs.size() >=3)
    {
      ql=QuatupletList(ll,subelements,cutoffs[2]);
    }

  int count=0;
  std::vector<double> pair_dist=pl.getUniqueDistances(cutoffs[0]);
  std::sort (pair_dist.begin(),pair_dist.end());
  std::cout<<"number of pairs "<< pair_dist.size()<<std::endl;
  std::ofstream outF;

  outF.open(fileName.c_str());

  //pritn offset
  outF<<std::setprecision(16)<<0<<" "<<parameters[0]; 
  count++;

  for(int jj=0; jj<subelements.size()-1; jj++)
    {
      outF<<std::endl;      
      outF<<1<< " "<<parameters[count];
      count++;
    }

  
  // int atat_pairs;
  // if(subelements.size()==2)
  //   {
  //     atat_pairs=1;
  //   }
  
  
  
  // else if(subelements.size()==3)
  //   {
  //     atat_pairs=3;
  //   }
  

  if(pl.getNbrOfPairs()>0)
    {
      
      
      for(int i=0; i<pair_dist.size(); i++)
	{
	  std::vector<double> tempPair;
	  tempPair.resize(1);
	  tempPair[0]=pair_dist[i];
	  std::vector<std::vector<int > > cluster_func = symmetric_cluster_function(tempPair,numberOfElements,sortAlphabetically);
	  for(int j=0; j<cluster_func.size(); j++)
	    {
	      outF<<std::endl;
	      outF<<2<< " "<<pair_dist[i]<<" "<<parameters[count];
	      count++;
	    }
	}
    }
	  


  // for(int i=0; i<pair_dist.size(); i++)
  // {
  //   for(int jj=0; jj<atat_pairs; jj++)
  //     {
  // 	outF<<std::endl;
  // 	outF<<pair_dist[i]<< " "<<parameters[count];
  // 	count++;
  //     }
  //   }


  if(cutoffs.size()>=2)
    {
      if(tl.getNbrOfTriplets()>0)
	{
	  
	  std::vector<std::vector<double> > tl_dists=tl.getUniqueDistances(cutoffs[1]);
	  
	  for(int i=0; i<tl_dists.size(); i++)
	    {
	      std::vector<std::vector<int> > cluster_func = symmetric_cluster_function(tl_dists[i],numberOfElements, sortAlphabetically);
	      for(int j=0; j<cluster_func.size(); j++)
		{
		  outF<<std::endl;	 
		  outF<<3<<" "<<tl_dists[i][0]<< " "<<tl_dists[i][1]<< " "<<tl_dists[i][2]<< " "<<parameters[count];
		  count++;

		}
	  
	    }
	}
    }

  if(cutoffs.size()>=3)
    {
      if(ql.getNbrOfQuatuplets()>0)
	{
	  std::vector<std::vector<double> > ql_dists = ql.getUniqueDistances(cutoffs[2]); //sort here?
	  
	  for(int i=0; i<ql_dists.size(); i++)
	    {
	      std::cout<<"======================= ql.dists.size()"<<std::endl;
	      std::cout<<ql_dists.size()<<" cutoff: "<<cutoffs[2]<<std::endl;
	      std::cout<<"======================= end"<<std::endl;
	      std::vector<std::vector<int> > cluster_func = symmetric_cluster_function(ql_dists[i],numberOfElements,sortAlphabetically);
	      for(int j=0; j<cluster_func.size(); j++)
		{
		  outF<<std::endl;
		  outF<<4<< " "<<ql_dists[i][0]<< " "<<ql_dists[i][1]<< " "<<ql_dists[i][2]<< " "<<ql_dists[i][3]<< " "<<ql_dists[i][4]<< " "<<ql_dists[i][5]<<" "<<parameters[count];
		  count++;
		}
	    }
	}
    }
  
  if(parameters.size() != count)
    {
      std::cout<<"ERROR: doesnt write out parameters correctly...."<<std::endl;
    }
  std::cout<<"done writing parameters"<<std::endl;
}

      

  
void printVector(std::vector<double> lista)
{
  for(int i=0; i<lista.size(); i++)
    {
      std::cout<<lista[i]<<std::endl;
    }
}
