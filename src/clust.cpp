#include "clust.hpp"
#include <math.h>
#include <iostream>
#include <string>
#define PI 3.14159265

double clusterFunction(int nbr_of_elements,int atom_type,int clustFunction)
{
 //odd: use cos

  if(nbr_of_elements==3)
    {
      if(clustFunction==0)
	{

	  if(atom_type==0)
	    {
	      return 0.980580676;
	    }
	  else if(atom_type==1)
	    {
	      return -1.372812946; // -1.4; //-7*0.2;// -7.0/5.0;
	    }
	  else if (atom_type==2)
	    {
	      return 0.39223227; //.4; //2.0/5.0;
	    }
	}



      // if(atom_type==1 || atom_type==2)
      //   {
      //     return 3.0/((double)sqrt(6.0));
      //   }
      // else
      //   {
      //     return -2.0*3.0/((double)sqrt(6.0));
      //   }
    
      else if(clustFunction==1)
	{

	  if(atom_type==0)
	    {
	      return 1.019049331;  //1.5;
	    }
	  else if(atom_type==1)
	    {
	      return 0.33968311; //0.5;
	    }
	  else if(atom_type==2)
	    {
	      return -1.358732441; //-2.0;
	    }
	}
    }
  

  
	  // if(atom_type==0)
	  //   {
	  //     return -3.0/((double)sqrt(2));
	  //   }
	  // else if(atom_type==1)
	  //   {
	  //     return 3.0/((double)sqrt(2));
	  //   }
	  // else if(atom_type==2)
	  //   {
	  //     return 0.0;
	  //   }
  

  else if(nbr_of_elements==2)
    {
      if(atom_type==0)
	{
	  return -1.0;
	}
      else if(atom_type==1)
	{
	  return 1.0;
	}
    }

  else
    {
      
      if(((clustFunction+2)%2==0))
	{

	  //     std::cout<<"cos "<<clustFunction<< " "<< atom_type<< " "<< nbr_of_elements<<" "<<-cos(2.0*PI * (double)( (int)(clustFunction+2)/2 ) * (double)atom_type/((double)nbr_of_elements))<<  std::endl;
	  return -cos(2.0*PI * (double)( (int)(clustFunction+2)/2 ) * (double)atom_type/((double)nbr_of_elements));
	}
      //even: use sin
      else 
	{
	  //    std::cout<<"sin "<<clustFunction<< " "<< atom_type<< " "<< nbr_of_elements<<" "<<-sin(2.0*PI * (double)( (int)(clustFunction+2)/2 ) * (double)atom_type/((double)nbr_of_elements))<<std::endl;
	  return -sin(2.0*PI * (double)( (int)(clustFunction+2)/2 ) * (double)atom_type/((double)nbr_of_elements));
	}
    }
}




std::vector<std::vector<std::string> > symmetric_cluster_function(std::vector<double> dists,std::vector<std::string> subelements)
{
  bool reverse_sort=false;
  //N.B +1 in number of elements
  std::vector<std::vector<int> > int_base = symmetric_cluster_function(dists,subelements.size()+1,reverse_sort);
  // std::cout<<"int base size[0] "<<int_base[0].size()<<std::endl;
  std::vector<std::vector<std::string> > ret;
  ret.resize(int_base.size());
  for(int i=0; i<int_base.size(); i++)
    {
      for(int j=0; j<int_base[i].size(); j++)
	{
	  //std::cout<<int_base[i][j]<<std::endl;
	  ret[i].push_back( subelements[ int_base[i][j] ] );
	}
    }
  return ret;
}



std::vector<std::vector<int> > symmetric_cluster_function(std::vector<double> dists, int nbr_of_elements,bool reverse_sort)
{


  int tuple_size=dists.size();

  std::vector<std::vector<int> >returnVector;
  
  int cluster_max=nbr_of_elements-2;
  
  if(tuple_size==0)
    {
      std::vector<int> one_ret;
      one_ret.resize(1);
      for(int i=0; i<cluster_max+1; i++)
	{
	  one_ret[0]=i;
	  returnVector.push_back(one_ret);
	}
    }  
  else if(tuple_size==1)
    {
      std::vector<int> two_ret;
      two_ret.resize(2);
      std::vector<std::string> temp_element;
      temp_element.resize(2);
      for(int i=0; i<cluster_max+1; i++)
	{
	  for(int j=0; j<cluster_max+1; j++)
	    {
	      // make combos, only keep symmetric similar in reduced form
	      temp_element[0]=std::to_string(i);
	      temp_element[1]=std::to_string(j);
	      
	      sortElements(temp_element,reverse_sort);
	      two_ret[0] = std::stoi(temp_element[0]);
	      two_ret[1] = std::stoi(temp_element[1]);
	      
	      if(add_new_cluster_controller(returnVector,two_ret))
		{
		  returnVector.push_back(two_ret);
		}

	    }
	}
    }

  
  else if(tuple_size==3)
    {
      std::vector<int> three_ret;
      three_ret.resize(3);
      std::vector<std::string> temp_element;
      std::vector<double> temp_dists=dists;
      temp_element.resize(3);
      
      for(int i=0; i<cluster_max+1; i++)
	{
	  for(int j=0; j<cluster_max+1; j++)
	    {
	      for( int k=0; k<cluster_max+1; k++)
		{
		  temp_dists=dists;
		  temp_element[0]=std::to_string(i);
		  temp_element[1]=std::to_string(j);
		  temp_element[2]=std::to_string(k);
		  
		  tuple_remodulator(temp_dists,temp_element,reverse_sort);

		  three_ret[0]=std::stoi(temp_element[0]);
		  three_ret[1]=std::stoi(temp_element[1]);
		  three_ret[2]=std::stoi(temp_element[2]);
		  
		  if(add_new_cluster_controller(returnVector,three_ret))
		    {
		      returnVector.push_back(three_ret);
		    }
		}
	    }
	}
    }
		  
      



  
  else if(tuple_size==6)
    {
      std::vector<int> four_ret;
      four_ret.resize(4);
      std::vector<std::string> temp_element;
      std::vector<double> temp_dists=dists;
      temp_element.resize(4);

      for(int i=0; i<cluster_max+1; i++)
	{
	  for(int j=0; j<cluster_max+1; j++)
	    {
	      for(int k=0; k<cluster_max+1; k++)
		{
		  for(int l=0; l<cluster_max+1; l++)
		    {
		      temp_dists=dists;
		      temp_element[0]=std::to_string(i);
		      temp_element[1]=std::to_string(j);
		      temp_element[2]=std::to_string(k);
		      temp_element[3]=std::to_string(l);

		      tuple_remodulator(temp_dists,temp_element,reverse_sort);

		      four_ret[0]=std::stoi(temp_element[0]);
		      four_ret[1]=std::stoi(temp_element[1]);
		      four_ret[2]=std::stoi(temp_element[2]);
		      four_ret[3]=std::stoi(temp_element[3]);

		      if(add_new_cluster_controller(returnVector,four_ret))
			{
			  returnVector.push_back(four_ret);
			}
		    }
		}
	    }
	}

    }

  clust_sort_return_vector(returnVector);
  return returnVector;
}


void clust_sort_return_vector(std::vector<std::vector<int> > &return_vector)
{
  bool swapped=true;
  bool previous_columns_equal;
  // std::cout<<"before"<<std::endl;
  // for(int ii=0; ii<return_vector.size(); ii++)
  //   {
  //     for(int jj=0; jj<return_vector[0].size(); jj++)
  // 	{
  // 	  std::cout<<return_vector[ii][jj]<<" ";
  // 	}
  //     std::cout<<std::endl;
  //   }

  // while(swapped)
  //   {
  //     swapped=false;
  //     for(int i=0; i<return_vector.size()-1; i++)
  // 	{
      
  // 	  for(int j=0; j<return_vector[0].size(); j++)
  // 	    {
  // 	      if(fabs(return_vector[i][j]-return_vector[i+1][j])<1e-4)
  // 		{
  // 		  continue;
  // 		}
  // 	      else if(return_vector[i][j]>return_vector[i+1][j])
  // 		{
  // 		  std::cout<<i<<" "<<j<<std::endl;
  // 		  return_vector[i].swap(return_vector[i+1]);
  // 		  for(int ii=0; ii<return_vector.size(); ii++)
  // 		    {
  // 		      for(int jj=0; jj<return_vector[0].size(); jj++)
  // 			{
  // 			  std::cout<<return_vector[ii][jj]<<" ";
  // 			}
  // 		      std::cout<<std::endl;
  // 		    }
	      
		  
  // 		  //swap
  // 		  swapped=true;
  // 		}
  // 	      else
  // 		{
  // 		  break;
  // 		}
  // 	    }
  // 	}
  //   }
  

  
  for(int j=0; j<return_vector[0].size(); j++)
    {

      while(swapped)
  	{
  	  swapped=false;
  	  //sort the j:th column
  	  for(int i=0; i<return_vector.size()-1; i++)
  	    {
  	      previous_columns_equal=true;
  	      for(int jj= 0; jj<j; jj++)
  		{
  		  if(fabs(return_vector[i][jj]-return_vector[i+1][jj])>1e-4)
  		    {
  		      previous_columns_equal=false;
		      
  		    }
  		}
  	      // (j==0 || return_vector[i][j-1]==return_vector[i+1][j-1]) )
	      
  	      if(fabs(return_vector[i][j]-return_vector[i+1][j])>1e-4)
  		{

  		  if(return_vector[i][j]>return_vector[i+1][j] && previous_columns_equal) 
  		    {
  		      return_vector[i].swap(return_vector[i+1]);
  		      //swap
  		      swapped=true;
  		    }
		  
  		}
  	    }
  	}
    }



  // std::cout<< "after: "<<std::endl;

  // for(int ii=0; ii<return_vector.size(); ii++)
  //   {
  //     for(int jj=0; jj<return_vector[0].size(); jj++)
  // 	{
  // 	  std::cout<<return_vector[ii][jj]<<" ";
  // 	}
  //     std::cout<<std::endl;
  //   }
}


void clust_sort_dists(std::vector<std::vector<double> > &dists)
{
  bool swapped=true;
  bool previous_columns_equal;




  // while(swapped)
  //   {
  //     swapped=false;
  //     for(int i=0; i<dists.size()-1; i++)
  // 	{
      
  // 	  for(int j=0; j<dists[0].size(); j++)
  // 	    {
  // 	      if(fabs(dists[i][j]-dists[i+1][j])<1e-4)
  // 		{
  // 		  continue;
  // 		}
  // 	      else if(dists[i][j]>dists[i+1][j])
  // 		{
  // 		  // std::cout<<i<<" "<<j<<std::endl;
  // 		  // for(int ii=0; ii<dists.size(); ii++)
  // 		  //   {
  // 		  //     for(int jj=0; jj<dists[0].size(); jj++)
  // 		  // 	{
  // 		  // 	  std::cout<<dists[ii][jj]<<" ";
  // 		  // 	}
  // 		  //     std::cout<<std::endl;
  // 		  //   }
  // 		  dists[i].swap(dists[i+1]);
				  
		  
  // 		  //swap
  // 		  swapped=true;
  // 		}
  // 	      else
  // 		{
  // 		  break;
  // 		}
  // 	    }
  // 	}
  //   }
  









  for(int j=0; j<dists[0].size(); j++)
    {

      while(swapped)
  	{
  	  swapped=false;
  	  //sort the j:th column
  	  for(int i=0; i<dists.size()-1; i++)
  	    {
  	      previous_columns_equal=true;
  	      for(int jj= j-1; jj>=0; jj--)
  		{
  		  if(dists[i][jj] != dists[i+1][jj])
  		    {
  		      previous_columns_equal=false;
		      
  		    }
  		}
  	      // (j==0 || return_vector[i][j-1]==return_vector[i+1][j-1]) )
	      
  	      if(dists[i][j]>dists[i+1][j] && previous_columns_equal) 
  		{
  		  dists[i].swap(dists[i+1]);
  		  //swap
  		  swapped=true;
  		}
  	    }
  	}
    }
}

bool add_new_cluster_controller(std::vector<std::vector<int> > bigger_vector, std::vector<int> prospect)
{
  if(bigger_vector.size()==0)
    {
      return true;
    }
  
  
  for(int i=0; i<bigger_vector.size(); i++)
    {
      int bool_sum=0;
      for(int j=0; j<prospect.size(); j++)
	{
	  
	  if(prospect[j] == bigger_vector[i][j])
	    {
	      bool_sum++;
	    }
	  if(bool_sum==prospect.size())
	    {
	      return false;
	    }
	}
    }
  return true;

}
	      




//this function takes distances and elements and make it into a symmetric equivalent reduced form if possible
// the distances and elements will be smallest first and alphabetical order
void tuple_remodulator(std::vector<double> &dists, std::vector<std::string> &elements,bool reverseSort)
{ 


  if((dists.size()==1 && 2!= elements.size()) || 
     (dists.size()==0 && elements.size()!=1) ||
     (dists.size()==3 && elements.size() != 3) ||
     (dists.size()==6 && elements.size() != 4))
    {
      std::cout<<"distance/elements mismatch in tuple_remodulator"<<std::endl;
      std::cout<<"Dists size: "<<dists.size()<< " elements size: "<<elements.size()<<std::endl;
      return;
    }
  int d_size = dists.size();
  
  if(d_size==1)
    {
      sortElements(elements,reverseSort);
    }
  if(d_size==3)
    {
      clust_sort_triplet(dists,elements,reverseSort);
    }
  if(d_size==6)
    {
      clust_sort_quatuplet(dists,elements,reverseSort);
    }
  
}


//this routine assumes that sites 1 and 2 are ordered and thus dists[0]
//it will then order site 3 and 4
//swapping site 3 and 4 will result in:
// old dist --> new dist
// r2 --> r3
// r4 -->r5
// r5--r4
// r6 -->r6
// so if r2 != r3 it is obvious what to do, else you have to compare higher order distances

void clust_sort_quatuplet_part2(std::vector<double>  &dists, std::vector<std::string> &elements,bool reverseSort)
{

  double reverseInt=1.0;
  if(reverseSort)
    {
      reverseInt=-1.0;
    }

  
  if(fabs(dists[1]-dists[2])<1e-3)
    {
      //r2==r3, check higher order
      //
      if(fabs(dists[3]-dists[4])<1e-3)
	{
	  // doesnt matter
	  //swap sites yo
	  if(reverseInt*elements[2].compare(elements[3])>0)
	    {
	      clust_swap_atom(3,4,elements,dists);
	    }
	  return;
	}
      else if(dists[3]<dists[4])
	{
	  //do not swap its oki doki
	  return;
	}
      else
	{
	  //swap it
	  clust_swap_atom(3,4,elements,dists);
	  return;
	}
    }
  else if(dists[1]<dists[2])
    {
      // do not swap, its oki doki
      return;
    }
  else
    {
      clust_swap_atom(3,4,elements,dists);
      return;
    }
  
}
//return true/false if dists1 is smaller/larger where first distances are compared first and if they are equal check higher order distance
//returns 2 if they are equal

int is_first_dist_smaller(std::vector<double> dists1, std::vector<double> dists2)
{
  if(dists1.size() != dists2.size())
    {
      std::cout<<"Error: dists1 and dists2 size doesn't match in 'is_first_dist_smaller' in clust.cpp"<<std::endl;
      exit(1);
    }
  int equal_count=0;
  int m=dists1.size();
  for(int i=0; i<m; i++)
    {
      if( fabs(dists1[i]-dists2[i])<1e-3 )
	{
	  equal_count++;
	}
      else if(dists1[i]-dists2[i]<1e-3)
	{
	  return true;
	}
      else if(dists1[i]-dists2[i]>1e-3)
	{
	  return false;
	}

    }

  if(equal_count==m)
    {
      return 2;
    }
  return false;
}
  
int trial_swap(std::vector<int> swap_order, std::vector<double> &temp_dists, std::vector<double> &min_dists, std::vector<std::string> &min_elements, std::vector<std::string> &temp_elements,bool reverseSort)
{
  if(swap_order.size()%2 !=0)
    {
      std::cout<<"Error: Swap order.size() is not a multiple of two in trial_swap in clust.cpp"<<std::endl;
      exit(1);
    }
  int smaller_dist_variable;

  for(int i=0; i<swap_order.size()/2; i++)
    {
      clust_swap_atom(swap_order[i*2],swap_order[i*2+1],temp_elements,temp_dists);
    }
  clust_sort_quatuplet_part2(temp_dists,temp_elements,reverseSort);

  smaller_dist_variable=is_first_dist_smaller(temp_dists,min_dists);
  
  if(smaller_dist_variable == true)
    {
      min_dists=temp_dists;
      min_elements=temp_elements;
      return true;
    }
  else if(smaller_dist_variable == 2)
    {
      if(is_elements_in_lower_order(temp_elements,min_elements,reverseSort))
	{
	  min_dists=temp_dists;
	  min_elements=temp_elements;
	  return true;
	}
      //check if elements get swapped
    }
  
  return false;
}


int is_elements_in_lower_order(std::vector<std::string> str1, std::vector<std::string> str2,bool reverseSort)
{
  
  if(str1.size() != str2.size())
    {
      std::cout<<"Error: str1 size and str2 size is not equal in is_elements_in_lower_order in clust.cpp"<<std::endl;
      exit(1);
    }

  double reverseInt=1.0;
  if(reverseSort)
    {
      reverseInt=-1.0;
    }
  int m= str1.size();
  for(int i=0; i<m; i++)
    {
      if(reverseInt*str1[i].compare(str2[i])<0)
	{
	  return true;
	}
      else if(reverseInt*str1[i].compare(str2[i])>0)
	{
	  return false;
	}
    }
  return 2;
  
}
  

  
  
  


void clust_sort_quatuplet(std::vector<double> &dists, std::vector<std::string> &elements,bool reverseSort)
{

  //new -more thought out stuff
  double min_dist=1e10; //==INTMAX or something
  
  std::vector<int> index_of_minimum_distances;
  std::vector<int> number_of_min_dists;
  for(int i=0; i<dists.size(); i++)
    {
      if(dists[i]<min_dist)
	{
	  min_dist=dists[i];
	}      
    }
  for(int i=0; i<dists.size(); i++)
    {
      if(fabs(dists[i]-min_dist)<1e-3)
	{
	  number_of_min_dists.push_back(i);
	}
    }

  //do part one which consists of finding site one, site two and also r1=dists[0]
  // number of candidates.
  //note that if m=1 there are two candidates for site1 either one of the ones involved in 
  int m=number_of_min_dists.size();
  // std::cout<<"================================+"<<std::endl;
  // for(int i=0; i<dists.size(); i++)
  //   {
  //     std::cout<<dists[i]<< " ";
  //   }
  //std::cout<<std::endl;
  std::vector<int> min_indices;
  min_indices.push_back(1);
  min_indices.push_back(2);
  std::vector<double> min_dists=dists;
  std::vector<std::string> min_elements=elements;
  std::vector<double> temp_dists=dists;

  std::vector<std::string> temp_elements=elements;
  for(int i=0; i<m; i++)
    {
      int j = number_of_min_dists[i];
      std::vector<int> swap_atom_order;
      temp_elements=elements;
      temp_dists=dists;
      
      if(j==0)
	{

	  //first try no swaps
	  // and then try 1,2 swap
	  swap_atom_order.clear();
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=1;
	      min_indices[1]=2;
	    }

	  swap_atom_order.clear();
	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(2);
	  
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=2;
	      min_indices[1]=1;
	    }
	}
      //if j==1, r2 is minimum candidate, meaning s1 is s1 and s2 is s3 or s1 is s3 and s2 is s1
      else if(j==1)
	{
	  // first case swap 2,3
	  swap_atom_order.clear();
	  swap_atom_order.push_back(2);
	  swap_atom_order.push_back(3);

	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=1;
	      min_indices[1]=3;
	    }	  
	  //------------------------------------
	  //now test swapping 1,2 (which is like putting s3 to s1 and s1 to s2 in the original configuration
	  swap_atom_order.clear();
	  swap_atom_order.push_back(2);
	  swap_atom_order.push_back(1);
	  
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=3;
	      min_indices[1]=1;
	    }
	}
      //if j==2 it is either swapping 2,4 or swapping 2,4 and 1,2. That is site 1 and site 4 need to occupy either site 1 or 2.
      else if(j==2)
	{
	  //first swap 2,4

	  swap_atom_order.clear();
	  swap_atom_order.push_back(2);
	  swap_atom_order.push_back(4);


	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=1;
	      min_indices[1]=4;
	    }	   

	  //----------------------
	  //now swap 1,2 that is put s4,s1 in s1,s2
	  swap_atom_order.clear();
	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(2);
	  
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=4;
	      min_indices[1]=1;
	    }	  
	  
	}

      //if j==3 s2 and s3 is candidates
      //first put s3 on s1
      //then swap s1,s2
      
      else if(j==3)
	{

	  
	  swap_atom_order.clear();
	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(3);


	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=3;
	      min_indices[1]=2;
	    }	   

	  //----------------------
	  //now swap 1,2 that is put s4,s1 in s1,s2
	  swap_atom_order.clear();
	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(2);
	  
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=2;
	      min_indices[1]=3;
	    }
	}

      //s2 and s4 are candidates
      //first trial_swap s1,s4 and then s1,s2
      else if(j==4)
	{
	  swap_atom_order.clear();
	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(4);


	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=4;
	      min_indices[1]=2;
	    }	   

	  //----------------------
	  //now swap 1,2 that is put s4,s1 in s1,s2
	  swap_atom_order.clear();

	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(2);
	  
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=2;
	      min_indices[1]=4;
	    }


	}
      //s3 and s4 are candidates
      //first swap s3,s1 and s2,s4 (dual swap action)
      //and then the usual s1,s2 swap
      else if(j==5)
	{

	  
	  swap_atom_order.clear();
	  swap_atom_order.push_back(3);
	  swap_atom_order.push_back(1);
	  
	  swap_atom_order.push_back(4);
	  swap_atom_order.push_back(2);


	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=3;
	      min_indices[1]=4;
	    }	   

	  //----------------------
	  //now swap 1,2 that is put s4,s1 in s1,s2
	  swap_atom_order.clear();
	  swap_atom_order.push_back(1);
	  swap_atom_order.push_back(2);
	  
	  if(trial_swap(swap_atom_order,temp_dists,min_dists,min_elements,temp_elements,reverseSort))
	    {
	      min_indices[0]=4;
	      min_indices[1]=3;
	    }
	}
      
    }
  
  dists=min_dists;
  elements=min_elements;
  //do part 2, final part..

  // for(int i=0; i<6; i++)
  //   {
  //     std::cout<<dists[i]<< " ";
  //   }
  // std::cout<<std::endl;

  clust_sort_quatuplet_part2(dists,elements,reverseSort);
  
  // for(int i=0; i<dists.size(); i++)
  //   {
  //     std::cout<<dists[i]<< " ";
  //   }
  // std::cout<<std::endl;
  // std::cout<<"================================+"<<std::endl;

  
  
  
  //end new .more thought out stuff
  
  
  
  
  
  
  
  /*
    a quatuplet is equal to [r1,r2,r3,r4,r5,r6] and [s1,s2,s3,s4]
    r1 = s1,s2
    r2 = s1,s3
    r3 = s1,s4
    r4 = s2,s3
    r5 = s2,s4
    r6 = s3,s4

    swapping s1<->s2:
    r2 <-> r4
    r3 <->r5

    swapping s2<->s3
    r1 <-> r2
    r5 <-> r6
    
    swapping s3 <-> s4
    r2 <-> r3
    r4 <-> r5
    
   


    swapping 1<->3 = 2<->3 + 1<->2 + 2<->3

    General strategy is to make first distance the smallest, r2<r3 < r4 and then the other two distances follow automatically
    


    if r2 < r1
    swap s2,s3


   */
  /*
  double DISTANCE_LIMIT=1e-4;

  bool swapped = true;
  double temp_distance;
  std::string temp_element;
  //first loop makes sure smallest distance is first
  int iterera=0;


  while(swapped)
    {
      swapped=false;
      iterera++;
      if(iterera>10)
	{
	  std::cout<<"iterear: "<<iterera<<std::endl;
	  for(int k=0; k<dists.size(); k++)
	    {
	      std::cout<<dists[k]<< " ";
	    }
	  std::cout<<std::endl;

	}
      for(int i=1; i<6; i++)
	{
	  if(dists[0]>dists[i])
	    {

	      //r2<r1 swap 
	      //swap s2,s3
	      if(i==1)
		{
		  clust_swap_atom(2,3,elements,dists);
		  swapped=true;
		}

	      //r3<r1
	      //swap s2,s4
	      if(i==2)
		{
		  clust_swap_atom(2,4,elements,dists);
		  swapped=true;
		}
	      //r4<r1
	      //swap s3,s1
	      if(i==3)
		{
		  clust_swap_atom(1,3,elements,dists);
		  swapped=true;
		}
	      //r5<r1
	      //swap s1,s4
	      if(i==4)
		{
		  clust_swap_atom(1,4,elements,dists);
		  swapped=true;
		}

	      //r6<r1
	      //swap s3,s1
	      //swap s4,s2
	      if(i==5)
		{
		  clust_swap_atom(1,3,elements,dists);
		  clust_swap_atom(2,4,elements,dists);
		  swapped=true;
		}
	    }
	}

      //rotations::

      


    }


  iterera=0;

  swapped=true;
	      
  std::vector<double> temp_distances;
  while(swapped)
    {
      swapped=false;

      for(int i=1; i<5; i++)
  	{
  	  for(int j=i+1; j<5; j++)
  	    {
	      // if(fabs(dists[i]-dists[j])<1e-4)
	      //	{
	      temp_distances=dists;
	      clust_swap_atom(i,j,elements,dists);
		  
	      for(int k=0; k<6; k++)
		{
		  if(fabs(temp_distances[k]-dists[k])<1e-4)
		    {
		      continue;
		    }

		  else if(temp_distances[k]<dists[k])
		    {
		      clust_swap_atom(i,j,elements,dists);
		      break;
		    }
		  else if(temp_distances[k]>dists[k])
		    {
		      swapped=true;
		      break;
		    }
		}
	      //	}
  	    }
  	}
      //rotate
      //-1 cw, -2=ccw

       if(!swapped)
	{
	  //	  std::cout<<"CW"<<std::endl;

	  //	  std::cout<<dists[0]<< " "<<dists[1]<< " "<<dists[2]<< " "<<dists[3]<< " "<<dists[4]<< " "<<dists[5]<< " "<<elements[0]<< " "<<elements[1]<< " "<<elements[2]<< " "<<elements[3]<<std::endl;

	  temp_distances=dists;

	  clust_swap_atom(-1,0,elements,dists);
	  // 	  std::cout<<dists[0]<< " "<<dists[1]<< " "<<dists[2]<< " "<<dists[3]<< " "<<dists[4]<< " "<<dists[5]<< " "<<elements[0]<< " "<<elements[1]<< " "<<elements[2]<< " "<<elements[3]<<std::endl;
	  for(int k=0; k<6; k++)
	    {
	      if(fabs(temp_distances[k]-dists[k])<1e-4)
		{
		  continue;
		}
	      else if(temp_distances[k]<dists[k])
		{
		  clust_swap_atom(-2,0,elements,dists);
		  break;
		}
	      else if(temp_distances[k]>dists[k])
		{
		  //  std::cout<<"good rotate"<<std::endl;
		  swapped=true;
		  break;
		}
	    }
	}
      if(!swapped)
	{
	  // std::cout<<"CCW"<<std::endl;
		  
	  // std::cout<<dists[0]<< " "<<dists[1]<< " "<<dists[2]<< " "<<dists[3]<< " "<<dists[4]<< " "<<dists[5]<< " "<<elements[0]<< " "<<elements[1]<< " "<<elements[2]<< " "<<elements[3]<<std::endl;

	  temp_distances=dists;
      
	  clust_swap_atom(-2,0,elements,dists);
	  //  std::cout<<dists[0]<< " "<<dists[1]<< " "<<dists[2]<< " "<<dists[3]<< " "<<dists[4]<< " "<<dists[5]<< " "<<elements[0]<< " "<<elements[1]<< " "<<elements[2]<< " "<<elements[3]<<std::endl;

	  for(int k=0; k<6; k++)
	    {
	      if(fabs(temp_distances[k]-dists[k])<1e-4)
		{
		  continue;
		}
	     else if(temp_distances[k]<dists[k])
		{
		  clust_swap_atom(-1,0,elements,dists);
		  break;
		}
	      else if(temp_distances[k]>dists[k])
		{
		  //  std::cout<<"good rotate"<<std::endl;

		  swapped=true;
		  break;
		}
	    }
	}

    }
	     
	      
	    

  iterera=0;
  swapped =false;
  while(swapped)
    {
      iterera++;
      if(iterera>10)
	{
	  std::cout<<"iterear: "<<iterera<<std::endl;
	  for(int k=0; k<dists.size(); k++)
	    {
	      std::cout<<dists[k]<< " ";
	    }
	  std::cout<<std::endl;

	}
      swapped=false;
      
      for(int i=0; i<2; i++)
	{
	  swapped=false;
	  if (dists[i]>dists[i+1])
	    {
	      if(i==0)
		{
		  //  std::cout<<"i=0"<<std::endl;
		  // r2 < r1
		  //swap s2,s3
		  //swap r1,r2 and r5,r6
		  clust_swap_atom(2,3,elements,dists);
		  swapped = true;
		}

	      if(i==1)
		{
		  // std::cout<<"i==1"<<std::endl;
		  //r3<r2
		  //swap s3,s4		   
		  // swap r2,r3
		  // swap r4,r5
		  clust_swap_atom(3,4,elements,dists);
		  swapped = true;
		}
	    }

	  // make sure r2 < r4 || r5
	  //if not swap s1,s2
	  //swapping s1,s2 means that
	  // r1=r1
	  //r2=r4
	  //r3=r5
	  //r6=r6
	  if(dists[1]> dists[3] && dists[2]>=dists[4])
	    {
	      //    std::cout<<"konstig swap h'r kanske"<<std::endl;
	      clust_swap_atom(1,2,elements,dists);
	      swapped=true;
	    }

	  
	}

      for(int i=1; i<5; i++)
	{
	  for(int j=i+1; j<5; j++)
	    {
	      if(fabs(dists[i]-dists[j])<1e-4)
		{
		  temp_distances=dists;
		  clust_swap_atom(i,j,elements,dists);
		  
		  for(int k=0; k<6; k++)
		    {
		      if(temp_distances[k]<dists[k])
			{
			  clust_swap_atom(i,j,elements,dists);
			  break;
			}
		      else if(temp_distances[k]>dists[k])
			{
			  swapped=true;
			  break;
			}
		    }
		}
	    }
	}


      //second part, sort alphabetically where possible
      // if r1==r2==r3 then sort r4,r5,r6
      double diff1 = fabs(dists[0]-dists[1]); // abs(r1-r2) //r1==r2
      double diff2 = fabs(dists[1]-dists[2]); // abs(r1-r3) //r2==r3
     
      double diff4 = fabs(dists[3]-dists[4]); // r4==r5
      double diff5 = fabs(dists[4]-dists[5]); //r5==6

      
      if(diff1<DISTANCE_LIMIT && diff2 < DISTANCE_LIMIT)
	{
	  // std::cout<<"in here?"<<std::endl;
	  swapped = true;
	  
	  while(swapped)
	    {
	      swapped=false;
	      for(int i=3; i<5; i++)
		{
		  if(dists[i]>dists[i+1])
		    {
		      if(i==3)
			{
			  //r5<r4
			  //swap r4,r5
			  //swap s3,s4
			  //r1,r2,r3 equal so no need to swap these
			  clust_swap_atom(3,4,elements,dists);
			  swapped = true;
			}

		      if(i==4)
			{
			  //r6<r5
			  //swap s2,s3
			  //swap r5,r6
			  clust_swap_atom(2,3,elements,dists);
			  swapped = true;
			}
		    }
		}
	    }
	  // diff4 = fabs(dists[3]-dists[4]);
	  // diff5 = fabs(dists[4]-dists[5]);
	  // if( diff4 < DISTANCE_LIMIT && diff5 < DISTANCE_LIMIT)
	  //   {
	  //     swapped = true;
	  //     while(swapped)
	  // 	{
	  // 	  swapped=false;
	  // 	  for(int i=1; i<3; i++)
	  // 	    {
	  // 	      if(reverseInt*(elements[i].compare(elements[i+1]))>0)
	  // 		{
	  // 		  temp_element=elements[i];
	  // 		  elements[i]=elements[i+1];
	  // 		  elements[i+1]=temp_element;
	  // 		  swapped=true;
	  // 		}
	  // 	    }
		  
	  // 	}
	      

	      
	}

      
    }

      
  std::vector<double> temp_dists;
  for(int i=0; i<dists.size(); i++)
    {
      temp_dists.push_back(dists[i]);
    }
  

  swapped=true;
  bool same_dists_after_swaps;
  while(swapped)
    {
      swapped=false;
      for(int i=0; i<3; i++)
  	{
	  for(int ii=i+1; ii<4; ii++)
	    {
	      // if(i==ii)
	      // 	{
	      // 	  continue;
	      // 	}
	      same_dists_after_swaps=false;
	      
	      if(reverseInt*elements[i].compare(elements[ii])>0)
		{

	      
		  clust_swap_atom(i+1,ii+1,elements,dists);
		  same_dists_after_swaps=true;
		  for(int j=0; j<temp_dists.size(); j++)
		    {
		      if(fabs(temp_dists[j]-dists[j])> 1e-4)
			{
			  same_dists_after_swaps=false;
		      
			}
		    }

		  if(same_dists_after_swaps)
		    {
		      swapped=true;
		    }
		  else
		    {
		      clust_swap_atom(i+1,ii+1,elements,dists);
		    }
		}

	    }
	  //rotation clockwise if 
	  // s2<=s1
	  // s4<=s2
	  //s3<=s4
	  //s1<=s3

	  //rotation c_clockwise is swap:
	  // desirable if the ones that go back are large and the ones that go forth are small

	  //s1 to s3
	  //s2 to s1
	  //s3 to s4
	  //s4 to s2
	  // old order, s1,s2,s3,s4
	  // new order, s2,s4,s1,s3
	  //first s2<=s1 s4<=s2, s3<=s1
	  //rotating ccw

	  bool rotate=false;
	  bool check_rotate=true;
	  while(check_rotate && swapped==false)
	    {

	      //first position
	      if(reverseInt*elements[1].compare(elements[0])!=0)
		{
		  if(reverseInt*elements[1].compare(elements[0])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      
	      //second position
	      else if(reverseInt*elements[3].compare(elements[1]) !=0)
		{
		  if(reverseInt*elements[3].compare(elements[1])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      //third position
	      else if(reverseInt*elements[0].compare(elements[2]) !=0)
		{
		  if(reverseInt*elements[0].compare(elements[2])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      //fourth position
	      else if(reverseInt*elements[2].compare(elements[3]) !=0)
		{
		  if(reverseInt*elements[2].compare(elements[3])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      check_rotate=false;
	    }
	  
	      //third postion
	      //if(

	      // else if(reverseInt*elements[1].compare(elements[0])==0)
	  //   {
	  //     if(reverseInt*elements[3].compare(elements[1])>0)
	  // 	{
	  // 	  rotate=true;
	  // 	}
	  //     else if(reverseInt*elements[3].compare(elements[1])==0)
	  // 	{
	  // 	  if(reverseInt*elements[0].compare(elements[2])<0)
	  // 	    {
	  // 	      rotate=true;
	  // 	    }
	  // 	  else if(reverseInt*elements[0].compare(elements[2])==0)
	  // 	    {
	  // 	      if(reverseInt*elements[2].compare(elements[3])<0)
	  // 		{
	  // 		  rotate=true;
	  // 		}
	  // 	    }
	  // 	}
	  //   }
	    



	  //changing back is the reverse
	     //   if(reverseInt*elements[0].compare(elements[3])>= 0 &&
	     // reverseInt*elements[1].compare(elements[0])<= 0 &&
	     // reverseInt*elements[2].compare(elements[1])<= 0 &&
	     // reverseInt*elements[3].compare(elements[2])<= 0 &&
	     // (reverseInt*elements[0].compare(elements[3])> 0 ||
	     //  reverseInt*elements[1].compare(elements[0])< 0 ||
	     //  reverseInt*elements[2].compare(elements[1])< 0 ||
	     //  reverseInt*elements[3].compare(elements[2])< 0 )

	  if(rotate)	    
	    {
	      //-1 cw, -2=ccw
	      clust_swap_atom(-1,0,elements,dists);
	      same_dists_after_swaps=true;

	      for(int j=0; j<temp_dists.size(); j++)
	  	{
	  	  if(fabs(temp_dists[j]-dists[j])> 1e-4)
	  	    {
	  	      same_dists_after_swaps=false;		      
	  	    }
	  	}

	      if(same_dists_after_swaps)
	  	{
	  	  swapped=true;
	  	}
	      else
	  	{
	  	  clust_swap_atom(-2,0,elements,dists);
	  	}
	    }
	

	  //clockwise rotation
	  // old: 1,2,3,4
	  // new: 3,1,4,2
	  rotate=false;
  while(check_rotate && swapped==false)
	    {

	      //first position
	      if(reverseInt*elements[2].compare(elements[0])!=0)
		{
		  if(reverseInt*elements[2].compare(elements[0])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      
	      //second position
	      else if(reverseInt*elements[0].compare(elements[1]) !=0)
		{
		  if(reverseInt*elements[0].compare(elements[1])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      //third position
	      else if(reverseInt*elements[3].compare(elements[2]) !=0)
		{
		  if(reverseInt*elements[3].compare(elements[2])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      //fourth position
	      else if(reverseInt*elements[1].compare(elements[3]) !=0)
		{
		  if(reverseInt*elements[1].compare(elements[3])<0)
		    {
		      rotate=true;
		      check_rotate=false;
		      break;
		    }
		  else
		    {
		      rotate=false;
		      check_rotate=false;
		      break;
		    }
		}
	      check_rotate=false;
	    }


  	  if(rotate)	    
	    {
	      //-1 cw, -2=ccw
	      clust_swap_atom(-2,0,elements,dists);
	      same_dists_after_swaps=true;

	      for(int j=0; j<temp_dists.size(); j++)
	  	{
	  	  if(fabs(temp_dists[j]-dists[j])> 1e-4)
	  	    {
	  	      same_dists_after_swaps=false;		      
	  	    }
	  	}

	      if(same_dists_after_swaps)
	  	{
	  	  swapped=true;
	  	}
	      else
	  	{
	  	  clust_swap_atom(-1,0,elements,dists);
	  	}
	    }

	}
    }


  //r2==r3, if swap s3,s4 it will also swap r4,r5 so only swap it if that doesn't change the distance order
  //   if(diff2 < DISTANCE_LIMIT && diff4 < DISTANCE_LIMIT)
  // 	{
  // 	  //alphabetically sort s3,s4
  // 	  if(reverseInt*(elements[2].compare(elements[3]))>0)
  // 	    {
  // 	      temp_element=elements[2];
  // 	      elements[2]=elements[3];
  // 	      elements[3]=temp_element;
  // 	    }
  // 	}
  //   // r1==r3, swaps s2,s3 and r5,r6 but only if r5==r6
  //   if(diff1 < DISTANCE_LIMIT && diff5  < DISTANCE_LIMIT)
  // 	{
  // 	  // alphabetically sort s2,s3
  // 	  if(reverseInt*(elements[1].compare(elements[2]))>0)
  // 	    {
  // 	      temp_element=elements[2];
  // 	      elements[2]=elements[3];
  // 	      elements[3]=temp_element;
  // 	    }    
  // 	}
  // }
      
  */
		      
}

  
void clust_swap_atom(int atom1, int atom2,std::vector<std::string> &element, std::vector<double> &dists)
{
  std::string temp_element;
  double tempDistance;
  


  


  if(atom2<atom1)
    {
      int temp_atom = atom1;
      atom1=atom2;
      atom2=temp_atom;
    }
  
  if(atom1==1)
    {

      //swap 1,2
      //swap r2,r4
      //swap r3,r5
      if(atom2==2)
	{

	  clust_swap_dist(1,3,dists);
	      
	  //r3==r5
	  clust_swap_dist(2,4,dists);
	      
	  //swap s1,s2
	  clust_swap_element(0,1,element);
	}


      //swap s1,s3
      // swap r1-r4
      // r2 same
      //r3-r6
     
      else if(atom2==3)
	{
	  clust_swap_dist(0,3,dists);      
	  //r3==r6
	  clust_swap_dist(2,5,dists);
	      
	  //swap s1,s3
	  clust_swap_element(0,2,element); 
	}

      //swap s1,s4
      //swap r1-r5
      //swap r2-r6
      else if(atom2==4)
	{
	  clust_swap_dist(0,4,dists);
	  clust_swap_dist(1,5,dists);
	  clust_swap_element(0,3,element);
	}
    }//atom1==1

  else if(atom1==2)
    {

      //swap s2,s3
      //swap r1,r2
      //swap r5,r6
      if(atom2==3)
	{
	  clust_swap_dist(0,1,dists);
	  clust_swap_dist(4,5,dists);
	  clust_swap_element(1,2,element);
	}
      //swap s2,s4
      //swap r1,r3
      //swap r4,r6
      else if(atom2==4)
	{
	  clust_swap_dist(0,2,dists);
	  clust_swap_dist(3,5,dists);
	  clust_swap_element(1,3,element);
	}

    }

  else if(atom1==3)
    {

      //swap s3,s4
      //swap r2,r3
      //swap r4,r5
      if(atom2==4)
	{

	  clust_swap_dist(1,2,dists);
	  clust_swap_dist(3,4,dists);
	  clust_swap_element(2,3,element);
	}
    }
  //clockwise rotation
	  //s1,s2
	  //s2,s4
	  //s3,s4
  else if(atom1==-1)
    {
      clust_swap_atom(1,2,element,dists);
      clust_swap_atom(2,4,element,dists);
      clust_swap_atom(3,4,element,dists);
    }
  //counter-clockwise rotation
  else if(atom1==-2)
    {
      clust_swap_atom(3,4,element,dists);
      clust_swap_atom(2,4,element,dists);
      clust_swap_atom(1,2,element,dists);

    }
      

}

void clust_swap_dist(int index1,int index2,std::vector<double> &dists)
{
  double temp_distance =dists[index1];

  dists[index1]=dists[index2];
  dists[index2]=temp_distance;
}

  

void clust_swap_element(int index1, int index2, std::vector<std::string> &elements)
{
  std::string temp_element=elements[index1];
  elements[index1]=elements[index2];
  elements[index2]=temp_element;
}
  
  







void clust_sort_triplet(std::vector<double> &dists, std::vector<std::string> &elements,bool reverseSort)
{
  bool swapped=true;
  double tempDr;
  std::string temp_element;
  

  // std::cout<<"before"<<std::endl;
  // std::cout<<dists[0]<< " "<< dists[1]<< " "<<dists[2]<< " "<< elements[0]<< " "<<elements[1]<< " "<<elements[2]<< " "<<reverseSort<<std::endl;
  //first step gives the most reduced form of the distances
  //second step consists of finding the alphabetical order of elements which are symmetric equivalent
  while(swapped)
    {
      swapped=false;
      for(int i=0; i<2; i++)
	{
	  if(dists[i]>dists[i+1])
	    {
	      if(i==0)
		{
		  //swap atoms 2 and 3
		  temp_element=elements[1];
		  elements[1]=elements[2];
		  elements[2]=temp_element;
		}
	      if(i==1)
		{
		  //swap atoms 1 and 2
		  temp_element=elements[0];
		  elements[0]=elements[1];
		  elements[1]=temp_element;
		}
	      
	      tempDr=dists[i];
	      dists[i]=dists[i+1];
	      dists[i+1]=tempDr;	      
	      swapped=true;
	    }
	}
    }
      //second part
      double diff1=fabs(dists[0]-dists[1]);
      double diff2=fabs(dists[1]-dists[2]);
      std::vector<std::string> types;
      double DISTANCE_LIMIT=1e-4;
      //all distances are equal
      if(diff1 < DISTANCE_LIMIT && diff2 < DISTANCE_LIMIT)
	{
	  //sort all elements
	  sortElements(elements,reverseSort);
	}
      else if(diff2 < DISTANCE_LIMIT)
	{
	  //distance2 == distance 3
	  // site 1 and site 2 equaivalent
	  
	  types.push_back(elements[0]);
	  types.push_back(elements[1]);
	  sortElements(types,reverseSort);
	  elements[0]=types[0];
	  elements[1]=types[1];
	}
      else if(diff1<DISTANCE_LIMIT)
	{
	  //d1==d2
	  //sort site2 site3
	  types.push_back(elements[1]);
	  types.push_back(elements[2]);
	  sortElements(types,reverseSort);
	  elements[1]=types[0];
	  elements[2]=types[1];
	}

     
	//  std::cout<<dists[0]<< " "<< dists[1]<< " "<<dists[2]<< " "<< elements[0]<< " "<<elements[1]<< " "<<elements[2]<<std::endl;

}
	


void sortElements (std::vector<std::string> &elements,bool reverseSort)
{

  double reverse_int=1;
  if(reverseSort)
    {
      reverse_int=-1;
    }




  bool swapped = true;
  std::string tempString;
  
  while(swapped)
    {
      swapped= false;
      for(int i=0; i<elements.size()-1; i++)
	{
	  if(reverse_int*elements[i].compare(elements[i+1])>0)
	    {
	      tempString = elements[i];
	      elements[i]=elements[i+1];
	      elements[i+1]=tempString;
	      swapped=true;
	    }
	}
    }
}




  
