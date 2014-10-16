#include "clust.hpp"
#include <math.h>
#include <iostream>
#include <string>
#define PI 3.14159265

double clusterFunction(int nbr_of_elements,int atom_type,int clustFunction)
{
   //odd: use cos
  if((clustFunction%2!=0))
    {
      return -cos(2.0*PI * (double)( (int)(clustFunction+2)/2 ) * (double)atom_type/((double)nbr_of_elements));
    }
  //even: use sin
  else 
    {
      return -sin(2.0*PI * (double)( (int)(clustFunction+2)/2 ) * (double)atom_type/((double)nbr_of_elements));
    }
}



std::vector<std::vector<std::string> > symmetric_cluster_function(std::vector<double> dists,std::vector<std::string> subelements)
{
  std::vector<std::vector<int> > int_base = symmetric_cluster_function(dists,subelements.size()+1);
  std::vector<std::vector<std::string> > ret;
  ret.resize(int_base.size());
  for(int i=0; i<int_base.size(); i++)
    {
      for(int j=0; j<int_base[i].size(); j++)
	{
	  ret[i].push_back( subelements[ int_base[i][j] ] );
	}
    }
  return ret;
}
  


std::vector<std::vector<int> > symmetric_cluster_function(std::vector<double> dists, int nbr_of_elements)
{
  int tuple_size=dists.size();
  std::vector<std::vector<int> >returnVector;
  
  int cluster_max=nbr_of_elements-2;
  
  if(tuple_size==1)
    {
      std::vector<int> one_ret;
      one_ret.resize(1);
      for(int i=0; i<cluster_max+1; i++)
	{
	  one_ret[0]=i;
	  returnVector.push_back(one_ret);
	}
    }  
  else if(tuple_size==2)
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
	      
	      sortElements(temp_element,true);
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
		  
		  tuple_remodulator(temp_dists,temp_element,true);

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

		      tuple_remodulator(temp_dists,temp_element,true);

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


  return returnVector;
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


  if(dists.size()<4 && dists.size() != elements.size())
    {
      std::cout<<"distance/elements mismatch in tuple_remodulator"<<std::endl;
      return;
    }
  int d_size = dists.size();
  
  if(d_size==2)
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



void clust_sort_quatuplet(std::vector<double> &dists, std::vector<std::string> &elements,bool reverseSort)
{

  double reverseInt =1.0;
  if(reverseSort)
    {
      reverseInt=-1.0;
    }
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
  double DISTANCE_LIMIT=1e-4;

  bool swapped = true;
  double temp_distance;
  std::string temp_element;


  while(swapped)
    {
      swapped=false;
      
      for(int i=0; i<2; i++)
	{
	  swapped=false;
	  if (dists[i]>dists[i+1])
	    {
	      if(i==0)
		{
		  // r2 < r1
		  //swap s2,s3
		  //swap r1,r2 and r5,r6
		  clust_swap_atom(2,3,elements,dists);

		  swapped = true;

		}

	      if(i==1)
		{
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
	  if(dists[1]> dists[3] || dists[1]>dists[4])
	    {
	      clust_swap_atom(1,2,elements,dists);
	      swapped=true;
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
  	  same_dists_after_swaps=false;
	      
  	  if(reverseInt*elements[i].compare(elements[i+1])>0)
  	    {

	      
  	      clust_swap_atom(i+1,i+2,elements,dists);
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
		  clust_swap_atom(i+1,i+2,elements,dists);
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




  
