#include "NeighbourList.hpp"
#include "Neighbour.hpp"
#include "LatticeList.hpp"
#include "ParameterList.hpp"
#include "Pair.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include "clust.hpp"
#include "Triplet.hpp"
#include "Quatuplet.hpp"
// should really not be used with empty constructor
NeighbourList::NeighbourList()
{
  thisIndex=0;
  distanceLimit=1e-4;
  offset=0;
  // empty
}
//check up if you should send references to latticeList or something like that
NeighbourList::NeighbourList(int i, LatticeList ll, ParameterList pl)
{
  thisIndex=i;
  distanceLimit=1e-4;
  findNeighbours(ll,pl);
  findTripletNeighbours(ll,pl);
  findQuatupletNeighbours(ll,pl);
  setOffset(pl.getOffsetValue());
}
void NeighbourList::addNeighbour(Neighbour nbr1)
{
  if(newNeighbour(nbr1))
    {
      nbrList.push_back(nbr1);
    }
}
int NeighbourList::newNeighbour(Neighbour nbr1)
{
  for(size_t i=0; i < nbrList.size(); i++)
    {
      if(nbr1==nbrList[i])
      {
	// std::cout<<nbrList.size()<<std::endl;
	// std::cout<<"Tried to add neighbour which already existed in neighbourList."<<std::endl<<"Neighbour tried to add: ";
	// nbr1.print();
	// std::cout<<"Index: "<<i<<std::endl;
	return false;
      }
    }
  return true;
}


void NeighbourList::printList()
{
  for(size_t i=0; i < manyNbrLists.size(); i++)
    {
      std::cout<<indexList[i]<<": ";
      printManyNbr(manyNbrLists[i]);
      std::cout<<"--------------------------"<<std::endl;
      std::cout<<" "<<std::endl;
    }
}

void NeighbourList::printManyNbr(std::vector<Neighbour> nl)
{
  for(size_t i=0; i < nl.size(); i++)
    {
      nl[i].print();
    }
}

//finds pairs
void NeighbourList::findNeighbours(LatticeList ll, ParameterList pl)
{
  double tempDist;
  bool paramAtIndex;
  int singletIndex=0;

  

  //commented out old singlet stuff. new is the new old.
  // while(pl.getPair(singletIndex).getDistance()<distanceLimit)
  //   {
  //     //std::cout<<"index: "<<singletIndex<<" dist: "<<pl.getPair(singletIndex).getDistance()<< " type "<<pl.getPair(singletIndex).getSite1()<< " energy "<<pl.getPair(singletIndex).getEnergy()<<std::endl;
  //     singletType.push_back(pl.getPair(singletIndex).getSite1());
  //     singletEnergy.push_back(pl.getPair(singletIndex).getEnergy());
  //     singletIndex++;
  //   }
  //std::cout<<"singlet index="<<singletIndex<<std::endl;
 
  ll.calculate_lookup_table();
  
  std::vector<double> temp_dist;
  temp_dist.resize(1);
  std::vector<std::string> ghost_atoms;
  ghost_atoms.resize(2);
  ghost_atoms[0]="A";
  ghost_atoms[1]="A";

  

  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      paramAtIndex=false;
      if(i == thisIndex)
	{
	  continue;
	}
      // tempDist=ll.getDistance(i,thisIndex);
      tempDist[0]=ll.fast_distance(i,thisIndex);
      
      tuple_remodulator(temp_dist,ghost_atoms,false);      
      std::vector<Pair> temp_pair_vector;
      
      for(int j=0; j<pl.getNbrOfPairs(); j++)
	{
	  if(fabs(tempDist-pl.getPair(j).getDistance())<distanceLimit)
	    {
	      paramAtIndex=true;
	      temp_pair_vector.push_back(pl.getPair(j));
	    }
	}

      if(paramAtIndex)
	{
	  pair_vector.push_back(temp_pair_vector);
	  pair_index.push_back(i);
	  	  
	}
	      

      // Neighbour tempNbr = Neighbour(pl.getPair(j).getSite1(),pl.getPair(j).getSite2(),pl.getPair(j).getEnergy(),pl.getPair(j).getDistance());
      // if( newNeighbour(tempNbr) )
      // 	{
      // 	  nbrList.push_back(tempNbr);
      // 	  paramAtIndex=true;		      
      // 	}


    }
}
	
//   if(!(nbrList.empty()))
// 	{
// 	  manyNbrLists.push_back(nbrList);
    // 	}
    //   nbrList.clear();
    //   if(paramAtIndex)
    // 	{
    // 	  indexList.push_back(i);
    // 	}
    // }// end i loop



  








double  NeighbourList::getMatchingEnergy(Neighbour nbr, std::vector<Neighbour> nl)
{
  for(size_t j=0; j<nl.size(); j++)
    {
      //std::cout<<"tempNbr: "<<std::endl;
      //nbr.print();
      //std::cout<<"neigbourlist  "<<std::endl;
      //nl[j].print();

      if(nbr==nl[j])
	{
	  return nl[j].getEnergy();
	}
    }
  return 0;
}

double NeighbourList::getLocalEnergy(LatticeList ll)
{

  ll.calculate_lookup_table();

  localEnergy=0.0;
  //add singlet Energy times two?
  std::cout<<"Getting local energy"<<std::endl;
  std::cout<<"pair index size: "<<indexList.size()<<std::endl;
  std::cout<<"Triplet index size "<<trip_index.size()<<std::endl;
  std::cout<<"Quatuplet index size "<<quat_index.size()<<std::endl;
  
  // for(int i=0; i<singletEnergy.size(); i++)
  //   {
  //     if(ll.getSite(thisIndex)==singletType[i])
  // 	{
  // 	  localEnergy+=singletEnergy[i];
  // 	}
  //   }
  // Neighbour tempNbr= Neighbour();
  // tempNbr.setSite1(ll.getSite(thisIndex));  
  // //TIMES A HALF since all pairs will get counted twice
  // for(size_t i=0; i < indexList.size(); i++)
  //   {
  //     tempNbr.setSite2(ll.getSite(indexList[i]));
  //     std::cout<<"found pair with energy "<< getMatchingEnergy(tempNbr,manyNbrLists[i])*0.5<< " total energy: "<<localEnergy<<std::endl;
  //     localEnergy += getMatchingEnergy(tempNbr,manyNbrLists[i])*0.5;
  //   }

  for(int i=0; i<singletList.size(); i++)
    {
      if(ll.getSite(thisIndex)==singletList[i].getType())
	{
	  localEnergy += singletList[i].getProperty();
	  break;
	}
    }
  std::vector<double> dists;
  std::vector<std::string> sites;

  
  if(pair_vector.size()>0)
    {
      dists.resize(1);
      sites.resize(2);
      Pair temp_pair = Pair();
      for(size_t i=0; i<pair_index.size(); i++)
	{
	  dists[0]=ll.fast_distance(thisIndex,pair_index[i]);
	  sites[0]=ll.getSite(i);
	  sites[1]=ll.getSite(thisIndex);
	  temp_pair.setDistance(dists[0]);
	  temp_pair.setSite1(sites[0]);
	  temp_pair.setSite2(sites[1]);
	  
	  for(size_t j=0; j<pair_vector[i].size(); j++)
	    {
	      localEnergy ++ pair_vector[i][j].getEnergy()*0.5;
	      break;
	    }
	}
    }
	  
  

  //triplets;;; times a third since each triplet gets counted thrice

  dists.resize(3);
  sites.resize(3);
  Triplet temp_trip=Triplet();
  for(size_t i=0; i<trip_index.size(); i++)
    {
      dists[0]=ll.fast_distance(thisIndex,trip_index[i][0]);
      dists[1]=ll.fast_distance(thisIndex,trip_index[i][1]);
      dists[2]=ll.fast_distance(trip_index[i][0],trip_index[i][1]);
      sites[0]=ll.getSite(thisIndex);
      sites[1]=ll.getSite(trip_index[i][0]);     
      sites[2]=ll.getSite(trip_index[i][1]);
      tuple_remodulator(dists,sites,false);
      temp_trip.setAll(dists,sites);

      for(size_t j=0; j<trip_vector[i].size(); j++)
	{
	  if(temp_trip==trip_vector[i][j])
	    {
	      localEnergy +=trip_vector[i][j].getEnergy()*0.333333333333;
	      // std::cout<<"found triplet with energy: "<<trip_vector[i][j].getEnergy()*0.333333333333<< " "<<"current total energy "<<localEnergy<<std::endl;
	      break;
	    }
	}
    }



  //quatuplets times a quarters since.... bla bla bla
  if(quat_index.size()>0)
    {
      dists.resize(6);
      sites.resize(4);
      Quatuplet temp_quat = Quatuplet();
      for(size_t i=0; i<quat_index.size(); i++)
	{
	  
	  dists[0]=ll.fast_distance(thisIndex,quat_index[i][0]);
	  dists[1]=ll.fast_distance(thisIndex,quat_index[i][1]);
	  dists[2]=ll.fast_distance(thisIndex,quat_index[i][2]);
	  dists[3]=ll.fast_distance(quat_index[i][0],quat_index[i][1]);
	  dists[4]=ll.fast_distance(quat_index[i][0],quat_index[i][2]);
	  dists[5]=ll.fast_distance(quat_index[i][1],quat_index[i][2]);
	  
	  sites[0]=ll.getSite(thisIndex);
	  sites[1]=ll.getSite(quat_index[i][0]);
      	  sites[2]=ll.getSite(quat_index[i][1]);
	  sites[3]=ll.getSite(quat_index[i][2]);
	  
	  tuple_remodulator(dists,sites,false);
	  temp_quat.setAll(dists,sites);
	  for(size_t j=0; j<quat_vector[i].size(); j++)
	    {
	      if(temp_quat==quat_vector[i][j])
		{
		  localEnergy += quat_vector[i][j].getEnergy()*0.25;
		  break;
		}
	    }
	}
    }
      




  
  
  return localEnergy;
}




int NeighbourList::isThisMatchingNeighbour(LatticeList ll, Neighbour n1)
{
  std::cout<<"ERROR ERROR ERROR ERROR"<<std::endl;
  std::cout<<" OMG IS THIS USED?!?!??! Contact Mattias"<<std::endl;
  std::cout<<"ERROR ERROR ERROR ERROR"<<std::endl;
  Neighbour tempNbr= Neighbour();
  tempNbr.setSite1(ll.getSite(thisIndex));
  int returnNumber=0;
  for(size_t i=0; i<indexList.size(); i++)
    {
      tempNbr.setDistance(ll.getDistance(thisIndex,indexList[i]));
      tempNbr.setSite2(ll.getSite(indexList[i]));			  
      if(tempNbr==n1 && fabs(tempNbr.getDistance()-n1.getDistance())<1e-3)
	{
	  // std::cout<<"==================== "<< thisIndex<< " "<<indexList[i]<<std::endl;
	  // tempNbr.print();
	  // n1.print();
	  returnNumber++;
	}
    }
  return returnNumber;
}
      
      
double NeighbourList::getCurrentLocalEnergy()
{
  return currentLocalEnergy;
}

void NeighbourList::setCurrentLocalEnergy(double newEnergy)
{
 currentLocalEnergy=newEnergy;
}

void NeighbourList::calcCurrentLocalEnergy(LatticeList ll)
{
  currentLocalEnergy = getLocalEnergy(ll);
}


/*
  loop over two atoms i,j. get the distances thisIndex,i,j
  then you need to find equivalent distances in parameterlist
  if you do then those will be possible clusters for lattice index: {thisIndex, i,j}
  

 */

void NeighbourList::findTripletNeighbours(LatticeList ll,ParameterList pl)
{
  ll.calculate_lookup_table();
  std::vector<double> temp_dist;
  temp_dist.resize(3);
  std::vector<std::string> ghost_atoms;
  ghost_atoms.resize(3);
  ghost_atoms[0]="A";
  ghost_atoms[1]="A";
  ghost_atoms[2]="A";
  // std::cout<<"3: ll.getNbrOfSites() "<<ll.getNbrOfSites()<<std::endl;
  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      if(i==thisIndex)
	{
	  continue;
	}
      for(int j=i+1; j<ll.getNbrOfSites(); j++)
	{
	  if(j==thisIndex)
	    {
	      continue;
	    }
	  std::vector<Triplet> tempTripVector;
	  bool somethingAtThisIndex=false;
	  temp_dist[0]=ll.fast_distance(thisIndex,i);
	  temp_dist[1]=ll.fast_distance(thisIndex,j);
	  temp_dist[2]=ll.fast_distance(i,j);
	  
	  tuple_remodulator(temp_dist,ghost_atoms,false);
	  //  std::cout<<" pl.getNbrOFTriplets() "<<pl.getNbrOfTriplets()<<std::endl;
	  //	  std::cout<< temp_dist[0]<< " "<<temp_dist[1]<< " "<<temp_dist[2]<<std::endl;
	  for(int k=0; k<pl.getNbrOfTriplets(); k++)
	    {
	      // pl.getTriplet(k).printTriplet();
	      //std::cout<< pl.getTriplet(k).getDistance1()-temp_dist[0]<< " "<< pl.getTriplet(k).getDistance2()-temp_dist[1]<< " "<<pl.getTriplet(k).getDistance3()-temp_dist[2]<<std::endl;
	      if(fabs(pl.getTriplet(k).getDistance1()-temp_dist[0])<1e-4 &&
		 fabs(pl.getTriplet(k).getDistance2()-temp_dist[1])<1e-4 &&
		 fabs(pl.getTriplet(k).getDistance3()-temp_dist[2])<1e-4)
		{
		  //  std::cout<<"Found something triplety in find triplets"<<std::endl;
		  somethingAtThisIndex=true;
		  tempTripVector.push_back(pl.getTriplet(k));		  
		}
	    } 
	  if(somethingAtThisIndex)
	    {
	      std::vector<int> temp_index;
	      temp_index.push_back(i);
	      temp_index.push_back(j);
	      trip_index.push_back(temp_index);
	      trip_vector.push_back(tempTripVector);
	    }
	  
	  
	}
    }
}

void NeighbourList::findQuatupletNeighbours(LatticeList ll, ParameterList pl)
{
  ll.calculate_lookup_table();
  std::vector<double> temp_dists;
  temp_dists.resize(6);
  std::vector<std::string> ghost_atom;
  for(int i=0; i<6; i++)
    {
      ghost_atom.push_back("A");
    }

  for(int i=0; i<ll.getNbrOfSites(); i++)
    {
      if(i==thisIndex)
	{
	  continue;
	}
      for(int j=i+1; j<ll.getNbrOfSites(); j++)
	{
	  if(j==thisIndex)
	    {
	      continue;
	    }

	  for(int k=j+1; k<ll.getNbrOfSites(); k++)
	    {
	      if(k==thisIndex)
		{
		  continue;
		}
	      bool somethingAtThisIndex=false;
	      std::vector<Quatuplet> temp_quat_vec;
	      temp_dists[0]=ll.fast_distance(thisIndex,i);
	      temp_dists[1]=ll.fast_distance(thisIndex,j);
	      temp_dists[2]=ll.fast_distance(thisIndex,k);
	      temp_dists[3]=ll.fast_distance(i,j);
	      temp_dists[4]=ll.fast_distance(i,k);
	      temp_dists[5]=ll.fast_distance(j,k);

	      tuple_remodulator(temp_dists,ghost_atom,false);
	      for(int l=0; l<pl.getNbrOfQuatuplets(); l++)
		{
		  if(fabs(pl.getQuatuplet(l).getDistance(0)-temp_dists[0])<1e-4 &&
		     fabs(pl.getQuatuplet(l).getDistance(1)-temp_dists[1])<1e-4 &&
		     fabs(pl.getQuatuplet(l).getDistance(2)-temp_dists[2])<1e-4 &&
		     fabs(pl.getQuatuplet(l).getDistance(3)-temp_dists[3])<1e-4 &&
		     fabs(pl.getQuatuplet(l).getDistance(4)-temp_dists[4])<1e-4 &&
		     fabs(pl.getQuatuplet(l).getDistance(5)-temp_dists[5])<1e-4 )
		    {
		      somethingAtThisIndex=true;
		      temp_quat_vec.push_back(pl.getQuatuplet(l));      
		    }
	      
		  if(somethingAtThisIndex)
		    {
		      std::vector<int> temp_index;
		      temp_index.push_back(i);
		      temp_index.push_back(j);
		      temp_index.push_back(k);
		      quat_index.push_back(temp_index);
		      quat_vector.push_back(temp_quat_vec);
		    }
		}
	    }
	}
    }

  
}



void NeighbourList::setOffset(double newOff)
{
  offset_value=newOff;
}

double NeighbourList::getOffset()
{
  return offset_value;
}

void NeighbourList::readSingletList(ParameterList pl)
{
  singletList = pl.returnSingletVector();
}
