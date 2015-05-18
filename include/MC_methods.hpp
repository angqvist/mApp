#pragma once
#include <vector>
#include <string>



double MC_totalEnergy(class LatticeList &,  std::vector<class NeighbourList> &);

double MC_orderValue(class LatticeList, double);

void MC_occupations(class LatticeList, double &, double &, double &,std::string);

//int MC_neighboursInInterval(class LatticeList,double,double);


std::vector<int> neighbourCount(class LatticeList,std::vector<std::string >,double,double);
// void neighbourCount(class LatticeList,std::vector<std::string >,double,double,std::string);
