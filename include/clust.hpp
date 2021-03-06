#pragma once
#include <vector>
#include <string>


double clusterFunction(int &, int &,  int &);
std::vector<std::vector<int> > symmetric_cluster_function(std::vector<double> & , int &,bool &);
std::vector<std::vector<std::string> > symmetric_cluster_function(std::vector<double> &,std::vector<std::string> &);
void clust_sort_return_vector(std::vector<std::vector<int> > &);
void clust_sort_dists(std::vector<std::vector<double> > &);
bool add_new_cluster_controller(std::vector<std::vector<int> > &, std::vector<int> &);
int is_first_dist_smaller(std::vector<double> &, std::vector<double> &);
int is_elements_in_lower_order(std::vector<std::string> &, std::vector<std::string> &,bool & );

int trial_swap(std::vector<int> &, std::vector<double> &, std::vector<double> &, std::vector<std::string> &, std::vector<std::string> &,bool &);

void tuple_remodulator(std::vector<double> &, std::vector<std::string> &,bool &);
void clust_sort_quatuplet_part2(std::vector<double> &, std::vector<std::string> &,bool &);
void clust_sort_quatuplet(std::vector<double> &, std::vector<std::string> &,bool &);

void clust_sort_triplet(std::vector<double> &,std::vector<std::string> &,bool &);
void sortElements(std::vector<std::string> &,bool &);
void clust_swap_atom(int ,int ,std::vector<std::string> &, std::vector<double> &);

void clust_swap_dist(int ,int ,std::vector<double> &);
void clust_swap_element(int ,int ,std::vector<std::string> &);
