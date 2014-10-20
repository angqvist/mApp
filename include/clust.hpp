#include <vector>
#include <string>


double clusterFunction(int, int, int);
std::vector<std::vector<int> > symmetric_cluster_function(std::vector<double> , int,bool);
std::vector<std::vector<std::string> > symmetric_cluster_function(std::vector<double> ,std::vector<std::string> );
void clust_sort_return_vector(std::vector<std::vector<int> > &);
bool add_new_cluster_controller(std::vector<std::vector<int> >, std::vector<int> );

void tuple_remodulator(std::vector<double> &, std::vector<std::string> &,bool);
void clust_sort_quatuplet(std::vector<double> &, std::vector<std::string> &,bool);

void clust_sort_triplet(std::vector<double> &,std::vector<std::string> &,bool);
void sortElements(std::vector<std::string> &,bool);
void clust_swap_atom(int,int,std::vector<std::string> &, std::vector<double> &);

void clust_swap_dist(int,int,std::vector<double> &);
void clust_swap_element(int,int,std::vector<std::string> &);
