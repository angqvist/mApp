#include <iostream>
#include <fstream>

#include <string>



int main()
{
  std::string fileName;
  fileName="ideal_pos";
  std::cout<<fileName<<std::endl;
  
  std::ifstream in(fileName.c_str());
  

  in.close(); 
  return 0;
}
