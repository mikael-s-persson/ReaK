
#include <iostream>

#include "mat_alg.hpp"



using namespace ReaK;

int main() {
  mat<double,mat_structure::rectangular> M1(4,4);
  
  M1(0,0) =  1; M1(0,1) =  2; M1(0,2) =  3; M1(0,3) =  4;
  M1(1,0) =  5; M1(1,1) =  6; M1(1,2) =  7; M1(1,3) =  8;
  M1(2,0) =  9; M1(2,1) = 10; M1(2,2) = 11; M1(2,3) = 12;
  M1(3,0) = 13; M1(3,1) = 14; M1(3,2) = 15; M1(3,3) = 16;
  
  std::cout << "M1 = " << M1 << std::endl << std::endl;
  
  auto M2 = ( ( sub(M1)(range(0,2),range(0,1)) & sub(M1)(range(0,2),range(2,3)) ) |
              ( sub(M1)(range(3,3),range(0,1)) & sub(M1)(range(3,3),range(2,3)) ) );
  
  std::cout << "M2 = " << M2 << std::endl << std::endl;
  
  return 0;
};













