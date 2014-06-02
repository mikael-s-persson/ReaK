
#include <iostream>

#include <ReaK/core/lin_alg/mat_alg.hpp>



using namespace ReaK;

int main() {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  
  mat<double,mat_structure::rectangular> M1(4,4);
  
  M1(0,0) =  1; M1(0,1) =  2; M1(0,2) =  3; M1(0,3) =  4;
  M1(1,0) =  5; M1(1,1) =  6; M1(1,2) =  7; M1(1,3) =  8;
  M1(2,0) =  9; M1(2,1) = 10; M1(2,2) = 11; M1(2,3) = 12;
  M1(3,0) = 13; M1(3,1) = 14; M1(3,2) = 15; M1(3,3) = 16;
  
  std::cout << "M1 = " << M1 << std::endl << std::endl;
  
  auto M2 = ( ( sub(M1)(range(0,3),range(0,2)) & sub(M1)(range(0,3),range(2,4)) ) |
              ( sub(M1)(range(3,4),range(0,2)) & sub(M1)(range(3,4),range(2,4)) ) );
  
  std::cout << "M2 = " << M2 << std::endl << std::endl;
  
  std::cout << "M2 0; 0 1; = " << ( ( M2 & mat_nil<double>(4,2) ) |
                                    ( mat_nil<double>(2,4) & mat_ident<double>(2) ) ) << std::endl;
  
  vect_n<double> V1(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0);
  
  std::cout << "V1 = " << V1 << std::endl;
  
  std::cout << "V1 as 4x3 col-major = " << make_mat(V1)(range(0,4),3) << std::endl;
  std::cout << "V1 as 4x3 row-major = " << make_mat(V1)(4,range(0,3)) << std::endl;
  
  
#endif
  return 0;
};













