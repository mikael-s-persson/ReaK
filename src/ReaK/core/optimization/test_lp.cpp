
/*
 *    Copyright 2011 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).  
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#include "serialization/xml_archiver.hpp"

#include "simplex_method.hpp"
#include "mehrotra_method.hpp"

#include <iostream>
#include <cmath>



using namespace ReaK;

void lp01(vect_n<double>& c, mat<double,mat_structure::rectangular>& A, 
	  vect_n<double>& b, vect_n<double>& l, vect_n<double>& u) {
  u = vect_n<double>(22,std::numeric_limits<double>::infinity());
  l = vect_n<double>(22,0);
  
  b = vect_n<double>(10);
  b[0] = -18;
  b[1] = -15;
  
};



int main(int argc, const char** argv) {
  vect_n<double> c(4);
  mat<double,mat_structure::rectangular> A(3,4);
  vect_n<double> b(3);
  
#if 0
  
  serialization::xml_oarchive out("lp_problems/lp_template.xml");
  out & RK_SERIAL_SAVE_WITH_NAME(c)
      & RK_SERIAL_SAVE_WITH_NAME(A)
      & RK_SERIAL_SAVE_WITH_NAME(b);
  
#else
  
  if(argc < 2)
    return 1;
  
  std::vector< vect_n<double> > probs_c;
  std::vector< mat<double,mat_structure::rectangular> > probs_A;
  std::vector< vect_n<double> > probs_b;
  
  for(int i = 1; i < argc; ++i) {
    serialization::xml_iarchive in(std::string("lp_problems/") + std::string(argv[i]));
    in & RK_SERIAL_LOAD_WITH_NAME(c)
       & RK_SERIAL_LOAD_WITH_NAME(A)
       & RK_SERIAL_LOAD_WITH_NAME(b);
       
    if((A.get_col_count() == c.size()) &&
       (A.get_row_count() == b.size())) {
      probs_c.push_back(c);
      probs_A.push_back(A);
      probs_b.push_back(b);
    };
  };
  
  
  for(std::size_t i = 0; i < probs_c.size(); ++i) {
    
    std::cout << "********************************************************************" << std::endl;
    std::cout << " Solving problem number " << i << std::endl;
    
    
    vect_n<double> x(probs_c[i].size(),1.0);
    vect_n<double> l(probs_c[i].size(),0.0);
    vect_n<double> u(probs_c[i].size(),std::numeric_limits<double>::infinity());
    x[0] = 20.0; x[1] = 20.0;
    
    try {
      ReaK::optim::simplex_method(probs_A[i],probs_b[i],-probs_c[i],x,l,u,1e-6);
      std::cout << "   Simplex method produced the following solution vector: " << x << " with cost: " << (probs_c[i] * x) << std::endl;
    } catch(std::exception& e) {
      std::cout << "   Simplex method failed with message: " << e.what() << std::endl;
    };
    
    try {
      ReaK::optim::mehrotra_method(probs_A[i],probs_b[i],probs_c[i],x,1e-6);
      std::cout << "   Mehrotra method produced the following solution vector: " << x << " with cost: " << (probs_c[i] * x) << std::endl;
    } catch(std::exception& e) {
      std::cout << "   Mehrotra method failed with message: " << e.what() << std::endl;
    };
  };
  
  
#endif
  
  return 0;
};








