
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

#include "base/defs.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

#include "mat_alg.hpp"

#include "mat_gaussian_elim.hpp"
#include "mat_cholesky.hpp"
#include "mat_jacobi_method.hpp"
#include "mat_qr_decomp.hpp"
#include "mat_svd_method.hpp"
#include "mat_schur_decomp.hpp"

#include "boost/date_time/posix_time/posix_time.hpp"

int main() {

  using namespace ReaK;

  unsigned int passed = 0;

  try {
  
    unsigned int inc(1);
    mat<double,mat_structure::symmetric> m_test(2.0,-1.0,0.0,2.0,-1.0,2.0);
    mat<double,mat_structure::symmetric> m_inc(2.0,-1.0,2.0);
    boost::posix_time::ptime t1;
    boost::posix_time::time_duration dt[11];

    std::ofstream out_stream;
    out_stream.open("performance_data.dat");
    out_stream << "N\tGauss\tPLU\tChol\tJac\tQR\tSymQR\tLDL\tSVD\tJac_E\tQR_E\tSVD_E" << std::endl;
    std::cout << "Recording performance..." << std::endl;
    for(unsigned int i=3;i<=1000;i += inc) {
      if(i == 50) {
	inc = 10; 
	m_inc = get_block(m_test,0,inc+1); 
      } else if(i == 100) {
	inc = 20; 
	m_inc = get_block(m_test,0,inc+1); 
      } else if(i == 500) {
	inc = 50; 
	m_inc = get_block(m_test,0,inc+1); 
      };

      m_test.set_row_count(i);
      set_block(m_test,m_inc,i-inc-1);
      
      mat<double,mat_structure::square> m_gauss_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_gaussian(m_test,m_gauss_inv,double(1E-15));
      dt[0] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::square> m_plu_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_PLU(m_test,m_plu_inv,double(1E-15));
      dt[1] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::symmetric> m_cholesky_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_Cholesky(m_test,m_cholesky_inv,double(1E-15));
      dt[2] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::symmetric> m_jacobi_pinv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_Jacobi(m_test,m_jacobi_pinv,double(1E-15));
      dt[3] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::square> m_qr_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_QR(m_test,m_qr_inv,double(1E-15));
      dt[4] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::square> m_symqr_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_SymQR(m_test,m_symqr_inv,double(1E-15));
      dt[5] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::square> m_ldl_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      invert_LDL(m_test,m_ldl_inv,double(1E-15));
      dt[6] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::square> m_svd_inv(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      pseudoinvert_SVD(m_test,m_svd_inv,double(1E-15));
      dt[7] = boost::posix_time::microsec_clock::local_time() - t1;;
      
      mat<double,mat_structure::diagonal> m_jacobi_E(i);
      mat<double,mat_structure::square> m_jacobi_Q(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      eigensolve_Jacobi(m_test,m_jacobi_E,m_jacobi_Q,double(1E-15));
      dt[8] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::diagonal> m_qr_E(i);
      mat<double,mat_structure::square> m_qr_Q(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      eigensolve_SymQR(m_test,m_qr_Q,m_qr_E,double(1E-15)); //this iterates forever!! (well, for a long time at least)
      dt[9] = boost::posix_time::microsec_clock::local_time() - t1;
      
      mat<double,mat_structure::diagonal> m_svd_E(i);
      mat<double,mat_structure::square> m_svd_U(i);
      mat<double,mat_structure::square> m_svd_V(i);
      t1 = boost::posix_time::microsec_clock::local_time();
      decompose_SVD(m_test,m_svd_U,m_svd_E,m_svd_V,double(1E-15));
      dt[10] = boost::posix_time::microsec_clock::local_time() - t1;
      
      out_stream << i << "\t" << dt[0].total_microseconds()
                      << "\t" << dt[1].total_microseconds()
                      << "\t" << dt[2].total_microseconds()
		      << "\t" << dt[3].total_microseconds()
		      << "\t" << dt[4].total_microseconds()
		      << "\t" << dt[5].total_microseconds()
		      << "\t" << dt[6].total_microseconds()
		      << "\t" << dt[7].total_microseconds()
		      << "\t" << dt[8].total_microseconds()
                      << "\t" << dt[9].total_microseconds()
                      << "\t" << dt[10].total_microseconds() << std::endl;
      std::cout << i << std::endl;

    };
    std::cout << "Done!" << std::endl;
    out_stream.close();
  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
  };
  
  RK_NOTICE(2,"There were " << passed << " successful tests passed on the math_gen library, out of 45 possible successes.");
  
  return 0;
};



