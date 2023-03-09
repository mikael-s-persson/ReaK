
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

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/mat_gaussian_elim.h"
#include "ReaK/math/lin_alg/mat_jacobi_method.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/lin_alg/mat_schur_decomp.h"
#include "ReaK/math/lin_alg/mat_svd_method.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "absl/log/log.h"

int main() {

  using namespace ReaK;

  using namespace std::chrono;

  unsigned int passed = 0;

  try {

    unsigned int inc(1);
    mat<double, mat_structure::symmetric> m_test(2.0, -1.0, 0.0, 2.0, -1.0,
                                                 2.0);
    mat<double, mat_structure::symmetric> m_inc(2.0, -1.0, 2.0);
    high_resolution_clock::time_point t1;
    std::array<high_resolution_clock::duration, 11> dt;

    std::ofstream out_stream;
    out_stream.open("performance_data.dat");
    out_stream
        << "N\tGauss\tPLU\tChol\tJac\tQR\tSymQR\tLDL\tSVD\tJac_E\tQR_E\tSVD_E"
        << std::endl;
    std::cout << "Recording performance..." << std::endl;
    for (unsigned int i = 3; i <= 1000; i += inc) {
      if (i == 50) {
        inc = 10;
        m_inc = get_block(m_test, 0, inc + 1);
      } else if (i == 100) {
        inc = 20;
        m_inc = get_block(m_test, 0, inc + 1);
      } else if (i == 500) {
        inc = 50;
        m_inc = get_block(m_test, 0, inc + 1);
      }

      m_test.set_row_count(i);
      set_block(m_test, m_inc, i - inc - 1);

      mat<double, mat_structure::square> m_gauss_inv(i);
      t1 = high_resolution_clock::now();
      invert_gaussian(m_test, m_gauss_inv, double(1E-15));
      dt[0] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::square> m_plu_inv(i);
      t1 = high_resolution_clock::now();
      invert_PLU(m_test, m_plu_inv, double(1E-15));
      dt[1] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::symmetric> m_cholesky_inv(i);
      t1 = high_resolution_clock::now();
      invert_Cholesky(m_test, m_cholesky_inv, double(1E-15));
      dt[2] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::symmetric> m_jacobi_pinv(i);
      t1 = high_resolution_clock::now();
      pseudoinvert_Jacobi(m_test, m_jacobi_pinv, double(1E-15));
      dt[3] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::square> m_qr_inv(i);
      t1 = high_resolution_clock::now();
      pseudoinvert_QR(m_test, m_qr_inv, double(1E-15));
      dt[4] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::square> m_symqr_inv(i);
      t1 = high_resolution_clock::now();
      pseudoinvert_SymQR(m_test, m_symqr_inv, double(1E-15));
      dt[5] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::square> m_ldl_inv(i);
      t1 = high_resolution_clock::now();
      invert_LDL(m_test, m_ldl_inv, double(1E-15));
      dt[6] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::square> m_svd_inv(i);
      t1 = high_resolution_clock::now();
      pseudoinvert_SVD(m_test, m_svd_inv, double(1E-15));
      dt[7] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::diagonal> m_jacobi_E(i);
      mat<double, mat_structure::square> m_jacobi_Q(i);
      t1 = high_resolution_clock::now();
      eigensolve_Jacobi(m_test, m_jacobi_E, m_jacobi_Q, double(1E-15));
      dt[8] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::diagonal> m_qr_E(i);
      mat<double, mat_structure::square> m_qr_Q(i);
      t1 = high_resolution_clock::now();
      eigensolve_SymQR(
          m_test, m_qr_Q, m_qr_E,
          double(
              1E-15));  // this iterates forever!! (well, for a long time at least)
      dt[9] = high_resolution_clock::now() - t1;

      mat<double, mat_structure::diagonal> m_svd_E(i);
      mat<double, mat_structure::square> m_svd_U(i);
      mat<double, mat_structure::square> m_svd_V(i);
      t1 = high_resolution_clock::now();
      decompose_SVD(m_test, m_svd_U, m_svd_E, m_svd_V, double(1E-15));
      dt[10] = high_resolution_clock::now() - t1;

      out_stream << i << "\t" << duration_cast<microseconds>(dt[0]).count()
                 << "\t" << duration_cast<microseconds>(dt[1]).count() << "\t"
                 << duration_cast<microseconds>(dt[2]).count() << "\t"
                 << duration_cast<microseconds>(dt[3]).count() << "\t"
                 << duration_cast<microseconds>(dt[4]).count() << "\t"
                 << duration_cast<microseconds>(dt[5]).count() << "\t"
                 << duration_cast<microseconds>(dt[6]).count() << "\t"
                 << duration_cast<microseconds>(dt[7]).count() << "\t"
                 << duration_cast<microseconds>(dt[8]).count() << "\t"
                 << duration_cast<microseconds>(dt[9]).count() << "\t"
                 << duration_cast<microseconds>(dt[10]).count() << std::endl;
      std::cout << i << std::endl;
    }
    std::cout << "Done!" << std::endl;
    out_stream.close();

  } catch (std::exception& e) {
    LOG(ERROR) << "An exception has occurred during the math_gen test: '"
               << e.what() << "'";
  } catch (...) {
    LOG(ERROR) << "An unexpected and unidentified exception has occurred "
                  "during the math_gen test.";
  }

  std::cout << "There were " << passed
            << " successful tests passed on the math_gen "
               "library, out of 45 possible successes.";

  return 0;
}
