
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

#include "rk_mat_num.hpp"


#if 0
namespace ReaK {



template <class T>
unsigned int GaussJordanNoPivot(T* A, unsigned int N, unsigned int M, T* b, bool rowMajor)
{
  unsigned int c=0;
  unsigned int r=0;
  unsigned int rank=0;
  T s;
  if(rowMajor) {
    while((c<M) && (r<N)) {
      if (A[r*M+c] != T(0.0)) {
        for(unsigned int k=c+1;k<M;k++)
          A[r*M+k] = A[r*M+k] / A[r*M+c];
        b[r] /= A[r*M+c];
        A[r*M+c] = 1.0;
        for(unsigned int j=0;j<N;j++) {
          if ((j != r) && (A[j*M+c] != 0.0)) {
            s = A[j*M+c];
            for(unsigned int k=0;k<M;k++)
              A[j*M+k] -= s * A[r*M+k];
            b[j] -= s * b[r];
          };
        };
        rank++;
        c=0;
        r++;
      } else {
        c++;
        if((c==M) && (r<N)) {
          c=0;
          r++;
        };
      };
    };
  } else {
    while((c<M) && (r<N)) {
      if (A[c*N+r] != T(0.0)) {
        for(unsigned int k=c+1;k<M;k++)
          A[k*N+r] = A[k*N+r] / A[c*N+r];
        b[r] /= A[c*N+r];
        A[c*N+r] = 1.0;
        for(unsigned int j=0;j<N;j++) {
          if ((j != r) && (A[c*N+j] != T(0.0))) {
            s = A[c*N+j];
            for(unsigned int k=0;k<M;k++)
              A[k*N+j] -= s * A[k*N+r];
            b[j] -= s * b[r];
          };
        };
        rank++;
        c=0;
        r++;
      } else {
        c++;
        if((c==M) && (r<N)) {
          c=0;
          r++;
        };
      };
    };
  };
  return rank;
};

template <class T>
bool IntersectLinearSpaces(T* A1, unsigned int c1, T* p1, T* A2, unsigned int c2, T* p2, unsigned int N, T** AR, unsigned int* cR, T* pR, bool rowMajor)
{
  if(rowMajor) {
    T* Ap = new T[(c1+c2)*N];
    T* bp = new T[N];

    for(unsigned int i=0;i<N;i++) {
      for(unsigned int j=0;j<c1;j++)
        Ap[i*(c1+c2)+j] = A1[i*c1+j];
      for(unsigned int j=0;j<c2;j++)
        Ap[i*(c1+c2)+(j+c1)] = -A2[i*c2+j];
      bp[i] = p1[i] - p2[i];
    };
    unsigned int rank = GaussJordanNoPivot<T>(Ap,N,c1+c2,bp,rowMajor);
    if (rank < N) {
      for(unsigned int i=0;i<N;i++) {
        bool ZeroRow = true;
        for(unsigned int j=0;j<c1+c2;j++)
          ZeroRow &= (Ap[i*(c1+c2)+j] == T(0.0));
        if ((ZeroRow) && (pR[i] != T(0.0))) {
          delete[] Ap;
          delete[] bp;
          AR[0] = NULL;
          cR[0] = 0;
          return false;
        };
      };
    };
    unsigned int* DepVars = new unsigned int[rank];
    cR[0] = c1+c2-rank;
    unsigned int* IndepVars = new unsigned int[cR[0]];
    unsigned int d = 0;
    for(unsigned int i=0;i<c1+c2;i++) {
      bool Dependant = false;
      for(unsigned int j=0;j<N;j++) {
        if(fabs(Ap[j*(c1+c2)+i]) > T(0.000001)) {
          if((Dependant) || (Ap[j*(c1+c2)+i] != T(1.0))) {
            Dependant = false;
            break;
          } else
            Dependant = true;
        };
      };
      if(Dependant)
        DepVars[d++] = i;
      else
        IndepVars[i-d] = i;
    };
    T* M = new T[N*N];
    for(unsigned int i=0;i<N*N;i++)
      M[i] = 0.0;
    for(unsigned int i=0;i<rank;i++) {
      for(unsigned int j=0;j<N;j++) {
        if(Ap[j*(c1+c2)+DepVars[i]] == T(1.0)) {
          d = j;
          break;
        };
      };
      for(unsigned int j=0;j<N;j++) {
        if(DepVars[i] < c1)
          M[j*N+d] = A1[j*c1+DepVars[i]];
        else
          M[j*N+d] = A2[j*c2+DepVars[i]-c1];
      };
    };
    for(unsigned int i=0;i<N;i++) {
      pR[i] = 0.0;
      for(unsigned int j=0;j<N;j++)
        pR[i] += M[i*N+j] * bp[j];
      pR[i] += p1[i] + p2[i];
      pR[i] /= T(2.0);
    };
    AR[0] = new T[cR[0]*N];
    for(unsigned int i=0;i<N;i++) {
      for(unsigned int j=0;j<cR[0];j++) {
        if(IndepVars[j] < c1)
          AR[0][i*cR[0]+j] = A1[i*c1+IndepVars[j]];
        else
          AR[0][i*cR[0]+j] = A2[i*c2+IndepVars[j]-c1];
        for(unsigned int k=0;k<N;k++)
          AR[0][i*cR[0]+j] -= M[i*N+k] * Ap[k*(c1+c2)+IndepVars[j]];
        AR[0][i*cR[0]+j] /= T(2.0);
      };
    };
    delete[] Ap;
    delete[] bp;
    delete[] M;
    delete[] DepVars;
    delete[] IndepVars;
  } else {
    T* Ap = new T[(c1+c2)*N];
    T* bp = new T[N];

    for(unsigned int i=0;i<N;i++) {
      for(unsigned int j=0;j<c1;j++)
        Ap[j*N+i] = A1[j*N+i];
      for(unsigned int j=0;j<c2;j++)
        Ap[(j+c1)*N+i] = -A2[j*N+i];
      bp[i] = p1[i] - p2[i];
    };
    unsigned int rank = GaussJordanNoPivot<T>(Ap,N,c1+c2,bp,rowMajor);
    if (rank < N) {
      for(unsigned int i=0;i<N;i++) {
        bool ZeroRow = true;
        for(unsigned int j=0;j<c1+c2;j++)
          ZeroRow &= (Ap[j*N+i] == T(0.0));
        if ((ZeroRow) && (pR[i] != T(0.0))) {
          delete[] Ap;
          delete[] bp;
          AR[0] = NULL;
          cR[0] = 0;
          return false;
        };
      };
    };
    unsigned int* DepVars = new unsigned int[rank];
    cR[0] = c1+c2-rank;
    unsigned int* IndepVars = new unsigned int[cR[0]];
    unsigned int d = 0;
    for(unsigned int i=0;i<c1+c2;i++) {
      bool Dependant = false;
      for(unsigned int j=0;j<N;j++) {
        if(fabs(Ap[i*N+j]) > T(0.000001)) {
          if((Dependant) || (Ap[i*N+j] != T(1.0))) {
            Dependant = false;
            break;
          } else
            Dependant = true;
        };
      };
      if(Dependant)
        DepVars[d++] = i;
      else
        IndepVars[i-d] = i;
    };
    T* M = new T[N*N];
    for(unsigned int i=0;i<N*N;i++)
      M[i] = 0.0;
    for(unsigned int i=0;i<rank;i++) {
      for(unsigned int j=0;j<N;j++) {
        if(Ap[DepVars[i]*N+j] == T(1.0)) {
          d = j;
          break;
        };
      };
      for(unsigned int j=0;j<N;j++) {
        if(DepVars[i] < c1)
          M[d*N+j] = A1[DepVars[i]*N+j];
        else
          M[d*N+j] = A2[(DepVars[i]-c1)*N+j];
      };
    };
    for(unsigned int i=0;i<N;i++) {
      pR[i] = 0.0;
      for(unsigned int j=0;j<N;j++)
        pR[i] += M[j*N+i] * bp[j];
      pR[i] += p1[i] + p2[i];
      pR[i] /= T(2.0);
    };
    AR[0] = new T[cR[0]*N];
    for(unsigned int i=0;i<N;i++) {
      for(unsigned int j=0;j<cR[0];j++) {
        if(IndepVars[j] < c1)
          AR[0][j*N+i] = A1[IndepVars[j]*N+i];
        else
          AR[0][j*N+i] = A2[(IndepVars[j]-c1)*N+i];
        for(unsigned int k=0;k<N;k++)
          AR[0][j*N+i] -= M[k*N+i] * Ap[IndepVars[j]*N+k];
        AR[0][j*N+i] /= T(2.0);
      };
    };
    delete[] Ap;
    delete[] bp;
    delete[] M;
    delete[] DepVars;
    delete[] IndepVars;
  };
  return true;
};

};

#endif

