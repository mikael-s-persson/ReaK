
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

#pragma hdrstop

#ifndef RK_OptimizationC
#define RK_OptimizationC

#if 0
//The following functions have been ported to the new framework.

template <class T>
void WeightJacobian2PtsForward(void* UserData, TBasisMap Map, T *u, T *hx, T *hxx, T delta, T *jac, bool rowMajor)
{
  int i, j;
  T tmp;
  T d;

  for(j=0; j<Map.basis_count; ++j){
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    d=(1E-04)*Map.weights[j]; // force evaluation
    d=FABS(d);
    if(d<delta)
      d=delta;

    tmp=Map.weights[j];
    Map.weights[j]+=d;

    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,hxx,Map.output_count);

    Map.weights[j]=tmp; /* restore */

    d=1.0/d; /* invert so that divisions can be carried out faster as multiplications */

    if (rowMajor) {
      for(i=0; i<Map.output_count; ++i)
        jac[i*Map.basis_count+j]=(hxx[i]-hx[i])*d;
    } else {
      for(i=0; i<Map.output_count; ++i)
        jac[j*Map.output_count+i]=(hxx[i]-hx[i])*d;
    };
  };
};

template <class T>
void WeightJacobian2PtsCentral(void* UserData, TBasisMap Map, T *u, T *hxm, T *hxp, T delta, T *jac, bool rowMajor)
{
  int i, j;
  T tmp;
  T d;

  for(j=0; j<Map.basis_count; ++j){
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    d=(1E-04)*Map.weights[j]; // force evaluation
    d=FABS(d);
    if(d<delta)
      d=delta;

    tmp=Map.weights[j];
    Map.weights[j]-=d;
    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,hxm,Map.output_count);

    Map.weights[j]=tmp+d;
    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,hxp,Map.output_count);
    Map.weights[j]=tmp; /* restore */

    d=0.5/d; /* invert so that divisions can be carried out faster as multiplications */
    if (rowMajor) {
      for(i=0; i<Map.output_count; ++i)
        jac[i*Map.basis_count+j]=(hxp[i]-hxm[i])*d;
    } else {
      for(i=0; i<Map.output_count; ++i)
        jac[j*Map.output_count+i]=(hxp[i]-hxm[i])*d;
    };
  };
};

template <class T>
void Jacobian5Pts(void* UserData,TBasisMap Map,T* u,T du,T* Jac, bool rowMajor)
{
  if((u == NULL) || (Jac == NULL) || (du == 0.0)) return;

  T* Y0 = new T[Map.output_count];
  T* Y1 = new T[Map.output_count];
  T* Y2 = new T[Map.output_count];
  T* Y3 = new T[Map.output_count];
  for(unsigned int i=0;i<Map.input_count;i++) {
    u[i] -= 2.0*du;
    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,Y0,Map.output_count);
    u[i] += du;
    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,Y1,Map.output_count);
    u[i] += 2.0*du;
    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,Y2,Map.output_count);
    u[i] += du;
    Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,Y3,Map.output_count);
    if (rowMajor) {
      for(unsigned int j=0;j<Map.output_count;j++)
        Jac[j*Map.input_count+i] = (Y0[j] - 8.0*Y1[j] + 8.0*Y2[j] - Y3[j])/(12.0*du);
    } else {
      for(unsigned int j=0;j<Map.output_count;j++)
        Jac[i*Map.output_count+j] = (Y0[j] - 8.0*Y1[j] + 8.0*Y2[j] - Y3[j])/(12.0*du);
    };
    u[i] -= 2.0*du;
  };
  delete[] Y0;
  delete[] Y1;
  delete[] Y2;
  delete[] Y3;
  return;
};

template <class T>
void Jacobian5Pts(void* UserData,TNonLinMap Map,T* u,T du,T* Jac, bool rowMajor)
{
  if((u == NULL) || (Jac == NULL) || (du == 0.0)) return;

  T* Y0 = new T[Map.output_count];
  T* Y1 = new T[Map.output_count];
  T* Y2 = new T[Map.output_count];
  T* Y3 = new T[Map.output_count];
  for(unsigned int i=0;i<Map.input_count;i++) {
    u[i] -= 2.0*du;
    Map.func(UserData,u,Map.input_count,Y0,Map.output_count);
    u[i] += du;
    Map.func(UserData,u,Map.input_count,Y1,Map.output_count);
    u[i] += 2.0*du;
    Map.func(UserData,u,Map.input_count,Y2,Map.output_count);
    u[i] += du;
    Map.func(UserData,u,Map.input_count,Y3,Map.output_count);
    if (rowMajor) {
      for(unsigned int j=0;j<Map.output_count;j++)
        Jac[j*Map.input_count+i] = (Y0[j] - 8.0*Y1[j] + 8.0*Y2[j] - Y3[j])/(12.0*du);
    } else {
      for(unsigned int j=0;j<Map.output_count;j++)
        Jac[i*Map.output_count+j] = (Y0[j] - 8.0*Y1[j] + 8.0*Y2[j] - Y3[j])/(12.0*du);
    };
    u[i] -= 2.0*du;
  };
  delete[] Y0;
  delete[] Y1;
  delete[] Y2;
  delete[] Y3;
  return;
};

template <class T>
int SimplexMethod(T* A,T* b,unsigned int M,T* c,unsigned int N,T* x0,T* l,T* u)
{
  T* x = new T[N+M];
  T* b_G = new T[M];
  T* c_G = new T[N+M];
  T* c_B = new T[M];
  T* l_G = new T[N+M];
  T* u_G = new T[N+M];
  T* A_G = new T[(N+M)*M];
  T* A_B = new T[M*M];
  unsigned int* i_b = new unsigned int[M];
  unsigned int* i_n = new unsigned int[N];
  T* y = new T[M];
  T* d = new T[M];
  T w;
  T t_max;
  unsigned int enter_var;
  unsigned int leave_var;
  memmove(x,x0,N*sizeof(T));
  memmove(b_G,b,M*sizeof(T));
  memmove(l_G,l,N*sizeof(T));
  memmove(u_G,u,N*sizeof(T));
  memmove(A_G,A,N*M*sizeof(T));
  for(unsigned int i=0;i<N;i++)
    i_n[i] = i;
  for(unsigned int i=0;i<M;i++)
    i_b[i] = i + N;

  for(unsigned int i=0;i<N;i++) {
    if (x[i] < l[i])
      x[i] = l[i];
    if (x[i] > u[i])
      x[i] = u[i];
    c_G[i] = 0.0;
  };

  bool FeasibleStart = true;
  for(unsigned int i=0;i<M;i++) {
    x[N+i] = b_G[i];
    for(unsigned int j=0;j<N;j++)
      x[N+i] -= A[j*M+i] * x[j];
    FeasibleStart = (x[N+i] != 0.0);
    if (x[N+i] >= 0.0) {
      l_G[N+i] = 0.0;
      u_G[N+i] = RK_F_INF;
    } else {
      l_G[N+i] = -RK_F_INF;
      u_G[N+i] = 0.0;
    };
    c_G[N+i] = -1.0;
    c_B[i] = -1.0;
  };

  for(unsigned int i=0;i<M;i++)
    for(unsigned int j=0;j<M;j++)
      A_G[(N+i)*M+j] = 0.0;
  for(unsigned int i=0;i<M;i++)
    A_G[(N+i)*M+i] = 1.0;
  memmove(A_B,&A_G[N*M],M*M*sizeof(T));

  if (!FeasibleStart) {
    // First-Phase

    // Step 1
    memmove(y,c_B,M*sizeof(T));
    PLUDecomposition<T>(A_B,M,y,NULL,true);
    while(true) {
      // Step 2
      bool Optimal = true;
      T sum;
      for(unsigned int i=0;i<N;i++) {
        sum = 0.0;
        for(unsigned int j=0;j<M;j++)
          sum += y[j] * A_G[i_n[i]*M + j];
        if (((sum < c_G[i_n[i]]) && (x[i_n[i]] < u_G[i_n[i]])) || ((sum > c_G[i_n[i]]) && (x[i_n[i]] > l_G[i_n[i]]))) {
          enter_var = i;
          Optimal = false;
          break;
        };
      };
      if (Optimal)
        break;
      // Step 3
      memmove(d,&A_G[i_n[enter_var]*M],M*sizeof(T));
      PLUDecomposition<T>(A_B,M,d,NULL,false);
      // Step 4
      if (sum < c_G[i_n[enter_var]]) {
        t_max = RK_F_INF;
        leave_var = 0;
        for(unsigned int i=0;i<M;i++) {
          if(d[i] > 0.0) {
            if (t_max * d[i] > x[i_b[i]] - l_G[i_b[i]]) {
              t_max = (x[i_b[i]] - l_G[i_b[i]])/d[i];
              leave_var = i;
            };
          } else if (d[i] < 0.0) {
            if (t_max * d[i] < x[i_b[i]] - u_G[i_b[i]]) {
              t_max = (x[i_b[i]] - u_G[i_b[i]])/d[i];
              leave_var = i;
            };
          };
        };
        if (t_max > u_G[i_n[enter_var]] - x[i_n[enter_var]]) {
          t_max = u_G[i_n[enter_var]] - x[i_n[enter_var]];
          x[i_n[enter_var]] = u_G[i_n[enter_var]];
          for (unsigned int i=0;i<M;i++)
            x[i_b[i]] -= t_max * d[i];
          continue;
        } else if (t_max != RK_F_INF) {
          x[i_n[enter_var]] += t_max;
          for (unsigned int i=0;i<M;i++)
            x[i_b[i]] -= t_max * d[i];
          memmove(&A_B[leave_var*M],&A_G[i_n[enter_var]*M],M*sizeof(T));
          c_B[leave_var] = c_G[i_n[enter_var]];
          unsigned int i = i_b[leave_var];
          i_b[leave_var] = i_n[enter_var];
          i_n[enter_var] = i;
        } else {
          delete[] x;
          delete[] b_G;
          delete[] c_G;
          delete[] c_B;
          delete[] l_G;
          delete[] u_G;
          delete[] A_G;
          delete[] A_B;
          delete[] i_b;
          delete[] i_n;
          delete[] y;
          delete[] d;

          return OPT_ERR_UNBOUNDED;
        };
      } else if (sum > c_G[i_n[enter_var]]) {
        t_max = RK_F_INF;
        leave_var = 0;
        for(unsigned int i=0;i<M;i++) {
          if(d[i] > 0.0) {
            if (t_max * d[i] > u_G[i_b[i]] - x[i_b[i]]) {
              t_max = (u_G[i_b[i]] - x[i_b[i]])/d[i];
              leave_var = i;
            };
          } else if (d[i] < 0.0) {
            if (t_max * d[i] < l_G[i_b[i]] - x[i_b[i]]) {
              t_max = (l_G[i_b[i]] - x[i_b[i]])/d[i];
              leave_var = i;
            };
          };
        };
        if (t_max > x[i_n[enter_var]] - l_G[i_n[enter_var]]) {
          t_max = x[i_n[enter_var]] - l_G[i_n[enter_var]];
          x[i_n[enter_var]] = l_G[i_n[enter_var]];
          for (unsigned int i=0;i<M;i++)
            x[i_b[i]] += t_max * d[i];
          continue;
        } else if (t_max != RK_F_INF) {
          x[i_n[enter_var]] -= t_max;
          for (unsigned int i=0;i<M;i++)
            x[i_b[i]] += t_max * d[i];
          memmove(&A_B[leave_var*M],&A_G[i_n[enter_var]*M],M*sizeof(T));
          c_B[leave_var] = c_G[i_n[enter_var]];
          unsigned int i = i_b[leave_var];
          i_b[leave_var] = i_n[enter_var];
          i_n[enter_var] = i;
        } else {
          delete[] x;
          delete[] b_G;
          delete[] c_G;
          delete[] c_B;
          delete[] l_G;
          delete[] u_G;
          delete[] A_G;
          delete[] A_B;
          delete[] i_b;
          delete[] i_n;
          delete[] y;
          delete[] d;

          return OPT_ERR_UNBOUNDED;
        };
      };

      // Step 1
      memmove(y,c_B,M*sizeof(T));
      PLUDecomposition<T>(A_B,M,y,NULL,true);

    };
  };

  // Did the first phase succeed?
  FeasibleStart = true;
  for (unsigned int i=0;i<M;i++) {
    if (x[N+i] != 0.0) {
      FeasibleStart = false;
      break;
    };
  };
  if (!FeasibleStart) {
    delete[] x;
    delete[] b_G;
    delete[] c_G;
    delete[] c_B;
    delete[] l_G;
    delete[] u_G;
    delete[] A_G;
    delete[] A_B;
    delete[] i_b;
    delete[] i_n;
    delete[] y;
    delete[] d;

    return OPT_ERR_INFEASIBLE;
  };

  // Getting Rid of the Artificial Variables
  for(unsigned int i=0;((i<M) && (i_b[i] >= N));i++) {
    if(i_b[i] >= N) {
      memmove(y,&A_G[i_b[i]*M],M*sizeof(T));
      PLUDecomposition<T>(A_B,M,y,NULL,true);
      for(unsigned int j=0;j<N;j++) {
        T sum = 0.0;
        for(unsigned int k=0;k<M;k++)
          sum += y[k] * A_G[i_n[j]*M+k];
        if ((sum != 0.0) && (i_n[j] < N)) {
          memmove(&A_B[i*M],&A_G[i_n[j]*M],M*sizeof(T));
          unsigned int l = i_b[i];
          i_b[i] = i_n[j];
          i_n[j] = l;
          break;
        };
      };
    };
  };
  unsigned int RedundantCount = 0;
  for(unsigned int i=0;i<M;i++) {
    if(i_b[i] >= N) {
      // Must Delete Redundant Equation
      RedundantCount++;
      M--;
      T* tempA = new T[N*M];
      T* tempB = new T[M*M];
      T* tempb = new T[M];
      for(unsigned int j=0;j<i;j++) {
        memmove(&tempB[j*M],&A_B[j*(M+1)],(i_b[i]-N)*sizeof(T));
        memmove(&tempB[j*M+i_b[i]-N],&A_B[j*(M+1)+i_b[i]-N+1],(M-i_b[i]+N)*sizeof(T));
      };
      for(unsigned int j=i;j<M;j++) {
        memmove(&tempB[j*M],&A_B[(j+1)*(M+1)],(i_b[i]-N)*sizeof(T));
        memmove(&tempB[j*M+i_b[i]-N],&A_B[(j+1)*(M+1)+i_b[i]-N+1],(M-i_b[i]+N)*sizeof(T));
      };
      for(unsigned int j=0;j<N;j++) {
        memmove(&tempA[j*M],&A_G[j*(M+1)],(i_b[i]-N)*sizeof(T));
        memmove(&tempA[j*M+i_b[i]-N],&A_G[j*(M+1)+i_b[i]-N+1],(M-i_b[i]+N)*sizeof(T));
      };
      memmove(tempb,b_G,(i_b[i]-N)*sizeof(T));
      memmove(&tempb[i_b[i]-N],&b_G[i_b[i]-N+1],(M+N-i_b[i])*sizeof(T));
      delete[] b_G;
      b_G = tempb;
      delete[] A_G;
      A_G = tempA;
      delete[] A_B;
      A_B = tempB;
    };
  };

  // Prepare the variables for the second phase
  {
  T* temp = new T[N];
  memmove(temp,x,N*sizeof(T));
  delete[] x;
  x = temp;

  temp = new T[N];
  memmove(temp,c,N*sizeof(T));
  delete[] c_G;
  c_G = temp;

  c_B = new T[M];
  for(unsigned int i=0;i<M;i++)
    c_B[i] = c_G[i_b[i]];
  l_G = l;
  u_G = u;

  unsigned int* tempIB = new unsigned int[M];
  unsigned int* tempIN = new unsigned int[N-M];
  {
  unsigned int j=0;
  for(unsigned int i=0;i<M+RedundantCount;i++) {
    if(i_b[i] < N) {
      tempIB[j] = i_b[i];
      j++;
    };
  };
  j=0;
  for(unsigned int i=0;i<N;i++) {
    if(i_n[i] < N) {
      tempIN[j] = i_n[i];
      j++;
    };
  };
  };
  delete[] i_b;
  delete[] i_n;
  i_b = tempIB;
  i_n = tempIN;
  };

  // Second-Phase
  memmove(y,c_B,M*sizeof(T));
  PLUDecomposition<T>(A_B,M,y,NULL,true);
  while(true) {
    // Step 2
    bool Optimal = true;
    T sum;
    for(unsigned int i=0;i<N-M;i++) {
      sum = 0.0;
      for(unsigned int j=0;j<M;j++)
        sum += y[j] * A_G[i_n[i]*M + j];
      if (((sum < c_G[i_n[i]]) && (x[i_n[i]] < u_G[i_n[i]])) || ((sum > c_G[i_n[i]]) && (x[i_n[i]] > l_G[i_n[i]]))) {
        enter_var = i;
        Optimal = false;
        break;
      };
    };
    if (Optimal)
      break;
    // Step 3
    memmove(d,&A_G[i_n[enter_var]*M],M*sizeof(T));
    PLUDecomposition<T>(A_B,M,d,NULL,false);
    // Step 4
    if (sum < c_G[i_n[enter_var]]) {
      t_max = RK_F_INF;
      leave_var = 0;
      for(unsigned int i=0;i<M;i++) {
        if(d[i] > 0.0) {
          if (t_max * d[i] > x[i_b[i]] - l_G[i_b[i]]) {
            t_max = (x[i_b[i]] - l_G[i_b[i]])/d[i];
            leave_var = i;
          };
        } else if (d[i] < 0.0) {
          if (t_max * d[i] < x[i_b[i]] - u_G[i_b[i]]) {
            t_max = (x[i_b[i]] - u_G[i_b[i]])/d[i];
            leave_var = i;
          };
        };
      };
      if (t_max > u_G[i_n[enter_var]] - x[i_n[enter_var]]) {
        t_max = u_G[i_n[enter_var]] - x[i_n[enter_var]];
        x[i_n[enter_var]] = u_G[i_n[enter_var]];
        for (unsigned int i=0;i<M;i++)
          x[i_b[i]] -= t_max * d[i];
        continue;
      } else if (t_max != RK_F_INF) {
        x[i_n[enter_var]] += t_max;
        for (unsigned int i=0;i<M;i++)
          x[i_b[i]] -= t_max * d[i];
        memmove(&A_B[leave_var*M],&A_G[i_n[enter_var]*M],M*sizeof(T));
        c_B[leave_var] = c_G[i_n[enter_var]];
        unsigned int i = i_b[leave_var];
        i_b[leave_var] = i_n[enter_var];
        i_n[enter_var] = i;
      } else {
        delete[] x;
        delete[] b_G;
        delete[] c_G;
        delete[] c_B;
        delete[] l_G;
        delete[] u_G;
        delete[] A_G;
        delete[] A_B;
        delete[] i_b;
        delete[] i_n;
        delete[] y;
        delete[] d;

        return OPT_ERR_UNBOUNDED;
      };
    } else if (sum > c_G[i_n[enter_var]]) {
      t_max = RK_F_INF;
      leave_var = 0;
      for(unsigned int i=0;i<M;i++) {
        if(d[i] > 0.0) {
          if (t_max * d[i] > u_G[i_b[i]] - x[i_b[i]]) {
            t_max = (u_G[i_b[i]] - x[i_b[i]])/d[i];
            leave_var = i;
          };
        } else if (d[i] < 0.0) {
          if (t_max * d[i] < l_G[i_b[i]] - x[i_b[i]]) {
            t_max = (l_G[i_b[i]] - x[i_b[i]])/d[i];
            leave_var = i;
          };
        };
      };
      if (t_max > x[i_n[enter_var]] - l_G[i_n[enter_var]]) {
        t_max = x[i_n[enter_var]] - l_G[i_n[enter_var]];
        x[i_n[enter_var]] = l_G[i_n[enter_var]];
        for (unsigned int i=0;i<M;i++)
          x[i_b[i]] += t_max * d[i];
        continue;
      } else if (t_max != RK_F_INF) {
        x[i_n[enter_var]] -= t_max;
        for (unsigned int i=0;i<M;i++)
          x[i_b[i]] += t_max * d[i];
        memmove(&A_B[leave_var*M],&A_G[i_n[enter_var]*M],M*sizeof(T));
        c_B[leave_var] = c_G[i_n[enter_var]];
        unsigned int i = i_b[leave_var];
        i_b[leave_var] = i_n[enter_var];
        i_n[enter_var] = i;
      } else {
        delete[] x;
        delete[] b_G;
        delete[] c_G;
        delete[] c_B;
        delete[] l_G;
        delete[] u_G;
        delete[] A_G;
        delete[] A_B;
        delete[] i_b;
        delete[] i_n;
        delete[] y;
        delete[] d;

        return OPT_ERR_UNBOUNDED;
      };
    };

    // Step 1
    memmove(y,c_B,M*sizeof(T));
    PLUDecomposition<T>(A_B,M,y,NULL,true);
  };

  memmove(x0,x,N*sizeof(T));
  delete[] x;
  delete[] b_G;
  delete[] c_G;
  delete[] c_B;
  delete[] l_G;
  delete[] u_G;
  delete[] A_G;
  delete[] A_B;
  delete[] i_b;
  delete[] i_n;
  delete[] y;
  delete[] d;

  return OPT_OPTIMAL;
};

template <class T>
void LinearLeastSquare(T* A,unsigned int N,unsigned int M,T* x,T* b)
{
  T* Q = new T[N*M];
  T* R = new T[M*M];
  QRDecomposition<T>(A,N,M,Q,R);
  T* y = new T[M];
  for(unsigned int i=0;i<M;i++) {
    y[i] = 0.0;
    for(unsigned int j=0;j<N;j++)
      y[i] += Q[i*N+j] * b[j];
  };
  for(int i=M-1;i>=0;i--) {
    if (R[i*M+i] != 0.0) {
      x[i] = y[i];
      for(unsigned int j=M-1;j>i;j--)
        x[i] += R[j*M+i] * x[j];
      x[i] /= R[i*M+i];
    } else
      x[i] = 0.0;
  };
  delete[] Q;
  delete[] R;
  delete[] y;
  return;
};

template <class T>
void LinearLeastSquare(void* UserData,TBasisMap Basis,T* input,T* output)
{
  T* A = new T[Basis.basis_count*Basis.output_count];
  for(unsigned int i=0;i<Basis.basis_count;i++)
    Basis.weights[i] = 0.0;
  for(unsigned int i=0;i<Basis.basis_count;i++) {
    Basis.weights[i] = 1.0;
    Basis.func(UserData,input,Basis.input_count,Basis.weights,Basis.basis_count,&A[i*Basis.output_count],Basis.output_count);
    Basis.weights[i] = 0.0;
  };
  LinearLeastSquare(A,Basis.output_count,Basis.basis_count,Basis.weights,output);
  delete[] A;
  return;
};

template <class T>
void CheckWeightJacobian(void* UserData, TBasisMap Map, FBasisFunction Jacf, T* u, T *err)
{
  T factor=100.0;
  T one=1.0;
  T zero=0.0;
  T *fvec, *fjac, *pp, *fvecp, *buf;

  int i, j;
  T eps, epsf, temp, epsmch;
  T epslog;
  int fvec_sz=Map.output_count, fjac_sz=Map.output_count*Map.basis_count, pp_sz=Map.basis_count, fvecp_sz=Map.output_count;

  epsmch=FLT_EPSILON;
  eps=sqrt(epsmch);

  buf= new T[fvec_sz + fjac_sz + pp_sz + fvecp_sz];
  fvec=buf;
  fjac=fvec+fvec_sz;
  pp=fjac+fjac_sz;
  fvecp=pp+pp_sz;

  /* compute fvec=func(p) */
  Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,fvec,Map.output_count);

  /* compute the jacobian at p */
  Jacf(UserData,u,Map.input_count,Map.weights,Map.basis_count,fjac,Map.output_count*Map.basis_count);

  /* compute pp */
  for(j=0; j<Map.basis_count; ++j){
    temp=eps*FABS(Map.weights[j]);
    if(temp==zero) temp=eps;
    pp[j]=Map.weights[j]+temp;
  };

  /* compute fvecp=func(pp) */
  Map.func(UserData,u,Map.input_count,pp,Map.basis_count,fvecp,Map.output_count);

  epsf=factor*epsmch;
  epslog=log10(eps);

  for(i=0; i<Map.output_count; ++i)
    err[i]=zero;

  for(j=0; j<Map.basis_count; ++j){
    temp=FABS(Map.weights[j]);
    if(temp==zero) temp=one;

    for(i=0; i<Map.output_count; ++i)
      err[i]+=temp*fjac[i*Map.basis_count+j];
  };

  for(i=0; i<Map.output_count; ++i){
    temp=one;
    if(fvec[i]!=zero && fvecp[i]!=zero && FABS(fvecp[i]-fvec[i])>=epsf*FABS(fvec[i]))
        temp=eps*FABS((fvecp[i]-fvec[i])/eps - err[i])/(FABS(fvec[i])+FABS(fvecp[i]));
    err[i]=one;
    if(temp>epsmch && temp<eps)
        err[i]=(log10(temp) - epslog)/epslog;
    if(temp>=eps) err[i]=zero;
  };

  delete[] buf;

  return;
};

template <class T>
int ComputeCovariance(T* JtJ, T* C, T sumsq, unsigned int N, unsigned int M)
{
  int rnk;
  T fact;

  SetToIdentity(C,M);
  //rnk=LUINVERSE(JtJ, C, M);
  T* U = new T[2*M*M + M];
  T* V = U + M*M;
  T* E = V + M*M;
  SVDecomposition<T>(JtJ,M,M,U,E,V);
  rnk = SVDNumericalRank<T>(E,M,M);
  PseudoInvertSVD<T>(U,E,V,M,M,C);
  delete[] U;

  if(!rnk) return 0;

  fact=sumsq/(N-rnk);
  for(int i=0; i<M*M; ++i)
    C[i]*=fact;

  return rnk;
};
#endif

template <class T>
int LMNonLinLsq(void* UserData, TBasisMap Map, FBasisFunction Jacf, T* u, T* x, unsigned int itmax, T tau, T epsj, T epsw, T epsx, T info[9], T* covar)
{
  /* Check if the problem is defined properly */
  if ((!Jacf) ||
      (!u) ||
      (!x) ||
      (Map.output_count < Map.basis_count) ||
      (itmax <= 1))
    return OPT_ERR_IMPROPER;

  int i, j, k, l;
  int worksz, issolved;
  /* temp work arrays */
  T *e,          /* nx1 */
        *work,
        *hx,         /* \hat{x}_i, nx1 */
        *jacTe,      /* J^T e_i mx1 */
        *jac,        /* nxm */
        *jacTjac,    /* mxm */
        *Dp,         /* mx1 */
        *diag_jacTjac,   /* diagonal of J^T J, mx1 */
        *pDp;        /* p + Dp, mx1 */

  T mu,  /* damping constant */
        tmp; /* mainly used in matrix & vector multiplications */
  T epsw_sq, p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
  T p_L2, Dp_L2=FLT_MAX, dF, dL;
  T init_p_eL2;
  int nu=2, nu2, stop, nfev, njev=0;
  const unsigned int nm=Map.output_count*Map.basis_count;

  mu=jacTe_inf=0.0; /* -Wall */

  if(tau <= 0.0) tau=LM_INIT_MU;
  if(epsj <= 0.0) epsj=LM_STOP_THRESH;
  if(epsw <= 0.0) epsw=LM_STOP_THRESH;
  epsw_sq = epsw*epsw;
  if(epsx <= 0.0) epsx=LM_STOP_THRESH;

  worksz= 2*Map.output_count+4*Map.basis_count + Map.output_count*Map.basis_count + Map.basis_count*Map.basis_count;
  work= new T[worksz]; /* allocate a big chunk in one step */

  /* set up work arrays */
  e=work;
  hx=e + Map.output_count;
  jacTe=hx + Map.output_count;
  jac=jacTe + Map.basis_count;
  jacTjac=jac + nm;
  Dp=jacTjac + Map.basis_count*Map.basis_count;
  diag_jacTjac=Dp + Map.basis_count;
  pDp=diag_jacTjac + Map.basis_count;

  /* compute e=x - f(p) and its L2 norm */
  Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,hx,Map.output_count); nfev=1;
  for(i=0, p_eL2=0.0; i<Map.output_count; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  };
  init_p_eL2=p_eL2;

  for(k=stop=0; k<itmax && !stop; ++k){
    /* Note that p and e have been updated at a previous iteration */

    if(p_eL2<=epsx){ /* error is small */
      stop=6;
      break;
    };

    /* Compute the jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
     * Since J^T J is symmetric, its computation can be speeded up by computing
     * only its upper triangular part and copying it to the lower part
     */

    Jacf(UserData,u,Map.input_count,Map.weights,Map.basis_count,jac,nm); njev++;
    //Note: Jacobian is calculated with (Map.basis_count) columns and (Map.output_count) rows

    /* J^T J, J^T e */
    //if(nm<__BLOCKSZ__SQ){ // this is a small problem
      /* This is the straightforward way to compute J^T J, J^T e. However, due to
       * its noncontinuous memory access pattern, it incures many cache misses when
       * applied to large minimization problems (i.e. problems involving a large
       * number of free variables and measurements), in which J is too large to
       * fit in the L1 cache. For such problems, a cache-efficient blocking scheme
       * is preferable.
       *
       * Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
       * performance problem.
       *
       * On the other hand, the straightforward algorithm is faster on small
       * problems since in this case it avoids the overheads of blocking.
       */

      for(i=0; i<Map.basis_count; ++i){
        for(j=i; j<Map.basis_count; ++j){

          for(l=0, tmp=0.0; l<Map.output_count; ++l)
            tmp+=jac[i*Map.output_count+l]*jac[j*Map.output_count+l];

		      /* store tmp in the corresponding upper and lower part elements */
          jacTjac[j*Map.basis_count+i]=jacTjac[i*Map.basis_count+j]=tmp;
        };

        /* J^T e */
        for(l=0, tmp=0.0; l<Map.output_count; ++l)
          tmp+=jac[i*Map.output_count+l]*e[l];
        jacTe[i]=tmp;
      };
    //}; else { // this is a large problem
      ///* Cache efficient computation of J^T J based on blocking
      // */
      //TRANS_MAT_MAT_MULT(jac, jacTjac, n, m);

      ///* cache efficient computation of J^T e */
      //for(i=0; i<Map.basis_count; ++i)
      //  jacTe[i].q=0.0;

      //for(i=0; i<Map.output_count; ++i){
      //  T *jacrow;

      //  for(l=0, jacrow=jac+i*Map.basis_count, tmp.q=e[i].q; l<Map.basis_count; ++l)
      //    jacTe[l].q+=jacrow[l].q*tmp.q;
      //}
    //}

	  /* Compute ||J^T e||_inf and ||p||^2 */
    for(i=0, p_L2=jacTe_inf=0.0; i<Map.basis_count; ++i){
      if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;

      diag_jacTjac[i]=jacTjac[i*Map.basis_count+i]; /* save diagonal entries so that augmentation can be later canceled */
      p_L2+=Map.weights[i]*Map.weights[i];
    };
    //p_L2=sqrt(p_L2);


    /*Could Add Here a callback for the current estimate*/


    /* check for convergence */
    if((jacTe_inf <= epsj)){
      Dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    };

   /* compute initial damping factor */
    if(k==0){
      for(i=0, tmp=FLT_MIN; i<Map.basis_count; ++i)
        if(diag_jacTjac[i]>tmp) tmp=diag_jacTjac[i]; /* find max diagonal element */
      mu=tau*tmp;
    };

    /* determine increment using adaptive damping */
    while(1){
      /* augment normal equations */
      for(i=0; i<Map.basis_count; ++i)
        jacTjac[i*Map.basis_count+i]+=mu;

      /* solve augmented equations */
      /* use the LU included with levmar */
      //issolved=AX_EQ_B_LU(jacTjac, jacTe, Dp, Map.basis_count);
      memmove(Dp,jacTe,Map.basis_count*sizeof(T));
      issolved=PLUDecomposition<T>(jacTjac,Map.basis_count,Dp,NULL);

      if(issolved){
        /* compute weight's new estimate and ||Dp||^2 */
        for(i=0, Dp_L2=0.0; i<Map.basis_count; ++i){
          pDp[i]=Map.weights[i] + (tmp=Dp[i]);
          Dp_L2+=tmp*tmp;
        };
        //Dp_L2=sqrt(Dp_L2);

        if(Dp_L2<=epsw_sq*p_L2){ /* relative change in p is small, stop */
        //if(Dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
          stop=2;
          break;
        };

        if(Dp_L2>=(p_L2+epsw)/(FLT_EPSILON*FLT_EPSILON)){ /* almost singular */
        //if(Dp_L2>=(p_L2+eps2)/CNST(EPSILON)){ /* almost singular */
          stop=4;
          break;
        };

        Map.func(UserData,u,Map.input_count,pDp,Map.basis_count,hx,Map.output_count); nfev++;
        for(i=0, pDp_eL2=0.0; i<Map.output_count; ++i){ /* compute ||e(pDp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pDp_eL2+=tmp*tmp;
        };

        for(i=0, dL=0.0; i<Map.basis_count; ++i)
          dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

        dF=p_eL2-pDp_eL2;

        if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
          tmp=(2.0*dF/dL-1.0);
          tmp=1.0-tmp*tmp*tmp;
          mu=mu*( (tmp>=ONE_THIRD)? tmp : ONE_THIRD );
          nu=2;

          for(i=0 ; i<Map.basis_count; ++i) /* update p's estimate */
            Map.weights[i]=pDp[i];

          for(i=0; i<Map.output_count; ++i) /* update e and ||e||_2 */
            e[i]=hx[i];
          p_eL2=pDp_eL2;
          break;
        };
      };

      /* if this point is reached, either the linear system could not be solved or
       * the error did not reduce; in any case, the increment must be rejected
       */

      mu*=nu;
      nu2=nu<<1; // 2*nu;
      if(nu2<=nu){ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
        stop=5;
        break;
      };
      nu=nu2;

      for(i=0; i<Map.basis_count; ++i) /* restore diagonal J^T J entries */
        jacTjac[i*Map.basis_count+i]=diag_jacTjac[i];
    }; /* inner loop */
  };

  if(k>=itmax) stop=3;

  for(i=0; i<Map.basis_count; ++i) /* restore diagonal J^T J entries */
    jacTjac[i*Map.basis_count+i]=diag_jacTjac[i];

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=jacTe_inf;
    info[3]=Dp_L2;
    for(i=0, tmp=FLT_MIN; i<Map.basis_count; ++i)
      if(tmp<jacTjac[i*Map.basis_count+i]) tmp=jacTjac[i*Map.basis_count+i];
    info[4]=mu/tmp;
    info[5]=(T)k;
    info[6]=(T)stop;
    info[7]=(T)nfev;
    info[8]=(T)njev;
  };

  /* covariance matrix */
  if(covar)
    ComputeCovariance<T>(jacTjac, covar, p_eL2, Map.output_count, Map.basis_count);

  delete[] work;

  return (stop!=4)?  k : -1;
};

template <class T>
int LMNonLinLsq(void* UserData, TBasisMap Map, T* u, T* x, int itmax, T tau, T epsj, T epsw, T epsx, T delta, T info[9], T *covar)
{
  /* Check if the problem is defined properly */
  if ((!u) ||
      (!x) ||
      (Map.output_count < Map.basis_count) ||
      (itmax <= 1))
    return OPT_ERR_IMPROPER;

  int i, j, k, l;
  int worksz, issolved;
  /* temp work arrays */
  T *work,
        *e,          /* nx1 */
        *hx,         /* \hat{x}_i, nx1 */
        *jacTe,      /* J^T e_i mx1 */
        *jac,        /* nxm */
        *jacTjac,    /* mxm */
        *Dp,         /* mx1 */
        *diag_jacTjac,   /* diagonal of J^T J, mx1 */
        *pDp,        /* p + Dp, mx1 */
        *wrk;        /* nx1 */

  int using_ffdif=1;
  T *wrk2=NULL; /* nx1, used for differentiating with central differences only */

  T mu,  /* damping constant */
        tmp; /* mainly used in matrix & vector multiplications */
  T p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
  T p_L2, Dp_L2=FLT_MAX, dF, dL;
  T epsw_sq;
  T init_p_eL2;
  int nu, nu2, stop, nfev, njap=0, K=(Map.basis_count>=10)? Map.basis_count: 10, updjac, updp=1, newjac;
  const int nm=Map.output_count*Map.basis_count;

  mu=jacTe_inf=p_L2=0.0; /* -Wall */
  stop=updjac=newjac=0; /* -Wall */

  if(tau <= 0.0) tau = LM_INIT_MU;
  if(epsj <= 0.0) epsj = LM_STOP_THRESH;
  if(epsw <= 0.0) epsw = LM_STOP_THRESH;
  epsw_sq = epsw * epsw;
  if(epsx <= 0.0) epsx = LM_STOP_THRESH;
  if(delta <= 0.0) delta = LM_DIFF_DELTA;

  worksz= 3*Map.output_count+4*Map.basis_count + Map.output_count*Map.basis_count + Map.basis_count*Map.basis_count;
  work= new T[worksz]; /* allocate a big chunk in one step */

  /* set up work arrays */
  e=work;
  hx=e + Map.output_count;
  jacTe=hx + Map.output_count;
  jac=jacTe + Map.basis_count;
  jacTjac=jac + nm;
  Dp=jacTjac + Map.basis_count*Map.basis_count;
  diag_jacTjac=Dp + Map.basis_count;
  pDp=diag_jacTjac + Map.basis_count;
  wrk=pDp + Map.basis_count;

  /* compute e=x - f(p) and its L2 norm */
  Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,hx,Map.output_count); nfev=1;
  for(i=0, p_eL2=0.0; i<Map.output_count; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  };
  init_p_eL2=p_eL2;

  nu=20; /* force computation of J */
  using_ffdif=1; /* force forward difference formula */

  for(k=0; k<itmax; ++k){
    /* Note that p and e have been updated at a previous iteration */

    if(p_eL2<=epsx){ /* error is small */
      stop=6;
      break;
    };

    /* Compute the jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
     * The symmetry of J^T J is again exploited for speed
     */

    if((updp && nu>16) || updjac==K){ /* compute difference approximation to J */
      if(using_ffdif){ /* use forward differences */
        WeightJacobian2PtsForward<T>(UserData, Map, u, hx, wrk, delta, jac);
        ++njap; nfev+=Map.basis_count;
      } else { /* use central differences */
        WeightJacobian2PtsCentral<T>(UserData, Map, u, wrk, wrk2, delta, jac);
        ++njap; nfev+=2*Map.basis_count;
      };
      nu=2; updjac=0; updp=0; newjac=1;
    };

    if(newjac){ /* jacobian has changed, recompute J^T J, J^t e, etc */
      newjac=0;

      /* J^T J, J^T e */
      //if(nm<=__BLOCKSZ__SQ){ // this is a small problem
        /* This is the straightforward way to compute J^T J, J^T e. However, due to
         * its noncontinuous memory access pattern, it incures many cache misses when
         * applied to large minimization problems (i.e. problems involving a large
         * number of free variables and measurements), in which J is too large to
         * fit in the L1 cache. For such problems, a cache-efficient blocking scheme
         * is preferable.
         *
         * Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
         * performance problem.
         *
         * On the other hand, the straightforward algorithm is faster on small
         * problems since in this case it avoids the overheads of blocking.
         */

        for(i=0; i<Map.basis_count; ++i) {
          for(j=i; j<Map.basis_count; ++j) {

            for(l=0, tmp=0.0; l<Map.output_count; ++l)
              tmp+=jac[i*Map.output_count+l]*jac[j*Map.output_count+l];

            jacTjac[j*Map.basis_count+i]=jacTjac[i*Map.basis_count+j]=tmp;
          };

          /* J^T e */
          for(l=0, tmp=0.0; l<Map.output_count; ++l)
            tmp+=jac[i*Map.output_count+l]*e[l];
          jacTe[i]=tmp;
        };
      //}; else{ // this is a large problem
        /* Cache efficient computation of J^T J based on blocking
         */
        //TRANS_MAT_MAT_MULT(jac, jacTjac, Map.output_count, Map.basis_count);

        /* cache efficient computation of J^T e */
        //for(i=0; i<Map.basis_count; ++i)
        //  jacTe[i].q=0.0;

        //for(i=0; i<Map.output_count; ++i){
        //  T *jacrow;

        //  for(l=0, jacrow=jac+i*Map.basis_count, tmp.q=e[i].q; l<Map.basis_count; ++l)
        //    jacTe[l].q+=jacrow[l].q*tmp.q;
        //};
      //};

      /* Compute ||J^T e||_inf and ||p||^2 */
      for(i=0, p_L2=jacTe_inf=0.0; i<Map.basis_count; ++i){
        if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;

        diag_jacTjac[i]=jacTjac[i*Map.basis_count+i]; /* save diagonal entries so that augmentation can be later canceled */
        p_L2+=Map.weights[i]*Map.weights[i];
      };
      //p_L2=sqrt(p_L2);
    };

    /* could add callback for current iteration */

    /* check for convergence */
    if((jacTe_inf <= epsj)){
      Dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    };

   /* compute initial damping factor */
    if(k==0){
      for(i=0, tmp=FLT_MIN; i<Map.basis_count; ++i)
        if(diag_jacTjac[i]>tmp) tmp=diag_jacTjac[i]; /* find max diagonal element */
      mu=tau*tmp;
    };

    /* determine increment using adaptive damping */

    /* augment normal equations */
    for(i=0; i<Map.basis_count; ++i)
      jacTjac[i*Map.basis_count+i]+=mu;

    /* solve augmented equations */
    /* use the LU included with levmar */
    //issolved=AX_EQ_B_LU(jacTjac, jacTe, Dp, Map.basis_count);
    memmove(Dp,jacTe,Map.basis_count*sizeof(T));
    issolved=PLUDecomposition<T>(jacTjac,Map.basis_count,Dp,NULL);

    if(issolved){
    /* compute p's new estimate and ||Dp||^2 */
      for(i=0, Dp_L2=0.0; i<Map.basis_count; ++i){
        pDp[i]=Map.weights[i] + (tmp=Dp[i]);
        Dp_L2+=tmp*tmp;
      };
      //Dp_L2=sqrt(Dp_L2);

      if(Dp_L2<=epsw_sq*p_L2){ /* relative change in p is small, stop */
      //if(Dp_L2<=eps2*(p_L2 + eps2)){ /* relative change in p is small, stop */
        stop=2;
        break;
      };

      if(Dp_L2>=(p_L2+epsw)/(FLT_EPSILON*FLT_EPSILON)){ /* almost singular */
      //if(Dp_L2>=(p_L2+eps2)/CNST(EPSILON)){ /* almost singular */
        stop=4;
        break;
      };

      /* evaluate function at p + Dp */
      Map.func(UserData,u,Map.input_count,pDp,Map.basis_count,wrk,Map.output_count); ++nfev;
      for(i=0, pDp_eL2=0.0; i<Map.output_count; ++i){ /* compute ||e(pDp)||_2 */
        tmp=x[i]-wrk[i];
        pDp_eL2+=tmp*tmp;
      };

      dF=p_eL2-pDp_eL2;
      if(updp || dF>0.0){ /* update jac */
        for(i=0; i<Map.output_count; ++i){
          for(l=0, tmp=0.0; l<Map.basis_count; ++l)
            tmp+=jac[l*Map.output_count+i]*Dp[l]; /* (J * Dp)[i] */
          tmp=(wrk[i] - hx[i] - tmp)/Dp_L2; /* (f(p+dp)[i] - f(p)[i] - (J * Dp)[i])/(dp^T*dp) */
          for(j=0; j<Map.basis_count; ++j)
            jac[j*Map.output_count+i]+=tmp*Dp[j];
        };
        ++updjac;
        newjac=1;
      };

      for(i=0, dL=0.0; i<Map.basis_count; ++i)
        dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

      if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
        dF=(2.0*dF/dL-1.0);
        tmp=dF*dF*dF;
        tmp=1.0-tmp*tmp*dF;
        mu=mu*( (tmp>=ONE_THIRD)? tmp : ONE_THIRD );
        nu=2;

        for(i=0 ; i<Map.basis_count; ++i) /* update p's estimate */
          Map.weights[i]=pDp[i];

        for(i=0; i<Map.output_count; ++i){ /* update e, hx and ||e||_2 */
          e[i]=x[i]-wrk[i];
          hx[i]=wrk[i];
        };
        p_eL2=pDp_eL2;
        updp=1;
        continue;
      };
    };

    /* if this point is reached, either the linear system could not be solved or
     * the error did not reduce; in any case, the increment must be rejected
     */

    mu*=nu;
    nu2=nu<<1; // 2*nu;
    if(nu2<=nu){ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
      stop=5;
      break;
    };
    nu=nu2;

    for(i=0; i<Map.basis_count; ++i) /* restore diagonal J^T J entries */
      jacTjac[i*Map.basis_count+i]=diag_jacTjac[i];
  };

  if(k>=itmax) stop=3;

  for(i=0; i<Map.basis_count; ++i) /* restore diagonal J^T J entries */
    jacTjac[i*Map.basis_count+i]=diag_jacTjac[i];

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=jacTe_inf;
    info[3]=Dp_L2;
    for(i=0, tmp=FLT_MIN; i<Map.basis_count; ++i)
      if(tmp<jacTjac[i*Map.basis_count+i]) tmp=jacTjac[i*Map.basis_count+i];
    info[4]=mu/tmp;
    info[5]=(T)k;
    info[6]=(T)stop;
    info[7]=(T)nfev;
    info[8]=(T)njap;
  };

  /* covariance matrix */
  if(covar){
    ComputeCovariance<T>(jacTjac, covar, p_eL2, Map.output_count, Map.basis_count);
  };


  delete[] work;

  if(wrk2) delete[] (T*)wrk2;

  //Add return of proper error codes
  return (stop!=4)?  k : -1;
};

template <class T>
void ProjectOnBox(T *p, T *lb, T *ub, unsigned int N)
{
  int i;

  if(lb == NULL){ /* no lower bounds */
    if(ub == NULL) /* no upper bounds */
      return;
    else { /* upper bounds only */
      for(i=0; i<N; ++i)
        if(p[i]>ub[i]) p[i]=ub[i];
    };
  } else {
    if(ub == NULL) { /* lower bounds only */
      for(i=0; i<N; ++i)
        if(p[i]<lb[i]) p[i]=lb[i];
    } else { /* box bounds */
      for(i=0; i<N; ++i) {
        if(p[i]>ub[i])
          p[i]=ub[i];
        else if(p[i]<lb[i])
          p[i]=lb[i];
      };
    };
  };
};

template <class T>
int CheckBoxConsistency(T *lb, T *ub, unsigned int N)
{
  int i;

  if((lb==NULL) || (ub==NULL)) return 1;

  for(i=0; i<N; ++i)
    if(lb[i]>ub[i]) return 0;

  return 1;
};

template <class T>
int LMBoxNonLinLsq(void* UserData, TBasisMap Map, FBasisFunction Jacf, T* u, T *x, T *lb, T *ub, int itmax, T tau, T epsj, T epsw, T epsf, T info[9], T *covar)
{
  if ((!Jacf) ||
      (!u) ||
      (!x) ||
      (!lb) ||
      (!ub) ||
      (!CheckBoxConsistency<T>(lb, ub, Map.basis_count)) ||
      (Map.output_count < Map.basis_count) ||
      (itmax <= 1))
    return OPT_ERR_IMPROPER;

  int i, j, k, l;
  int issolved;
/* temp work arrays */
  T *work,
        *e,          /* nx1 */
        *hx,         /* \hat{x}_i, nx1 */
        *jacTe,      /* J^T e_i mx1 */
        *jac,        /* nxm */
        *jacTjac,    /* mxm */
        *Dp,         /* mx1 */
        *diag_jacTjac,   /* diagonal of J^T J, mx1 */
        *pDp;        /* p + Dp, mx1 */

  T mu,  /* damping constant */
        tmp; /* mainly used in matrix & vector multiplications */
  T p_eL2, jacTe_inf, pDp_eL2; /* ||e(p)||_2, ||J^T e||_inf, ||e(p+Dp)||_2 */
  T p_L2, Dp_L2=FLT_MAX, dF, dL;
  T epsw_sq;
  T init_p_eL2;
  int nu=2, nu2, stop, nfev, njev=0;
  const int nm=Map.output_count*Map.basis_count;

  /* variables for constrained LM */
  T alpha=1e-4, beta=0.9, gamma=0.99995, gamma_sq=gamma*gamma, rho=1e-8;
  T t, t0;
  T jacTeDp;
  T tmin=1e-12, tming=1e-18; /* minimum step length for LS and PG steps */
  T tini=1.0; /* initial step length for LS and PG steps */
  int nLMsteps=0, nLSsteps=0, nPGsteps=0, gprevtaken=0;
  int numactive;

  mu=jacTe_inf=t=0.0; /* -Wall */


  if(tau <= 0.0) tau = LM_INIT_MU;
  if(epsj <= 0.0) epsj = LM_STOP_THRESH;
  if(epsw <= 0.0) epsw = LM_STOP_THRESH;
  epsw_sq = epsw*epsw;
  if(epsf <= 0.0) epsf = LM_STOP_THRESH;

  work= new T[2*Map.output_count+4*Map.basis_count+Map.output_count*Map.basis_count+Map.basis_count*Map.basis_count]; /* allocate a big chunk in one step */

  /* set up work arrays */
  e=work;
  hx=e + Map.output_count;
  jacTe=hx + Map.output_count;
  jac=jacTe + Map.basis_count;
  jacTjac=jac + nm;
  Dp=jacTjac + Map.basis_count*Map.basis_count;
  diag_jacTjac=Dp + Map.basis_count;
  pDp=diag_jacTjac + Map.basis_count;

  /* see if starting point is within the feasile set */
  memmove(pDp,Map.weights,Map.basis_count*sizeof(T));
  ProjectOnBox<T>(Map.weights, lb, ub, Map.basis_count); /* project to feasible set */

  /* compute e=x - f(p) and its L2 norm */
  Map.func(UserData,u,Map.input_count,Map.weights,Map.basis_count,hx,Map.output_count); nfev=1;
  for(i=0, p_eL2=0.0; i<Map.output_count; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  };
  init_p_eL2=p_eL2;

  for(k=stop=0; k<itmax && !stop; ++k){
 //printf("%d  %.15g\n", k, 0.5*p_eL2);
    /* Note that p and e have been updated at a previous iteration */

    if(p_eL2<=epsf){ /* error is small */
      stop=6;
      break;
    };

    /* Compute the jacobian J at p,  J^T J,  J^T e,  ||J^T e||_inf and ||p||^2.
     * Since J^T J is symmetric, its computation can be speeded up by computing
     * only its upper triangular part and copying it to the lower part
     */

    Jacf(UserData,u,Map.input_count,Map.weights,Map.basis_count,jac,nm); ++njev;

    /* J^T J, J^T e */
    //if(nm<__BLOCKSZ__SQ){ // this is a small problem
      /* This is the straightforward way to compute J^T J, J^T e. However, due to
       * its noncontinuous memory access pattern, it incures many cache misses when
       * applied to large minimization problems (i.e. problems involving a large
       * number of free variables and measurements), in which J is too large to
       * fit in the L1 cache. For such problems, a cache-efficient blocking scheme
       * is preferable.
       *
       * Thanks to John Nitao of Lawrence Livermore Lab for pointing out this
       * performance problem.
       *
       * On the other hand, the straightforward algorithm is faster on small
       * problems since in this case it avoids the overheads of blocking.
       */

      for(i=0; i<Map.basis_count; ++i){
        for(j=i; j<Map.basis_count; ++j){

          for(l=0, tmp=0.0; l<Map.output_count; ++l){
            tmp+=jac[i*Map.output_count+l]*jac[j*Map.output_count+l];
          };

		      /* store tmp in the corresponding upper and lower part elements */
          jacTjac[i*Map.basis_count+j]=jacTjac[j*Map.basis_count+i]=tmp;
        };

        /* J^T e */
        for(l=0, tmp=0.0; l<Map.output_count; ++l)
          tmp+=jac[i*Map.output_count+l]*e[l];
        jacTe[i]=tmp;
      };
    //} else { // this is a large problem
    //  /* Cache efficient computation of J^T J based on blocking
    //   */
    //  TRANS_MAT_MAT_MULT(jac, jacTjac, Map.output_count, Map.basis_count);
    //
    //  /* cache efficient computation of J^T e */
    //  for(i=0; i<Map.basis_count; ++i)
    //    jacTe[i].q=0.0;
    //
    //  for(i=0; i<Map.output_count; ++i){
    //    T *jacrow;
    //
    //    for(l=0, jacrow=jac+i*Map.basis_count, tmp.q=e[i].q; l<Map.basis_count; ++l)
    //      jacTe[l].q+=jacrow[l].q*tmp.q;
    //  };
    //};

	  /* Compute ||J^T e||_inf and ||p||^2. Note that ||J^T e||_inf
     * is computed for free (i.e. inactive) variables only.
     * At a local minimum, if p[i]==ub[i] then g[i]>0;
     * if p[i]==lb[i] g[i]<0; otherwise g[i]=0
     */
    for(i=j=numactive=0, p_L2=jacTe_inf=0.0; i<Map.basis_count; ++i){
      if(ub && Map.weights[i]==ub[i]){ ++numactive; if(jacTe[i]>0.0) ++j; }
      else if(lb && Map.weights[i]==lb[i]){ ++numactive; if(jacTe[i]<0.0) ++j; }
      else if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;

      diag_jacTjac[i]=jacTjac[i*Map.basis_count+i]; /* save diagonal entries so that augmentation can be later canceled */
      p_L2+=Map.weights[i]*Map.weights[i];
    };
    //p_L2=sqrt(p_L2);

    // could add code for each iterations

    /* check for convergence */
    if(j==numactive && (jacTe_inf <= epsj)){
      Dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    };

   /* compute initial damping factor */
    if(k==0){
      if(!lb && !ub){ /* no bounds */
        for(i=0, tmp=FLT_MIN; i<Map.basis_count; ++i)
          if(diag_jacTjac[i]>tmp) tmp=diag_jacTjac[i]; /* find max diagonal element */
        mu=tau*tmp;
      } else
        mu=0.5*tau*p_eL2; /* use Kanzow's starting mu */
    };

    /* determine increment using a combination of adaptive damping, line search and projected gradient search */
    while(1){
      /* augment normal equations */
      for(i=0; i<Map.basis_count; ++i)
        jacTjac[i*Map.basis_count+i]+=mu;

      /* solve augmented equations */
      /* use the LU included with levmar */
      //issolved=AX_EQ_B_LU(jacTjac, jacTe, Dp, Map.basis_count);
      memmove(Dp,jacTe,Map.basis_count*sizeof(T));
      issolved=PLUDecomposition<T>(jacTjac,Map.basis_count,Dp,NULL);

      if(issolved){
        for(i=0; i<Map.basis_count; ++i)
          pDp[i]=Map.weights[i] + Dp[i];

        /* compute p's new estimate and ||Dp||^2 */
        ProjectOnBox<T>(pDp, lb, ub, Map.basis_count); /* project to feasible set */
        for(i=0, Dp_L2=0.0; i<Map.basis_count; ++i){
          Dp[i]=tmp=pDp[i]-Map.weights[i];
          Dp_L2+=tmp*tmp;
        };
        //Dp_L2=sqrt(Dp_L2);

        if(Dp_L2<=epsw_sq*p_L2){ /* relative change in p is small, stop */
          stop=2;
          break;
        };

        if(Dp_L2>=(p_L2+epsw)/(FLT_EPSILON*FLT_EPSILON)){ /* almost singular */
          stop=4;
          break;
        };


        Map.func(UserData,u,Map.input_count,pDp,Map.basis_count,hx,Map.output_count); ++nfev;
        for(i=0, pDp_eL2=0.0; i<Map.output_count; ++i){ /* compute ||e(pDp)||_2 */
          hx[i]=tmp=x[i]-hx[i];
          pDp_eL2+=tmp*tmp;
        };

        if(pDp_eL2<=gamma_sq*p_eL2){
          for(i=0, dL=0.0; i<Map.basis_count; ++i)
            dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

          if(dL>0.0){
            dF=p_eL2-pDp_eL2;
            tmp=(2.0*dF/dL-1.0);
            tmp=1.0-tmp*tmp*tmp;
            mu*=( (tmp>=ONE_THIRD)? tmp : ONE_THIRD );
          } else
            mu=(mu>=pDp_eL2)? pDp_eL2 : mu; /* pDp_eL2 is the new pDp_eL2 */

          nu=2;

          memmove(Map.weights,pDp,Map.basis_count*sizeof(T));

          memmove(e,hx,Map.output_count*sizeof(T));
          p_eL2=pDp_eL2;
          ++nLMsteps;
          gprevtaken=0;
          break;
        };
      } else {

      /* the augmented linear system could not be solved, increase mu */

        mu*=nu;
        nu2=nu<<1; // 2*nu;
        if(nu2<=nu){ /* nu has wrapped around (overflown). Thanks to Frank Jordan for spotting this case */
          stop=5;
          break;
        };
        nu=nu2;

        for(i=0; i<Map.basis_count; ++i) /* restore diagonal J^T J entries */
          jacTjac[i*Map.basis_count+i]=diag_jacTjac[i];

        continue; /* solve again with increased nu */
      };

      /* if this point is reached, the LM step did not reduce the error;
       * see if it is a descent direction
       */

      /* negate jacTe (i.e. g) & compute g^T * Dp */
      for(i=0, jacTeDp=0.0; i<Map.basis_count; ++i){
        jacTe[i]=-jacTe[i];
        jacTeDp+=jacTe[i]*Dp[i];
      };

      if(jacTeDp<=-rho*pow(Dp_L2, 1.05)){
        /* Dp is a descent direction; do a line search along it */
        int mxtake, iretcd;

        tmp=sqrt(p_L2);

        /* use the simpler (but slower!) line search described by Kanzow */
        for(t=tini; t>tmin; t*=beta){
          for(i=0; i<Map.basis_count; ++i) {
            pDp[i]=Map.weights[i] + t*Dp[i];
            //pDp[i]=__MEDIAN3(lb[i], pDp[i], ub[i]); /* project to feasible set */
          };

          Map.func(UserData,u,Map.input_count,pDp,Map.basis_count,hx,Map.output_count); ++nfev;
          for(i=0, pDp_eL2=0.0; i<Map.output_count; ++i){ /* compute ||e(pDp)||_2 */
            hx[i]=tmp=x[i]-hx[i];
            pDp_eL2+=tmp*tmp;
          };
          if(pDp_eL2<=p_eL2 + 2.0*t*alpha*jacTeDp) break;
        };

        ++nLSsteps;
        gprevtaken=0;

        /* NOTE: new estimate for p is in pDp, associated error in hx and its norm in pDp_eL2.
         * These values are used below to update their corresponding variables
         */
      } else {
gradproj: /* Note that this point can also be reached via a goto when LNSRCH() fails */

        /* jacTe is a descent direction; make a projected gradient step */

        /* if the previous step was along the gradient descent, try to use the t employed in that step */
        /* compute ||g|| */
        for(i=0, tmp=0.0; i<Map.basis_count; ++i)
          tmp=jacTe[i]*jacTe[i];
        tmp=100.0/(1.0+sqrt(tmp));
        t0=(tmp<=tini)? tmp : tini; /* guard against poor scaling & large steps; see (3.50) in C.T. Kelley's book */

        for(t=(gprevtaken)? t : t0; t>tming; t*=beta){
          for(i=0; i<Map.basis_count; ++i)
            pDp[i]=Map.weights[i] - t*jacTe[i] ;
          ProjectOnBox<T>(pDp, lb, ub, Map.basis_count); /* project to feasible set */
          for(i=0; i<Map.basis_count; ++i)
            Dp[i]=pDp[i]-Map.weights[i];

          Map.func(UserData,u,Map.input_count,pDp,Map.basis_count,hx,Map.output_count); ++nfev;
          for(i=0, pDp_eL2=0.0; i<Map.output_count; ++i){ /* compute ||e(pDp)||_2 */
            hx[i]=tmp=x[i]-hx[i];
            pDp_eL2+=tmp*tmp;
          };
          for(i=0, tmp=0.0; i<Map.basis_count; ++i) /* compute ||g^T * Dp|| */
            tmp+=jacTe[i]*Dp[i];

          if(gprevtaken && pDp_eL2<=p_eL2 + 2.0*0.99999*tmp){ /* starting t too small */
            t=t0;
            gprevtaken=0;
            continue;
          };
          if(pDp_eL2<=p_eL2 + 2.0*alpha*tmp) break;
        };

        ++nPGsteps;
        gprevtaken=1;
        /* NOTE: new estimate for p is in pDp, associated error in hx and its norm in pDp_eL2 */
      };

      /* update using computed values */

      for(i=0, Dp_L2=0.0; i<Map.basis_count; ++i){
        tmp=pDp[i]-Map.weights[i];
        Dp_L2+=tmp*tmp;
      };
      //Dp_L2=sqrt(Dp_L2);

      if(Dp_L2<=epsw_sq*p_L2){ /* relative change in p is small, stop */
        stop=2;
        break;
      };

      memmove(Map.weights,pDp,Map.basis_count*sizeof(T));

      memmove(e,hx,Map.output_count*sizeof(T));
      p_eL2=pDp_eL2;
      break;
    }; /* inner loop */
  };

  if(k>=itmax) stop=3;

  for(i=0; i<Map.basis_count; ++i) /* restore diagonal J^T J entries */
    jacTjac[i*Map.basis_count+i]=diag_jacTjac[i];

  if(info){
    info[0]=init_p_eL2;
    info[1]=p_eL2;
    info[2]=jacTe_inf;
    info[3]=Dp_L2;
    for(i=0, tmp=FLT_MIN; i<Map.basis_count; ++i)
      if(tmp<jacTjac[i*Map.basis_count+i]) tmp=jacTjac[i*Map.basis_count+i];
    info[4]=mu/tmp;
    info[5]=(T)k;
    info[6]=(T)stop;
    info[7]=(T)nfev;
    info[8]=(T)njev;
  };

  /* covariance matrix */
  if(covar){
    ComputeCovariance<T>(jacTjac, covar, p_eL2, Map.output_count, Map.basis_count);
  };

  delete[] work;

  // add code to return info
  return (stop!=4)?  k : -1;
};

template <class T>
void LBFGSmcstep(T* stx,
     T* fx,
     T* dx,
     T* sty,
     T* fy,
     T* dy,
     T* stp,
     T fp,
     T dp,
     bool* brackt,
     T stmin,
     T stmax,
     int* info)
{
  bool bound;
  T gamma;
  T p;
  T q;
  T r;
  T sgnd;
  T stpc;
  T stpf;
  T stpq;
  T theta;

  info[0] = 0;
  if( (brackt[0] && ((stp[0] <= MIN(stx[0],sty[0])) || (stp[0] >= MAX(stx[0],sty[0])))) || (dx[0]*(stp[0]-stx[0]) >= 0) || (stmax < stmin) )
    return;
  sgnd = dp*(dx[0]/fabs(dx[0]));
  if( fp > fx[0] ) {
    info[0] = 1;
    bound = true;
    theta = 3.0*(fx[0]-fp)/(stp[0]-stx[0])+dx[0]+dp;
    gamma = sqrt(theta*theta-dx[0]*dp);
    if( stp[0] < stx[0] )
      gamma = -gamma;
    p = gamma - dx[0] + theta;
    q = 2.0*gamma - dx[0] + dp;
    r = p / q;
    stpc = stx[0] + r*(stp[0] - stx[0]);
    stpq = stx[0] + dx[0]/((fx[0] - fp)/(stp[0] - stx[0]) + dx[0])/2.0*(stp[0] - stx[0]);
    if( fabs(stpc-stx[0])<fabs(stpq-stx[0]) )
      stpf = stpc;
    else
      stpf = stpc+(stpq-stpc)/2.0;
    brackt[0] = true;
  } else {
    if( sgnd < 0 ) {
      info[0] = 2;
      bound = false;
      theta = 3.0*(fx[0]-fp)/(stp[0]-stx[0])+dx[0]+dp;
      gamma = sqrt(theta*theta-dx[0]*dp);
      if( stp[0] > stx[0] )
        gamma = -gamma;
      p = gamma-dp+theta;
      q = 2.0*gamma-dp+dx[0];
      r = p/q;
      stpc = stp[0]+r*(stx[0]-stp[0]);
      stpq = stp[0]+dp/(dp-dx[0])*(stx[0]-stp[0]);
      if( fabs(stpc-stp[0])>fabs(stpq-stp[0]) )
        stpf = stpc;
      else
        stpf = stpq;
      brackt[0] = true;
    } else {
      if( fabs(dp) < fabs(dx[0]) ) {
        info[0] = 3;
        bound = true;
        theta = 3.0*(fx[0]-fp)/(stp[0]-stx[0])+dx[0]+dp;
        gamma = sqrt(MAX(0.0, theta*theta - dx[0]*dp));
        if( stp[0] > stx[0] )
          gamma = -gamma;
        p = gamma - dp + theta;
        q = 2.0*gamma + dx[0] - dp;
        r = p/q;
        if( (r < 0) && (gamma != 0) )
          stpc = stp[0] + r*(stx[0] - stp[0]);
        else {
          if( stp[0] > stx[0] )
            stpc = stmax;
          else
            stpc = stmin;
        };
        stpq = stp[0] + dp/(dp-dx[0])*(stx[0]-stp[0]);
        if( brackt[0] ) {
          if( fabs(stp[0]-stpc)<fabs(stp[0]-stpq) )
            stpf = stpc;
          else
            stpf = stpq;
        } else {
          if( fabs(stp[0]-stpc)>fabs(stp[0]-stpq) )
            stpf = stpc;
          else
            stpf = stpq;
        };
      } else {
        info[0] = 4;
        bound = false;
        if( brackt[0] ) {
          theta = 3.0*(fp-fy[0])/(sty[0]-stp[0])+dy[0]+dp;
          gamma = sqrt(theta*theta-dy[0]*dp);
          if( stp[0] > sty[0] )
            gamma = -gamma;
          p = gamma-dp+theta;
          q = 2.0*gamma-dp+dy[0];
          r = p/q;
          stpf = stpc = stp[0]+r*(sty[0]-stp[0]);
        } else {
          if( stp[0] > stx[0] )
            stpf = stmax;
          else
            stpf = stmin;
        };
      };
    };
  };
  if( fp > fx[0] ) {
    sty[0] = stp[0];
    fy[0] = fp;
    dy[0] = dp;
  } else {
    if( sgnd < 0.0 ) {
      sty[0] = stx[0];
      fy[0] = fx[0];
      dy[0] = dx[0];
    };
    stx[0] = stp[0];
    fx[0] = fp;
    dx[0] = dp;
  };
  stpf = MAX(stmin, MIN(stmax, stpf));
  stp[0] = stpf;
  if( brackt[0] && bound ) {
    if( sty[0] > stx[0] )
      stp[0] = MIN(stx[0]+0.66*(sty[0]-stx[0]), stp[0]);
    else
      stp[0] = MAX(stx[0]+0.66*(sty[0]-stx[0]), stp[0]);
  };
};

template <class T>
int LBFGSmcsrch(void* UserData,
     TBasisMap Map,
     FBasisFunction Grad,
     T* x,
     T* f,
     T* g,
     T* s,
     T* stp,
     T ftol,
     T xtol,
     int maxfev,
     int* nfev,
     T* wa,
     T gtol,
     T stpmin,
     T stpmax)
{
  int infoc;
  int j;
  bool brackt;
  bool stage1;
  T dg;
  T dgm;
  T dginit;
  T dgtest;
  T dgx;
  T dgxm;
  T dgy;
  T dgym;
  T finit;
  T ftest1;
  T fm;
  T fx;
  T fxm;
  T fy;
  T fym;
  T p5;
  T p66;
  T stx;
  T sty;
  T stmin;
  T stmax;
  T width;
  T width1;
  T xtrapf;
  T mytemp;

  p5 = 0.5;
  p66 = 0.66;
  xtrapf = 4.0;
  Map.func(UserData,x,Map.input_count,Map.weights,Map.basis_count,f,1);
  Grad(UserData,x,Map.input_count,Map.weights,Map.basis_count,g,Map.input_count);
  infoc = 1;
  dginit = DotProduct<T>(Map.input_count,g,s);
  if( dginit>=0 )
    return 0;
  brackt = false;
  stage1 = true;
  nfev[0] = 0;
  finit = f[0];
  dgtest = ftol*dginit;
  width = stpmax-stpmin;
  width1 = 2.0*width;
  memmove(wa,x,Map.input_count*sizeof(T));
  stx = 0.0;
  fx = finit;
  dgx = dginit;
  sty = 0.0;
  fy = finit;
  dgy = dginit;
  while(true) {
    if( brackt ) {
      if( stx < sty ) {
        stmin = stx;
        stmax = sty;
      } else {
        stmin = sty;
        stmax = stx;
      };
    } else  {
      stmin = stx;
      stmax = stp[0]+xtrapf*(stp[0]-stx);
    };
    if( stp[0] > stpmax )
      stp[0] = stpmax;
    if( stp[0] < stpmin )
      stp[0] = stpmin;
    if( ((brackt) && ((stp[0] <= stmin) || (stp[0] >= stmax) || (stmax-stmin <= xtol*stmax))) || (nfev[0] >= maxfev-1) || (infoc==0) )
      stp[0] = stx;
    for(j=0;j<Map.input_count;j++)
      x[j] = wa[j]+stp[0]*s[j];
    Map.func(UserData,x,Map.input_count,Map.weights,Map.basis_count,f,1);
    Grad(UserData,x,Map.input_count,Map.weights,Map.basis_count,g,Map.input_count);
    nfev[0]++;
    dg = 0.0;
    for(j=0;j<Map.input_count;j++)
      dg += g[j]*s[j];
    ftest1 = finit+stp[0]*dgtest;
    if( (brackt && ((stp[0] <= stmin) || (stp[0] >= stmax))) || (infoc == 0) )
      return 6;
    if( (stp[0] == stpmax) && (f[0] <= ftest1) && (dg <= dgtest) )
      return 5;
    if( (stp[0] == stpmin) && ((f[0] > ftest1) || (dg >= dgtest)) )
      return 4;
    if( nfev[0] >= maxfev )
      return 3;
    if( brackt && (stmax-stmin <= xtol*stmax) )
      return 2;
    if( (f[0] <= ftest1) && (fabs(dg) <= -gtol*dginit) )
      return 1;
    mytemp = ftol;
    if( gtol < ftol )
      mytemp = gtol;
    if( stage1 && (f[0] <= ftest1) && (dg >= mytemp*dginit) )
      stage1 = false;
    if( stage1 && (f[0] <= fx) && (f[0] > ftest1) ) {
      fm = f[0]-stp[0]*dgtest;
      fxm = fx-stx*dgtest;
      fym = fy-sty*dgtest;
      dgm = dg-dgtest;
      dgxm = dgx-dgtest;
      dgym = dgy-dgtest;
      LBFGSmcstep<T>(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, fm, dgm, &brackt, stmin, stmax, &infoc);
      fx = fxm+stx*dgtest;
      fy = fym+sty*dgtest;
      dgx = dgxm+dgtest;
      dgy = dgym+dgtest;
    } else
      LBFGSmcstep<T>(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, f[0], dg, &brackt, stmin, stmax, &infoc);
    if( brackt ) {
      if( fabs(sty-stx)>=p66*width1 )
        stp[0] = stx+p5*(sty-stx);
      width1 = width;
      width = fabs(sty-stx);
    };
  };
};

template <class T>
int LBFGSMinimize(void* UserData,
     TBasisMap Map,
     FBasisFunction Grad,
     unsigned int M,
     T* x,
     T epsg,
     T epsf,
     T epsx,
     unsigned int maxits)
{
  if((Map.output_count > 1) || (maxits < 2) || (Grad == NULL) || (epsg < 0.0) || (epsf < 0.0) || (epsx < 0.0)) return OPT_ERR_IMPROPER;
  if(M > Map.input_count) M = Map.input_count;

  T* w = (T*) new T[Map.input_count*(2*M+1)+2*M];
  T f;
  T fold;
  T tf;
  T v;
  T* xold = (T*) new T[Map.input_count];
  T* tx = (T*) new T[Map.input_count];
  T* g = (T*) new T[Map.input_count];
  T* diag = (T*) new T[Map.input_count];
  T* ta = (T*) new T[Map.input_count];
  T gnorm;
  T stp1;
  T ftol;
  T stp;
  T ys;
  T yy;
  T sq;
  T yr;
  T beta;
  int iter;
  int nfun;
  int point;
  int ispt;
  int iypt;
  int maxfev;
  int bound;
  int npt;
  int cp;
  int nfev;
  unsigned int i;
  int inmc;
  int iycn;
  int iscn;
  T xtol;
  T gtol;
  T stpmin;
  T stpmax;

  Map.func(UserData,x,Map.input_count,Map.weights,Map.basis_count,&f,1);
  Grad(UserData,x,Map.input_count,Map.weights,Map.basis_count,g,Map.input_count);
  fold = f;
  iter = 0;

  nfun = 1;
  point = 0;
  npt = 0;
  for(i=0;i<Map.input_count;i++)
    diag[i] = 1.0;
  xtol = 0.0000001; // for single precision Ting point value
  gtol = 0.9;
  stpmin = pow(10, -20.0);
  stpmax = pow(10, 20.0);
  ispt = Map.input_count+2*M;
  iypt = ispt+Map.input_count*M;
  for(i=0;i<Map.input_count;i++)
    w[ispt+i] = -g[i]*diag[i];
  gnorm = sqrt(DotProduct<T>(Map.input_count, g, g));
  stp1 = 1.0/gnorm;
  ftol = 0.0001;
  maxfev = 20;
  while(true)
  {
    memmove(xold,x,Map.input_count*sizeof(T));
    bound = ++iter;
    if( iter!=1 )
    {
      if(iter>M)
        bound = M;
      ys = DotProduct<T>(Map.input_count,&w[iypt+npt],&w[ispt+npt]);
      yy = DotProduct<T>(Map.input_count,&w[iypt+npt],&w[iypt+npt]);
      for(i=0;i<Map.input_count;i++)
        diag[i] = ys/yy;
      cp = point;
      if( point<0 )
        cp = M-1;
      w[Map.input_count+cp] = 1.0/ys;
      for(i=0;i<Map.input_count;i++)
        w[i] = -g[i];
      cp = point;
      for(cp=point-1;cp>=0;cp--)
      {
        sq = DotProduct<T>(Map.input_count,&w[ispt+cp*Map.input_count],w);
        inmc = Map.input_count+M+cp;
        iycn = iypt+cp*Map.input_count;
        w[inmc] = w[Map.input_count+cp]*sq;
        VectLinComb<T>(Map.input_count, w, 1.0, &w[iycn], -w[inmc]);
      };
      for(i=0;i<Map.input_count;i++)
        w[i] *= diag[i];
      for(cp=0;cp<bound;cp++)
      {
        yr = DotProduct<T>(Map.input_count,&w[iypt+cp*Map.input_count],w);
        beta = w[Map.input_count+cp]*yr;
        inmc = Map.input_count+M+cp;
        beta = w[inmc]-beta;
        iscn = ispt+cp*Map.input_count;
        VectLinComb<T>(Map.input_count, w, 1.0, &w[iscn], beta);
      };
      for(i=0;i<Map.input_count;i++)
        w[ispt+point*Map.input_count+i] = w[i];
    }
    nfev = 0;
    stp = 1.0;
    if( iter==1 )
      stp = stp1;
    memmove(w,g,Map.input_count*sizeof(T));
    int info = LBFGSmcsrch<T>(UserData,Map,Grad, x, &f, g, &w[ispt+point*Map.input_count], &stp, ftol, xtol, maxfev, &nfev, diag, gtol, stpmin, stpmax);
    if( info!=1 )
      if( info==0 )
        return OPT_ERR_IMPROPER;
    nfun = nfun+nfev;
    npt = point*Map.input_count;
    for(i=0;i<Map.input_count;i++)
    {
      w[ispt+npt+i] = stp*w[ispt+npt+i];
      w[iypt+npt+i] = g[i]-w[i];
    };
    point++;
    if( point==M )
      point = 0;
    if(iter>maxits) {
      delete[] w;
      delete[] xold;
      delete[] tx;
      delete[] g;
      delete[] diag;
      delete[] ta;
      return OPT_MAX_ITERATION;
    };

    //could add code here for reporting each iterations.

    gnorm = sqrt(DotProduct<T>(Map.input_count,g, g));
    tf = MAX(fabs(fold), MAX(fabs(f), 1.0));
    memmove(tx,xold,Map.input_count*sizeof(T));
    VectSub<T>(Map.input_count,tx, x);
    v = sqrt(DotProduct<T>(Map.input_count, tx, tx));
    if(( v<=epsx ) || ( gnorm<=epsg ) || ( fold-f<=epsf*tf )) {
      delete[] w;
      delete[] xold;
      delete[] tx;
      delete[] g;
      delete[] diag;
      delete[] ta;
      return OPT_OPTIMAL;
    };

    fold = f;
    memmove(xold, x,Map.input_count*sizeof(T));
  };
};

#endif //RK_OptimizationC


//---------------------------------------------------------------------------
#pragma package(smart_init)

