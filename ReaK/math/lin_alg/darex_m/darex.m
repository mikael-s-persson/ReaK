function [A,B,Q,R,S,X,parout,G,C,Q0]=darex(index,parin)
%DAREX
%
% Test examples for the discrete-time algebraic Riccati equation (DARE) 
%
%                                                   -1
% (I)   0 = DR(X) = A'XA - X - (A'XB + S) (R + B'XB)  (B'XA + S') + Q
%
% Here, A,Q, and X are n-by-n matrices, B and S are n-by-m, and R is 
% m-by-m. Q and R are symmetric and X is the required solution matrix. 
% One common approach to solve DAREs is to compute a deflating subspace
% of the symplectic pencil 
%
%                           -1
%                  ( A - B R  S    0 )     ( I        G        )
% (II)  L - s M := (    -1           ) - s (             -1    )
%                  ( S R  S' - Q   I )     ( 0   (A - B R  S')')
%
%               -1  
% where  G = B R  B' is a symmetric n-by-n matrix. Q may also be given 
% in factored form, Q = C' Q0 C, where C is a p-by-n and Q0 is a p-by-p
% matrix. 
% NOTE that for DAREs, R being a singular matrix is not uncommon. In this 
% case, the symplectic pencil cannot be formed as in (II), but a solution 
% of the DARE can be computed via a deflating subspace of the extended
% pencil 
%
%                    ( A   0   B )     ( I   0   0 )
% (III) LL - s MM := ( Q  -I   S ) - s ( 0  -A'  0 ) . 
%                    ( S'  0   R )     ( 0  -B'  0 )     
%
% For examples with singular R-matrix, G can not be computed and is thus
% not returned.  
%
% Input:
%  - index: number of example to generate, indices refer to example 
%           numbers in [1].
%  - parin: input parameters (optional, defaults values given in [1]). 
%           For Example number 
%           +  1-11: not referenced ([1], Section 2).
%           + 12-13: parin(1) = real-valued scalar.
%           + 14   : parin(1:4) = [tau, D, K, r], real-valued scalars.
%           + 15   : parin(1) = n = problem size.
%                    parin(2) = r = real-valued scalar.
%
% Output:
%  - A, B, Q, R, S: coefficient matrices of DARE as in (I).
%  - X            : exact solution of DARE (if available), usually the
%                   stabilizing solution.  If an exact solution is not
%                   available, the empty matrix is returned.
%  - parout       : vector with system properties,
%                   parout(1:3) = [n, m, p].
%                   parout(4)   = 2-norm condition number of A.
%                   The following parameters are only returned if an 
%                   solution of the DARE is available:
%                   parout(5)   = absolute value of closed-loop eigenvalue
%                                 of largest modulus. 
%                   parout(6)   = 2-norm of X.
%                   parout(7)   = 2-norm condition number of X.
%                   parout(8)   = condition number of DARE as defined in [2].
%  - G, C, Q0     : optional output matrices as defined above. NOTE that
%                   G can only be computed if R is nonsingular. Otherwise,
%                   G contains on output the empty matrix.
%
% References:
%
% [1] P.BENNER, A.J. LAUB, V. MEHRMANN: 'A Collection of Benchmark 
%     Examples for the Numerical Solution of Algebraic Riccati 
%     Equations II: Discrete-Time Case', Tech. Report SPC 95_23, 
%     Fak. f. Mathematik, TU Chemnitz-Zwickau (Germany), December 1995.
% [2] T. GUDMUNDSSON, C. KENNEY, A.J. LAUB: 'Scaling of the Discrete-Time
%     Algebraic Riccati Equation to Enhance Stability of the Schur 
%     Solution Method', IEEE Transactions on Automatic Control, vol. 37, 
%     no. 4, pp. 513-518, 1992. 

%  Peter Benner, Volker Mehrmann (TU Chemnitz-Zwickau, Germany),
%  Alan J. Laub (University of California at Santa Barbara)
%  12-14-1995, 02-28-1996 

nex   = 15;
nfsex = 14;

error(nargchk(1,2,nargin))
if (nargin > 1  &  index <= nex) 
  if (index > nfsex & index < nex) & (abs(round(parin(1))) ~= parin(1)),
    msg = sprintf('Parameter for example #%i must be positive integer.',index);
    error(msg)
  end
end    

X = [];
S = [];

if index == 1,
  A = [4, 3; -4.5, -3.5];
  B = [1; -1];  
  R = 1;
  Q = [9, 6; 6, 4];
  X = (1+sqrt(5))/2*[9, 6; 6, 4];
  parout = [2, 1, 2];
  if nargout > 6
      G = [1, -1; -1, 1];  C = eye(2);  Q0 = Q;
  end

elseif index == 2,
  A = [0.9512, 0; 0, 0.9048];
  B = [4.877, 4.877; -1.1895, 3.569];
  R = [1/3, 0; 0, 3];
  Q = [0.005, 0; 0, 0.02];
  parout = [2, 2, 2];
  if nargout > 6,
      G = B/R*B';  C = eye(2);  Q0 = Q;
  end

elseif index == 3,
  A = [2, -1; 1, 0];
  B = [1; 0];
  R = 0;
  Q = [0, 0; 0, 1];
  X = eye(2);
  parout = [2, 1, 1];
  if nargout > 6,
    G = [];  C = [0 1];  Q0 = 1;
  end

elseif index == 4,
  A = [0, 1; 0, -1];
  B = [1, 0; 2, 1];
  R = [9, 3; 3, 1];
  d = -4/11;
  Q = [d, d; d, 7/11];
  S = [3, 1; -1,  7];
  parout = [2, 2, 2];
  if nargout > 6,
    G = [];  C = eye(2);  Q0 = Q;
  end

elseif index == 5,
  A = [0, 1; 0, 0];
  B = [0; 1];
  R = 1;
  Q = [1, 2; 2, 4];
  X = [1, 2; 2, 2+sqrt(5)]; 
  parout = [2, 1, 2];
  if nargout > 6,
    G = B*B';  C = eye(2);  Q0 = Q;
  end

elseif index == 6,
  A = 0.998*eye(4);
  A(1,2) = 0.067;   A(2,1) = -A(1,2);
  A(3,4) = 0.153;   A(4,3) = -A(3,4);
  B = [0.0033, 0.02; 0.1, -0.0007; 0.04, 0.0073; -0.0028, 0.1];
  R = eye(2);
  Q = [1.87, 0, 0, -0.244; 0, 0.744, 0.205, 0; 0, 0.205, 0.589, 0;...
       -0.244, 0, 0, 1.048];
  parout = [4, 2, 4];
  if nargout > 6,
    G = B*B';  C = eye(4);  Q0 = Q;
  end

elseif index == 7,
  A = [ 0.98475,  -0.079903,  0.0009054, -0.0010765;... 
        0.041588,  0.99899,  -0.035855,   0.012684;...
       -0.54662,   0.044916, -0.32991,    0.19318;...
        2.6624,   -0.10045,  -0.92455,   -0.26325];
  B = [0.0037112, 0.0007361; -0.087051, 9.3411e-6; -1.19844, -4.1378e-4;...
       -3.1927, 9.2535e-4];
  R = eye(2);
  Q = 0.01*eye(4);
  parout = [4, 2, 4];
  if nargout > 6,
    G = B*B';  C = eye(4);  Q0 = Q;
  end

elseif index == 8,
  B  = -triu(ones(4)) + diag(2*ones(4,1));
  R  = eye(4);
  C  = triu(ones(4)) + diag(ones(2,1),2) + diag(3,3);
  Q0 = [2, -1, 0, 0; -1, 2, -1, 0; 0, -1, 2, 0; 0, 0, 0, 0];
  Q  = C'*Q0*C;
  A  = B*[0.4, 0, 0, 0; 1, 0.6, 0, 0; 0, 1, 0.8, 0; 0, 0, 0, -0.999982]*C;
  parout = [4, 4, 4];
  if nargout > 6,
    G = B*B';  
  end

elseif index == 9,
  A = 0.01*[ 95.407,   1.9643,  0.3597,  0.0673,  0.019;... 
             40.849,  41.317,  16.084,   4.4679,  1.1971;...
             12.217,  26.326,  36.149,  15.93,   12.383;...
              4.1118, 12.858,  27.209,  21.442,  40.976;...
              0.1305,  0.5808,  1.875,   3.6162, 94.28];
  B = 0.01*[ 0.0434,  2.6606,  3.753,  3.6076,  0.4617;...
            -0.0122, -1.0453, -5.51,  -6.6,    -0.9148]';
  R = eye(2);  
  Q = eye(5);
  parout = [5, 2, 5];
  if nargout > 6,
    G = B*B';  C = eye(5);  Q0 = Q;
  end

elseif index == 10,
  A = kron(eye(2), diag([1,1], 1));
  B = [0, 0, 1, 0, 0, 0;  0, 0, 0, 0, 0, 1]'; 
  C = [1, 1, 0, 0, 0, 0;  0, 0, 0, 1, -1, 0];
  D = [1, 0;  1, 0];
  Q = C'*C;
  R = eye(2) + D'*D;
  S = C'*D;
  parout = [6, 2, 2];
  if nargout > 6,  
    G = B/R*B';  
    Q0 = eye(2);  
  end
  
elseif index == 11,
  A  = [870.1,  135.0   11.59   0.5014 -37.22  0.3484 0.0     4.242  7.249;...
         76.55 897.4   12.72   0.5504 -40.16  0.3743 0.0     4.53   7.499;...
       -127.2  357.5  817.0    1.455 -102.8   0.987  0.0    11.85  18.72;...
       -363.5  633.9   74.91 796.6   -273.5   2.653  0.0    31.72  48.82;...
       -960.0 1645.9 -128.9   -5.597   71.42  7.108  0.0    84.52 125.9;...
       -664.4  112.96 -88.89  -3.854   84.47 13.6    0.0   144.3  101.6;...
       -410.2  693.0  -54.71  -2.371   66.49 12.49   0.1063 99.97  69.67;...
       -179.9  301.7  -23.93  -1.035   60.59 22.16   0.0   213.9   35.54;...
       -345.1  580.4  -45.96  -1.989  105.6  19.86   0.0   219.1  215.2 ];
  A  = 0.001*A;
  B  = [  4.76    0.879   1.482  3.892 10.34   7.203  4.454  1.971  3.773;...
         -0.5701 -4.773 -13.12 -35.13 -92.75 -61.59 -36.83 -15.54 -30.28;...
        -83.68   -2.73    8.876 24.8   66.8   38.34  20.29   6.937 14.69 ]';
  B  = 0.0001*B;
  R  = eye(3);
  C  = [1, zeros(1,8); zeros(1,4), 1, zeros(1,4)];
  Q0 = 50*eye(2);
  Q = C'*Q0*C;
  parout = [9, 3, 2];
  if nargout > 6,  G = B*B';  end

elseif index == 12,
  if nargin < 2,   d = 10^6;   else,   d = parin(1);   end
  A = [0, d; 0, 0];
  B = [0; 1];  
  R = 1;
  Q = eye(2);
  X = [1, 0; 0, 1+d^2];
  parout = [2, 1, 2];
  if nargout > 6,
    G = [0, 0; 0, 1];  C =eye(2);  Q0 = Q;
  end

elseif index == 13;
  if nargin < 2,   d = 10^6;   else,   d = parin(1);   end
  U = eye(3) - (2/3)*ones(3);
  A = U*diag([0, 1, 3])*U;
  B = eye(3);   
  R = d*eye(3);
  Q = d*eye(3); 
  X = d*U*diag([1, (1+sqrt(5))/2, (9+sqrt(85))/2])*U;
  parout = [3, 3, 3];
  if nargout > 6,
    G = eye(3)/d;  C = U;  Q0 = d*eye(3);
  end

elseif index == 14,
  if nargin < 2,  lp = 0;  else,  lp = length(parin);  end
  if lp < 4,  R   = 0.25;  else,  R   = parin(4);  end
  if lp < 3,  K   = 1;     else,  K   = parin(3);  end
  if lp < 2,  D   = 1;     else,  D   = parin(2);  end  
  if lp < 1,  tau = 10^8;  else,  tau = parin(1);  end
  alpha = 1 - D/tau;
  beta  = K*D/tau;
  A = diag(ones(3,1),-1);  A(1,1) = alpha;
  B = [beta; 0; 0; 0];
  Q = diag(flipud(eye(4,1)));
  d  = (alpha + 1) * (alpha - 1);
  dd = (R*tau/K/D)*tau/K/D;
  X = eye(4);
  X(1,1) = dd*d/2 + 0.5 + sqrt(0.25*(dd*d + 1)^2 + dd);
  parout = [4, 1, 1];
  if nargout > 6,
    G = diag([beta*beta/R, 0, 0 ,0]);  C = [0, 0, 0, 1];  Q0 = 1;
  end

elseif index == 15,
  if nargin < 2,  lp = 0;  else,  lp = length(parin);  end
  if lp < 2,  R = 1;     else,  R = parin(2);  end
  if lp < 1,  n = 100;   else,  n = parin(1);  end
  A = diag(ones(n-1,1),1);
  B = flipud(eye(n,1));
  Q = eye(n);
  X = diag(1:n);
  parout = [n, 1, n];
  if nargout > 6,
    G = B*B';  C = eye(n);  Q0 = Q;
  end

else
  error('This example is not available!');
end  