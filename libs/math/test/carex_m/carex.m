function [A,G,Q,X,parout,B,R,C,Q0]=carex(index,parin)
%CAREX
%
% Test examples for the continuous-time algebraic Riccati equation    
%
%  (CARE)   0 = Q  +  A' X  +  X A  -  X G X.
%
% Here, A, G, Q, X are n-by-n matrices,  G and Q are symmetric.  
% G, Q may be in factored form  G = B R^(-1) B',   Q = C' Q0 C.
% Then  B is n-by-m, R m-by-m, C p-by-n, and Q0 is p-by-p. The
% corresponding Hamiltonian matrix is defined as
%
%                          ( A  -G )   (    A    -B/R B')
%      H := Ham(A,G,Q) := (       ) = (                       ).
%                          (-Q  -A')   (-C' Q0 C    -A' )
%
% Input:
%  - index : number of example to generate, indices refer to example 
%            numbers in [1].
%  - parin : input parameters (optional, default values given in [1]). 
%      For Example number 
%      +  1- 6: not referenced ([1], Section 2).
%      +  7-14: parin(1) = real-valued scalar, ([1], Section 3).
%      + 15   : parin(1) = N,  n = 2*N-1,  m = N,  p = N-1.
%      + 16-18: parin(1) = n = problem size = size of A, G, Q, X.
%      + 17   : parin(2) = Q0 (real scalar).  
%               parin(3) = R  (real scalar).
%      + 18   : parin(2:8) = real-valued scalars, where
%                parin(2:4) = [a, b, c].
%               parin(5:6) = [beta_1,beta_2].
%                parin(7:8) = [gamma_1,gamma_2].
%      + 19   : parin(1) = l = number of springs (n = 2l).
%               parin(2:4) = [mu, kappa, delta].
%      + 20   : parin = name of data file containing the number of masses
%                       l (n = 2*l-1), and the vectors mu, delta, gamma,
%                       kappa described in [1].
%                         
% Output:
%  - A, G, Q : system matrices from CARE.
%  - X             : exact stabilizing solution of CARE (if available).
%              If an exact (analytical) solution is not available, the
%              empty matrix is returned.
%              For Example 17,  X = X(1,n) = X(n,1) = sqrt(Q0*R), which
%              is the only available information.  
%  - parout  : vector with system properties,
%              parout(1:3) = [n, m, p]
%               parout(4)   = norm(H) = 2-norm of H = Ham(A,G,Q)
%               parout(5)   = 2-norm condition number of H.
%              The following parameters are only returned if an analytical 
%              solution of the CARE is available:
%               parout(6)   = 2-norm of X
%               parout(7)   = 2-norm condition number of X
%               parout(8:9) = [KU,KL] = upper and lower bound for condition 
%                             number of CARE (see [2]) (not available for
%                            Example 11).  
%                   
%  - B,R,C,Q0: optional output matrices if factored form is required.
% 
% References:
% [1] P.BENNER, A.J. LAUB, V. MEHRMANN: 'A Collection of Benchmark 
%     Examples for the Numerical Solution of Algebraic Riccati 
%     Equations I: Continuous-Time Case', Tech. Report SPC 95_23, 
%     Fak. f. Mathematik, TU Chemnitz-Zwickau (Germany), October 1995.
% [2] C. KENNEY, G. HEWER: 'The sensitivity of the algebraic and 
%     differential Riccati equations', SIAM J. Control Optim., 
%     vol. 28 (1990), pp.50-69.
 
%  Peter Benner, Volker Mehrmann (TU Chemnitz-Zwickau, Germany),
%  Alan J. Laub (University of California at Santa Barbara)
%  10-09-1995 
%
%  For questions concerning this M-file, send e-mail to
%
%        benner@mathematik.tu-chemnitz.de

nex   = 20;
nfsex = 14;

if (nargin > 1  &  index <= nex) 
  if (index > nfsex & index < nex) & (abs(round(parin(1))) ~= parin(1)),
    msg = sprintf('Parameter for example #%i must be positive integer.',index);
    error(msg)
  end
end    

X = [];

if index == 1,
  A = [0 1; 0 0];
  G = [0 0; 0 1];
  Q = [1 0; 0 2];
  X = [2 1; 1 2];
  parout = [2, 1, 2];
  if nargout > 5
    B = [0; 1];  R = 1;  C = eye(2);  Q0 = Q;
  end

elseif index == 2,
  A = [4 3; -4.5 -3.5];
  G = [1 -1; -1 1];
  Q = [9 6; 6 4];
  X = (1+sqrt(2))*[9 6; 6 4];
  parout = [2, 1, 2];
  if nargout > 5
    B = [1; -1];  R = 1;  C = eye(2);  Q0 = Q;
  end

elseif index == 3,
  A = [0      1       0     0; ...
       0     -1.89    0.39 -5.53; ... 
       0     -0.034  -2.98  2.43; ...
       0.034 -0.0011 -0.99 -0.21];
  B = [0 0; 0.36 -1.6; -0.95 -0.032; 0.03 0];
  G = B*B';
  Q = [2.313 2.727 0.688 0.023; ...
       2.727 4.271 1.148 0.323; ...
       0.688 1.148 0.313 0.102; ...
       0.023 0.323 0.102 0.083];
  parout = [4, 2, 4];
  if nargout > 5,
    R = eye(2);  C = eye(4);  Q0 = Q;
  end

elseif index == 4,
  A = [-0.991   0.529   0       0       0       0       0       0     ; ...
        0.522  -1.051   0.596   0       0       0       0       0     ; ...
        0       0.522  -1.118   0.596   0       0       0       0     ; ...
        0       0       0.522  -1.548   0.718   0       0       0     ; ...
        0       0       0       0.922  -1.640   0.799   0       0     ; ...
        0       0       0       0       0.922  -1.721   0.901   0     ; ...
        0       0       0       0       0       0.922  -1.823   1.021 ; ...
        0       0       0       0       0       0       0.922  -1.943];
  B = 0.001*[ 3.84  4.00 37.60  3.08  2.36  2.88  3.08  3.00; ...
             -2.88 -3.04 -2.80 -2.32 -3.32 -3.82 -4.12 -3.96]';
  G = B*B';
  Q = [0.5  0  0  0.1; 0.1  0  0  0; 0  0.5  0  0; 0  0  0  0];
  Q = [eye(4), Q; Q', 0.1*eye(4)];
  parout = [8, 2, 8];
  if nargout > 5,
    R = eye(2);  C = eye(8);  Q0 = Q;
  end

elseif index == 5,
  A = [-4.019  5.12   0      0    -2.082     0       0     0    0.87; ...
       -0.346  0.986  0      0    -2.34      0       0     0    0.97; ...
       -7.909 15.407 -4.069  0    -6.45      0       0     0    2.68; ...
      -21.816 35.606 -0.339 -3.87 -17.8      0       0     0    7.39; ...
      -60.196 98.188 -7.907  0.34 -53.008    0       0     0    20.4; ...
        0      0      0      0     94.0   -147.2     0    53.2   0; ...
        0      0      0      0      0       94.0  -147.2   0     0; ...
        0      0      0      0      0       12.8     0   -31.6   0; ...
        0      0      0      0     12.8      0       0    18.8 -31.6];
  B = [ 0.010  0.003  0.009  0.024  0.068 0 0 0 0; ...
       -0.011 -0.021 -0.059 -0.162 -0.445 0 0 0 0; ...
       -0.151  0      0      0      0     0 0 0 0]';
  G = B*B';
  Q = eye(9);
  parout = [9, 3, 9];
  if nargout > 5,
    R = eye(3);  C = eye(9);  Q0 = Q;
  end

elseif index == 6,
  load carex6
  A = [A1, A121, zeros(16,1), A122, zeros(16,2), A123, zeros(16,8);...
       zeros(2,16), A21, zeros(2,12); ...
       zeros(3,18), A22, zeros(3,9); ...
       zeros(3,21), A23, zeros(3,6); ...
       A31, zeros(6,8), A3];
  B = [zeros(16,3);  B21, zeros(2,2);  zeros(3,1), B22, zeros(3,1); ...
       zeros(3,2), B23;  zeros(6,3)];
  G = B*B';
  Q = C'*C;
  parout = [30, 3, 5];
  if nargout > 5,
    R = eye(3);   Q0 = eye(5);
  end

elseif index == 7,
  if nargin < 2,   d = 10^(-6);   else,   d = parin(1);   end
  dd = d^2;
  A = [1 0; 0 -2];
  G = [dd 0; 0 0];
  Q = ones(2);
  s  = sqrt(1 + dd);
  s2 = 2 + s;
  X = [(1+s)/dd  1/s2; 1/s2  (1-dd/(s2^2))/4];
  parout = [2, 1, 1];
  if nargout > 5,
    B = [d; 0];  R = 1;  C = [1 1];  Q0 = 1;
  end

elseif index == 8,
  if nargin < 2,   d = 10^(-8);   else,   d = parin(1);   end
  A = [-0.1 0; 0 -0.02];
  B = [0.1 0; 0.001 0.01];
  R = [1+d 1; 1 1];
  G = B/R*B';
  Q = [100 1000; 1000 10000];
  parout = [2, 2, 1];
  if nargout > 5,  
    C = [10 100];   Q0 = 1;
  end

elseif index == 9,
  if nargin < 2,   d = 10^6;   else,   d = parin(1);   end
  A = [0 d; 0 0];
  G = [0 0; 0 1]; 
  Q = eye(2);
  s = sqrt(1+2*d);
  X = [s/d 1; 1 s];
  parout = [2, 1, 2];
  if nargout > 5,
    B = [0; 1];   R = 1;   C = eye(2);   Q0 = Q; 
  end

elseif index == 10,
  if nargin < 2,   d = 10^(-7);   else,   d = parin(1);   end
  A = [d+1 1; 1 d+1];
  G = eye(2);
  Q = d^2*eye(2);
  d1 = d + 1;
  x1 = (2*d1 + sqrt(2*d1^2 + 2) + sqrt(2)*d)/2;
  x2 = x1/(x1 - d1);
  X = [x1 x2; x2 x1];
  parout = [2, 2, 2];
  if nargout > 5,
    B = G;  R = G;  C = eye(2);  Q0 = Q;
  end

elseif index == 11,
  if nargin < 2,   d = 0.0;   else,   d = parin(1);   end
  A = [3-d 1; 4 2-d];
  G = ones(2);
  Q = [4*d-11 2*d-5; 2*d-5 2*d-2];
  X = [2 1; 1 1];
  parout = [2, 1, 2];
  if nargout > 5,
    B = ones(2,1);  R = 1;  C = eye(2);  Q0 = Q;
  end

elseif index == 12,
  if nargin < 2,   d = 10^6;   else,   d = parin(1);   end
  C = eye(3) - (2/3)*ones(3);
  A = d*(C*diag([1:3])*C);
  G = eye(3)/d;
  Q = C*diag([1/d 1 d])*C;        Q = (Q + Q')/2;
  dd = d^2;
  X = [dd+sqrt(1+dd^2), 2*dd+sqrt(d+4*dd^2), 3*dd+d*sqrt(1+9*dd)];   
  X = C*diag(X)*C;
  parout = [3, 3, 3];
  if nargout > 5,
    B = eye(3);    R = d*eye(3);    Q0 = diag([1/d, 1, d]);
  end

elseif index == 13,
  if nargin < 2,   d = 10^(-6);   else,   d = parin(1);   end
  A = [0 0.4 0 0; 0 0 0.345 0; 0 -0.524/d -0.465/d 0.262/d; 0 0 0 -1/d];
  B = [0; 0; 0; 1/d];
  G = B*B';
  Q = diag([1, 0, 1, 0]);
  parout = [4, 1, 2];
  if nargout > 5,
    R = 1;  C = [1 0 0 0; 0 0 1 0];  Q0 = eye(2);
  end

elseif index == 14,
  if nargin < 2,   d = 10^(-6);   else,   d = parin(1);   end
  A = [-d 1 0 0; -1 -d 0 0; 0 0 d 1; 0 0 -1 d];
  G = ones(4);
  Q = ones(4);
  parout = [4, 1, 1];
  if nargout > 5,
    B = ones(4,1);  R = 1;  C = ones(1,4);  Q0 = 1;
  end

elseif index == 15,
  if nargin < 2,   l = 20;   else,   l = parin(1);   end 
  n = 2*l - 1;
  A = zeros(n);
  B = zeros(n,l);
  C = zeros(l-1,n);
  for i=1:l-1
    j = 2*i;
    k = j-1;  
    A(k,k)   = -1;
    A(j,k)   = 1;
    A(j,j+1) = -1;
    B(k,i)   = 1;
    C(i,j)   = 1;
  end
  A(n,n) = -1;
  B(n,l) = 1;
  G  = B*B';
  Q  = 10*C'*C;
  parout = [n, l, l-1];
  if nargout > 5,
    R = eye(l);  Q0 = 10*eye(l-1);
  end

elseif index == 16,
  if nargin < 2,   n = 64;   else,   n = parin(1);    end
  A = diag(ones(n-1,1),-1) + diag(-2*ones(n,1)) + diag(ones(n-1,1),1);
  A(n,1) = 1;
  A(1,n) = 1;
  G = eye(n);
  Q = eye(n);
  cosvec = cos(2*pi*[0:n-1]/n);
  su = -2*ones(1,n) + 2*cosvec + sqrt(5*ones(1,n) - 8*cosvec + 4*cosvec.^2);  
  X = real(toeplitz(fft(su)/n));
  parout = [n, n, n];
  if nargout > 5,
    B = G;  R = G;  C = Q;  Q0 = Q;
  end

elseif index == 17,
  if nargin < 2,  
    n  = 21;
    lp = 0;  
  else,  
    n  = parin(1);  
    lp = length(parin);    
  end
  if lp < 3,  R  = 1;  else,  R  = parin(3);  end
  if lp < 2,  Q0 = 1;  else,  Q0 = parin(2);  end  
  if n > 1,  A = diag(ones(n-1,1),1);  else,  A = 0;  end
  B = flipud(eye(n,1));
  G = B/R*B';
  C = eye(1,n);
  Q = C'*Q0*C;
  X = sqrt(R*Q0);
  parout = [n, 1, 1];

elseif index == 18,
  if nargin < 2,
    n  = 100;
    lp = 0;  
  else
    n  = parin(1);  
    lp = length(parin);
  end
  if lp < 8,  c_int = [0.2, 0.3];  else,  c_int = parin(7:8);        end
  if lp < 6,  b_int = [0.2, 0.3];  else,  b_int = parin(5:6);        end
  if lp < 4,  gamma = 1.0;          else,  gamma = parin(4);        end
  if lp < 3,  beta  = 1.0;          else,  beta  = parin(3);        end
  if lp < 2,  alpha = 0.01;            else,  alpha = parin(2);        end
  evec = ones(n-1,1);

% Gram matrix
  M = (diag(evec,-1) + diag(4*ones(n,1)) + diag(evec,1)) / (6*(n+1));

% generate system matrices
  A = M \ (alpha*(n+1)*(diag(evec,-1) - diag(2*ones(n,1)) + diag(evec,1)));
  b = zeros(n,1);
  C = zeros(1,n);

% evaluate integrals
  for i = 1:n,
    b1 = max((i-1)/(n+1),b_int(1));        b2 = min((i+1)/(n+1),b_int(2));
    c1 = max((i-1)/(n+1),c_int(1));        c2 = min((i+1)/(n+1),c_int(2));
    if b1 >= b2,   
      b(i) = 0;
    else
      b(i) = b2 - b1;
      temp = min(b2,i/(n+1));
      if b1 < temp,
        b(i) = b(i) + (n+1)*(temp^2 - b1^2)/2 + i*(b1 - temp);
      end
      temp = max(b1,i/(n+1));
      if temp < b2,
        b(i) = b(i) - (n+1)*(b2^2 - temp^2)/2 - i*(temp - b2);
      end
      b(i) = beta * b(i);
    end 
    if c1 >= c2, 
      C(i) = 0;
    else
      C(i) = c2 - c1;
      temp = min(c2,i/(n+1));
      if c1 < temp,
        C(i) = C(i) + (n+1)*(temp^2 - c1^2)/2 + i*(c1 - temp);
      end
      temp = max(c1,i/(n+1));
      if temp < c2,
        C(i) = C(i) - (n+1)*(c2^2 - temp^2)/2 - i*(temp - c2);
      end
      C(i) = gamma * C(i);
    end 
  end
  B = M\b;
  G = B*B';
  Q = C'*C;
  parout = [n, 1, 1];
  if nargout > 5,
    R = 1.0;    Q0 = 1.0;
  end

elseif index == 19,
  if nargin < 2,   
    l  = 30;
    lp = 0;
  else
    l  = parin(1);
    lp = length(parin);
  end; 
  if lp < 4,        delta = 4.0;        else,        delta  = parin(4);        end
  if lp < 3,        kappa = 1.0;        else,        kappa  = parin(3);        end  
  if lp < 2,        mu    = 4.0;        else,        mu     = parin(2);        end
  n = 2*l;
  m = 2;
  p = 2*l;  
  temp = [1; 2*ones(l-2,1); 1];  
  A = [zeros(l), eye(l); ...
       -kappa*(diag(temp) - diag(ones(l-1,1),1) - diag(ones(l-1,1),-1))/mu,...
       -delta/mu*eye(l) ];
  B = [zeros(l,m); eye(l,1)/mu, flipud(-eye(l,1))/mu];
  G = B*B';
  C = [eye(l), zeros(l); zeros(l), eye(l)];
  Q = C;
  parout = [n, m, p];
  if nargout > 5,
    R = eye(m);   Q0 = eye(n);
  end

elseif index == 20,
  if nargin < 2,
    load carex20.mat
  else
    eval(['load ' parin]);
    if (size(mu,2) > 1),     mu = mu';        end
    if (size(delta,2) > 1),  delta = delta';  end
    if (size(gamma,2) > 1),  gamma = gamma';  end
    if (size(kappa,2) > 1),  kappa = kappa';  end
  end
  n = 2*l-1;
  A          = zeros(n,n);
  temp       = gamma./mu(1:l-1);
  temp2      = gamma./mu(2:l);
  A(1,2)     = -temp(1);
  A(2:l,2:l) = -(diag(temp+temp2) - diag(temp(2:l-1),1) - ...
                 diag(temp2(1:l-2),-1));
  temp       = delta./mu; 
  A(1:l,1:l) = A(1:l,1:l) - diag(temp);
  A(2:l,1)   = temp(2:l) - temp(1:l-1);
  for i = 2:l-1,
    A(i+1:l,i) = A(i+1:l,i) + temp(i:l-1) - temp(i+1:l);
  end
  temp     = kappa./mu(1:l-1);
  temp2    = kappa./mu(2:l);
  A(1,l+1) = -temp(1);
  A(2:l,l+1:n) = -(diag(temp + temp2) - diag(temp(2:l-1),1) - ...
                   diag(temp2(1:l-2),-1));
  A(l+1:n,:) = [zeros(l-1,1), eye(l-1), zeros(l-1)];
  temp = ones(l,1)./mu;
  B = [diag([temp(1);-temp(2:l)]) + diag(temp(1:l-1),-1); zeros(l-1,l)];
  G = B*B';
  C = [diag([1; gamma]), [zeros(1,l-1); diag(kappa)]]; 
  for i = 1:l,
    Q0(i) = sqrt(C(i,i)*C(i,i) + C(i,l+i-1)*C(i,l+i-1));
  end
  Q0 = diag(1./Q0);
  Q = Q0*C;  Q0 = Q0'*Q0;          
  Q = Q'*Q;       
  parout = [n, l, l];
  if nargout > 5,
    R = eye(l);  
  end        

else
  error('This example is not available!');
end  


