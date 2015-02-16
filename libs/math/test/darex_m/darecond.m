function kappa = darecond(X,A,B,Q,R,S)
%DARECOND
%
% Computes the condition number of discrete-time algebraic Riccati
% equations (DARE) as described in [1].
% The DARE is given by
%                                                  -1
% (I)  0 = DR(X) = A'XA - X - (A'XB + S) (R + B'XB)  (B'XA + S') + Q
%
% Since the computation of the condition number involves Kronecker
% products, sparse matrix functions are used as far as possible.
% For a dimension of n > 10 (i.e., the resulting linear systems are
% of order n*n > 100), 2-norms are estimated via power iteration.
%
% Input: 
%   - X            : (approximate) solution of the DARE.
%   - A, B, Q, R, S: coefficient matrices of the DARE as in (I). (S is
%                    optional, on default, S = 0.)
%
%  Output:     
%   - kappa        : The DARE condition number DARE as described in [1]. 
%
%
%  Reference:
%
%  [1] T. GUDMUNDSSON, C. KENNEY, A.J. LAUB: 'Scaling of the Discrete-Time
%      Algebraic Riccati Equation to Enhance Stability of the Schur 
%      Solution Method', IEEE Transactions on Automatic Control, vol. 37, 
%      no. 4, pp. 513-518, 1992.
 
%  Peter Benner, Volker Mehrmann (TU Chemnitz-Zwickau, Germany),
%  Alan J. Laub (University of California at Santa Barbara)
%  12-15-1995 

error(nargchk(5,6,nargin))
[na,ma] = size(A);
if (na~=ma),  error('Input matrix A must be square.'),  end
n = na;
[nb,mb] = size(B);
if (nb~=n),
  error('Input matrices A and B must have same number of rows.'),  
end  
m = mb;
[nq,mq] = size(Q);
if (nq~=mq),  error('Input matrix Q must be square.'),  end
if (nq~=n),   error('Input matrices A and Q must have same dimension.'),  end
[nr,mr] = size(R);
if (nr~=mr),  error('Input matrix R must be square.'),  end
if (mr~=m),
  error('Input matrices B and R must have same number of columns.'),  
end
if (nargin > 5),
  if size(S) ~= [n,m],  
    error('Input matrices S and B must have the same size.')
  end
end
if rank(R) < m,
  error('Input matrix R must be nonsingular.')
end

G    = sparse(B/R*B');
Anrm = norm(A,'fro');
Gnrm = norm(G,'fro');
Qnrm = norm(Q,'fro');
G    = speye(n) + sparse(G*X);
AS   = sparse(G \ A);

P  = sparse(speye(n*n) - kron(AS',AS'));
P  = sparse(speye(n*n)/P);
for i = 1:n, 
  for j = 0:n-1, 
    perm = [perm,i+j*n]; 
  end, 
end
S = speye(n*n);
S = S(perm,:);

L  = sparse(AS'*X);
DA =  Anrm*P*(kron(speye(n),L) + kron(L,speye(n))*S);
L  =  sparse(A'*X / G);
DG = -Gnrm*P*sparse(kron(L,L));;
DQ = -Qnrm*P;  

if (n < 11),
  Lnrm = norm(full([DA,DG,DQ]));
else

% 2-norm estimation by power iteration 
% (equivalent to  Lnrm  = nrm2est([DA,DG,DQ],1e-4); )
  x1    = sum(abs(DA))';
  x2    = sum(abs(DG))';
  x3    = sum(abs(DQ))';
  Lnrm  = norm(x1)^2;
  Lnrm  = Lnrm + norm(x2)^2;
  Lnrm  = sqrt(Lnrm + norm(x3)^2);
  x1    = x1/Lnrm;
  x2    = x2/Lnrm;
  x3    = x3/Lnrm;
  dummy = 0;
  while abs(Lnrm - dummy) > 1e-4*Lnrm,
    dummy = Lnrm;
    y     = DA*x1;
    y     = y + DG*x2;
    y     = y + DQ*x3;
    Lnrm  = norm(y);
    x1    = DA'*y;
    x2    = DG'*y;
    x3    = DQ'*y;
    xnrm  = norm(x1)^2;
    xnrm  = xnrm + norm(x2)^2;
    xnrm  = sqrt(xnrm + norm(x3)^2);
    x1    = x1/xnrm;
    x2    = x2/xnrm;
    x3    = x3/xnrm;
  end
end

kappa = Lnrm/norm(X,'fro');
