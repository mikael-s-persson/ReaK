function [KU,KL] = carecond(A,G,Q,X)
%CARECOND
%
% Computes upper and lower bounds for the condition number of
% continuous-time algebraic Riccati equations (CARE)
%
% (1)  0 =  Q  +  A' X  +  X A  -  X G X.  
%
% The condition number for CARE defined in [1] is given by
%
%                     ||DX||    ||DA||       ||DG||       ||DQ||
%   K =  lim    sup{ ------- ; ------ <= d, ------ <= d, ------ <= d }
%       d -> 0       d*||X||    ||A||        ||G||        ||Q||
%
% where DA = A - AA, DG = G - GG, DQ = Q - QQ,  DX = X - XX and X, XX are
% the positive semidefinite solutions of the CARE (1) resp. the CARE with
% disturbed data  
%
%   0 =  QQ  +  AA' XX  +  XX AA  -  XX GG XX.
%
% All matrices are square n-by-n matrices and G,Q,X as well as GG,QQ,XX 
% are symmetric positive semidefinite. ||.|| denotes the 2-norm of a
% matrix.
% In [1], it is proved that  KL/3 <= K <= KU, where the upper (KU) and 
% lower (KL) bounds are given by  
%   
%          ||H_0||*||Q|| + 2*sqrt(||H_0||*||H_2||)*||A|| + ||H_2||*||G||  
%   KU  = --------------------------------------------------------------- ,
%                                     ||X||
%   
%          ||H_0||*||Q|| + 2*||H_1||*||A|| + ||H_2||*||G||  
%   KL  = ------------------------------------------------- .
%                               ||X||
%
% Here the H_k, k = 0,1,2, are the solutions of the Lyapunov equations
%
%   (A - GX)' H_k + H_k (A - GX) = -X^k,   k = 0,1,2.
%
% CALLING SEQUENCE:
%
% [KU,KL] = carecond(A,G,Q,X)
%
% where KU,KL,A,G,Q,X are defined above. 
%
% Reference:
% [1] C. KENNEY, G. HEWER: 'The sensivity of the algebraic and differential
%     Riccati equation', SIAM J. Control. Optim., vol. 28 (1990), pp.50-69.
 
% Peter Benner (TU Chemnitz-Zwickau, Germany),  06-06-1995 
%
%  For questions concerning this M-file, send e-mail to
%
%        benner@mathematik.tu-chemnitz.de

error(nargchk(1,4,nargin))
[na,ma] = size(A);
if na~=ma,        error('Input matrix A must be square.'),          end
n = na;
[ng,mg] = size(G);
if ng~=mg,        error('Input matrix G must be square.'),          end
if ng~=n,         error('Incorrectly sized matrix G.'),                end
[nq,mq] = size(Q);
if nq~=mq,        error('Input matrix Q must be square.'),          end
if nq~=n,         error('Incorrectly sized matrix Q.'),                end
[nx,mx] = size(X);
if nx~=mx,        error('Input matrix X must be square.'),          end
if nx~=n,         error('Incorrectly sized matrix X.'),                end

Anorm = norm(A);
Gnorm = norm(G);
Qnorm = norm(Q);
Xnorm = norm(X);

F = A - G*X;
 
H0 = lyap(F',eye(n));
H1 = lyap(F',X);
H2 = lyap(F',X*X);

H0norm = norm(H0);
H1norm = norm(H1);
H2norm = norm(H2);

KU = (H0norm*Qnorm + 2*sqrt(H0norm*H2norm)*Anorm + H2norm*Gnorm)/Xnorm; 
KL = (H0norm*Qnorm + 2*H1norm*Anorm + H2norm*Gnorm)/Xnorm; 
