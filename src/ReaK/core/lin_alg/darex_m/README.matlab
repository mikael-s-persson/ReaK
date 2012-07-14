          DAREX - A Benchmark Collection for DARE
         =========================================

The MATLAB function DAREX is designed to generate all examples of
discrete-time algebraic Riccati equations (DARE)

                  T            T               T    -1  T       T
(I)  0 = DR(X) = A X A - X - (A X B + S) (R + B X B)  (B X A + S ) + Q

collected in [1]. Here, A, Q, and X are n-by-n matrices, B and S are
n-by-m, and R is m-by-m. The coefficient matrices Q and R are symmetric
and usually, the required solution matrix X is symmetric, too. The
coefficient matrix Q is often given in factored form as 

                T 
(II)     Q  =  C  Q0 C,

where C is p-by-n and Q0 is p-by-p. This factorization often arises
in (but is not limited to) DARE coming from control theory. Also, if R
is nonsingular, the DARE (I) is equivalent to 

                          T             T             -1 
(III)   0  =  DR(X)  =  AA X AA - X - AA X G (I + X G)  X AA + QQ 

where
                     -1 T                     -1 T
        AA =  A - B R  S ,       QQ =  Q - S R  S , 

and
                 -1  T
(IV)    G  =  B R   B .

The required solution X often has some particular or extremal
properties, for instance in control theory one is usually concerned
with the "stabilizing" solution of (I), i.e., a solution such that the
"feedback gain matrix" 
                                 -1  T       T
(V)	F  =  A  -  B (R + B X B)  (B X A + S )
 
has all its eigenvalues inside the unit circle.

The benchmark collection [1] consists of examples which can be used for
testing purposes in the construction of new numerical methods to solve
DAREs or as a reference set for the comparison of methods. Although
the presented benchmark examples have a control-theoretic background,
they can be considered as examples for general DARE of the form (I).
  
---------------------------------------------------------------------------

IMPLEMENTATION:
===============

All benchmark examples from [1] can be created by the matlab function
DAREX. The calling sequence is

	[A,B,Q,R,S,X,parout,G,C,Q0]=darex(index,parin)

where any number of output parameters less than or equal to 10 is
allowed. The input parameter 'index' is the number of the required
example corresponding to [1] and 'parin' may contain optional
parameters used to generate the coefficient matrices. 
On return, A,B,Q,R, and S contain the coefficient matrices of the DARE
as given in (I) whereas G is the matrix G from (IV) (only available
if R is nonsingular) and C,Q0 contain the factors of the matrix Q
as defined by (II). 

A detailed description of all input and output parameters of DAREX
is given in the prolog of the MATLAB function DAREX which is shown,
for instance, by the command

	>> help darex

in a MATLAB environment. Some examples of how to use DAREX are given
in the Section EXAMPLES below.

If an analytical solution of the DARE is available, it will be returned,
too (in the output argument X). This returned solution is generally
the "stabilizing" solution in the control-theoretic sense, i.e., the
feedback gain matrix F as in (V) will have all its eigenvalues inside the 
unit circle. If such an analytical solution is available, the DARE
condition number given in [2] is computed by the MATLAB function DARECOND
provided in a separate function file together with darex.m.

Some of the examples have one integer and/or one or several real input
parameters. These can be supplied by the user or they are set by default.

---------------------------------------------------------------------------

CONTENTS:
=========

You can receive the file darex_m.tar.Z by anonymous ftp at

	ftp.tu-chemnitz.de

from the directory

	pub/Local/mathematik/Benner

(observe the capital "L" in Local !) where you can also receive a
postscript version of the preprint [1] (this file is called blm2.ps.Z).

Then by decompressing darex_m.tar.Z via

	uncompress darex_m.tar.Z

or

	compress -d darex_m.tar.Z

and extracting the resulting file darex.tar by

	tar xf darex_m.tar

a directory darex_m is created containing the following files : 

darecond.m    - The MATLAB function DARECOND for computing the condition 
                number of discrete-time algebraic Riccati equations given
                in [2]. 
darex.m       - The MATLAB function DAREX for generating all the
                benchmark examples presented in [1].
README.matlab - This file.

---------------------------------------------------------------------------

EXAMPLES:
=========

In the following we will describe how to use DAREX by showing some
sample calls to DAREX and the resulting output. 
We used MATLAB Version 4.2a on an HP Apollo 705 work station with
operating system HP-UX 9.01.

A) We generate Example 1 of [1] where neither G nor C or Q0 are needed
and show the resulting ouptut.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>> [A,B,Q,R,S,X,parout]=darex(1)
Condition is infinite

A =

   4.0000e+00   3.0000e+00
  -4.5000e+00  -3.5000e+00


B =

     1
    -1


Q =

     9     6
     6     4


R =

     1


S =

     0
     0


X =

   1.4562e+01   9.7082e+00
   9.7082e+00   6.4721e+00


parout =

  Columns 1 through 6 

   2.0000e+00   1.0000e+00   2.0000e+00   1.1499e+02   5.0000e-01   2.1034e+01

  Columns 7 through 8 

          Inf   1.8845e+01

>> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The message 'Condition is infinite' is caused by the computation of the
condition number of A since A is singular.
The first three components of the vector 'parout' contain in consecutive
order n, m, and p. The next two components are the 2-norm condition number 
of the coefficient matrix A and the absolute value of the eigenvalue of F
closest to the unit circle and parout(6:7) contain the 2-norm and 2-norm
condition number of the solution matrix X. The last component of 'parout'
is the DARE condition number as given in (II). Here, X is singular to
working precision and thus, parout(7) is set to infinity.

___________________________________________________________________________

B) Now we want to generate Example 4 of [1]. This example has a nonzero S
matrix but no analytical solution is known and the matrices G,C,Q0 from
(II) and (IV) are to be returned.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>> [A,B,Q,R,S,X,parout,G,C,Q0]=darex(4)
A =

     0     1
     0    -1


B =

     1     0
     2     1


Q =

   -0.3636   -0.3636
   -0.3636    0.6364


R =

     9     3
     3     1


S =

     3     1
    -1     7


X =

     []


parout =

     2     2     2   Inf


G =

     [] 


C =

     1     0
     0     1


Q0 =

   -0.3636   -0.3636
   -0.3636    0.6364

>> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Since X is not known, parout contains only 4 parameters because the other
values would require the knowledge of X.   
Note that here, R is singular. Thus, G cannot be computed and the empty
matrix is returned in G. 

___________________________________________________________________________

C) Next, we create the coefficient matrices corresponding to Example 14
from [1]. This example has four real parameters. The first two parameters
are set to 1.57 and 12.34 whereas the other two parameters are set by
default. None of the matrices G,C,Q0 is required.

>> [A,B,Q,R,S,X,parout]=darex(14,[1.57,12.34])
Condition is infinite

A =

  -6.8599e+00            0            0            0
   1.0000e+00            0            0            0
            0   1.0000e+00            0            0
            0            0   1.0000e+00            0


B =

   7.8599e+00
            0
            0
            0


Q =

     0     0     0     0
     0     0     0     0
     0     0     0     0
     0     0     0     1


R =

   2.5000e-01


S =

     0
     0
     0
     0


X =

   1.1898e+00            0            0            0
            0   1.0000e+00            0            0
            0            0   1.0000e+00            0
            0            0            0   1.0000e+00


parout =

  Columns 1 through 6 

   4.0000e+00   1.0000e+00   1.0000e+00          Inf   2.3253e-02   1.1898e+00

  Columns 7 through 8 

   1.1898e+00   2.6484e+02

>> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

___________________________________________________________________________

D) The last example given here has the number 15 in [1]. This example has 
one integer parameter which defines the order of the DARE (1) and one
real parameter defining the matrix R which, in this example, is 1-by-1.
Setting n to 3 and R to 4 can be achieved by

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

>> [A,B,Q,R,S,X,parout,G]=darex(15,[3,4])
Condition is infinite

A =

     0     1     0
     0     0     1
     0     0     0


B =

     0
     0
     1


Q =

     1     0     0
     0     1     0
     0     0     1


R =

     4


S =

     0
     0
     0


X =

     1     0     0
     0     2     0
     0     0     3


parout =

  Columns 1 through 6 

   3.0000e+00   1.0000e+00   3.0000e+00          Inf            0   3.0000e+00

  Columns 7 through 8 

   3.0000e+00   1.9854e+00


G =

     0     0     0
     0     0     0
     0     0     1

>> 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

---------------------------------------------------------------------------

HELP AND BUGS:
==============

If you have trouble using the functions or if you find any bugs please send
an e-mail message reporting the problem or bug to  

        benner@mathematik.tu-chemnitz.de

We will get in touch with you as soon as possible. We will also appreciate
any proposal for examples to be included in further releases of DAREX.

---------------------------------------------------------------------------

REFERENCES:
===========

[1] P. Benner, A.J. Laub, and V.Mehrmann
    'A Collection of Benchmark Examples for the Numerical Solution of
     Algebraic Riccati Equations II: Discrete-Time Case',
    Technical Report SPC 95_23, TU Chemnitz-Zwickau (FRG), 1995.
    (A postscript version is available by anonymous FTP on
     ftp.tu-chemnitz.de in directory /pub/Local/mathematik/Benner,
     get blm2.ps.Z .)

[2] T. Gudmundsson, C. Kenney, A.J. Laub 
    'Scaling of the Discrete-Time Algebraic Riccati Equation to Enhance
     Stability of the Schur Solution Method', 
    IEEE Transactions on Automatic Control, vol. 37, no. 4, pp. 513-518, 
    1992. 

---------------------------------------------------------------------------
 
CONTRIBUTORS:
=============

Peter Benner and Volker Mehrmann
Fakultaet fuer Mathematik
Technische Universitaet Chemnitz-Zwickau 
09107 Chemnitz (FRG)
e-mail: benner@mathematik.tu-chemnitz.de
        mehrmann@mathematik.tu-chemnitz.de

Alan J. Laub
Department of Electrical and Computer Engineering
University of California
Santa Barbara, CA 93106-9560 (USA)
e-mail: laub@ece.ucsb.edu

---------------------------------------------------------------------------

Peter Benner, January 31, 1996.
