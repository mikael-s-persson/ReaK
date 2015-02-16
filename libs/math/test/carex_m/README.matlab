          CAREX - A Benchmark Collection for CARE
         =========================================

The MATLAB function CAREX is designed to generate all examples of
continuous-time algebraic Riccati equations (CARE)

                     T 
(I)     0  =  Q  +  A  X  +  X A  -  X G X

collected in [1]. Here, A, G, Q, and X are real n-by-n matrices. G and Q
are symmetric and usually, the required solution matrix X is symmetric,
too. The coefficient matrices G and Q are often given in factored form as

                 -1  T            T 
(II)    G  =  B R   B ,    Q  =  C  Q0 C,

where B is n-by-m, R is m-by-m, C is p-by-n, and Q0 is p-by-p. This
factorization often arises in (but is not limited to) CARE coming from
control theory. The required solution X often has some particular or
extremal properties, for instance in control theory one is usually
concerned with the "stabilizing" solution of (I), i.e., a solution
such that the "feedback gain matrix"

	F  =  A  -  G X
 
has all its eigenvalues in the open left half plane.
The benchmark collection [1] consists of examples which can be used for
testing purposes in the construction of new numerical methods to solve
CAREs or as a reference set for the comparison of methods. Although
the presented benchmark examples have a control-theoretic background,
they can be considered as examples for general CARE of the form (I).
  
---------------------------------------------------------------------------

IMPLEMENTATION:
===============

All benchmark examples from [1] can be created by the matlab function
CAREX. The calling sequence is

	[A,G,Q,X,parout,B,R,C,Q0]=carex(index,parin)

where any number of output parameters less than 10 is allowed. 'index'
is the number of the required example corresponding to [1].
On return, G and Q contain the coefficient matrices G and Q from (I)
whereas B,R, and C, Q0, respectively, contain the factors from (II).

A detailed description of all input and output parameters of CAREX
is given in the prolog of the MATLAB function CAREX which is shown,
for instance, by the command

	>> help carex

in a MATLAB environment. Some examples of how to use CAREX are given
in the Section EXAMPLES below.

If an analytical solution of the CARE is available, it will be returned,
too (in the output argument X). This returned solution is generally
the "stabilizing" solution in the control-theoretic sense, i.e., the
feedback gain matrix

	F  =  A  -  G X

has all its eigenvalues in the open left half plane.
If such an analytical solution is available, the upper and lower bounds of
the CARE condition number given in [2] are computed which requires the
solution of three Lyapunov equations. These equations are solved using the
MATLAB function LYAP from the MATLAB CONTROL TOOLBOX. If this is not
available in your computer environment, please contact the authors
(see the HELP AND BUGS section in this file).     

Most of the examples have one integer and/or one or several real input
parameters. These can be supplied by the user or they are set by default.
Two of the examples require data files which are supplied together with
carex.m. 

---------------------------------------------------------------------------

CONTENTS:
=========

You can receive the file carex_m.tar.Z by anonymous ftp at

	ftp.tu-chemnitz.de

from the directory

	pub/Local/mathematik/Benner

(observe the capital "L" in Local !) where you can also receive a
postscript version of the preprint [1] (this file is called blm1.ps.Z).

Then by decompressing carex_m.tar.Z via

	uncompress carex_m.tar.Z

or

	compress -d carex_m.tar.Z

and extracting the resulting file carex.tar by

	tar xf carex_m.tar

a directory carex_m is created containing the following files : 

carecond.m    - The MATLAB function CARECOND for computing lower and
                upper bounds for the condition number of algebraic
                Riccati equations given in [2].
carex.m       - The MATLAB function CAREX for generating all the
                benchmark examples presented in [1].
carex6.mat    - Data file in internal binary MATLAB format required by
                CAREX for generating Example 6 of [1].    
carex20.mat   - Data file in internal binary MATLAB format required by
                CAREX for generating Example 20 of [1].    
README.matlab - This file.

---------------------------------------------------------------------------

EXAMPLES:
=========

In the following we will describe how to use CAREX by showing some
sample calls to CAREX and, in part, the resulting output. 
We used MATLAB Version 4.2a on an HP Apollo 705 work station with
operating system HP-UX 9.01.

A) We generate Example 1 of [1] where Q0 and C are not needed and show
the resulting ouptut.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
>> [A,G,Q,X,parout,B,R]=carex(1)

A =

     0     1
     0     0


G =

     0     0
     0     1


Q =

     1     0
     0     2


X =

     2     1
     1     2


parout =

  Columns 1 through 7 

    2.0000    1.0000    2.0000    2.4142    5.8284    3.0000    3.0000

  Columns 8 through 9 

    5.0399    4.7333


B =

     0
     1


R =

     1  
	
>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The first three components of the vector 'parout' contain in consecutive
order n, m, and p. The next two components are the 2-norm and 2-norm
condition number of the Hamiltonian matrix 

	    (  A  -G  )
	H = (         )  
            ( -Q  -A' )

corresponding to the CARE (I). parout(6:7) contains the 2-norm and 2-norm
condition number of the solution matrix X (which are both equal to 3.0000
in this example). The last two components of 'parout' are the lower and
upper bound of the CARE condition number as given in (II). 


B) Now we want to generate Example 8 of [1]. This example has as
parameter one real number. We want to set this number to 3.14159. 
This is achieved by
 
>> [A,G,Q,X,parout,B,R,C,Q0]=carex(8,3.14159);

Here we suppress printing the output. From this calling sequence, we also
obtain the factors C and Q0 as in (II).


C) Next, we create the coefficient matrices corresponding to Example 15
from [1]. This example has one integer parameter which we choose equal to
6. Only G and Q are to be returned in product form as in (I). 

>> [A,G,Q,X,parout]=carex(15,6);

Since for this example, no analytical solution is known, X contains an
empty matrix :

>> X

X =

     []

The output vector 'parout' has therefore only five components,

>> parout

parout =

   11.0000    6.0000    5.0000   10.0000   53.0309

The order of the problem is given by n = 2*N - 1 = 11 (where N is the input 
parameter) and m = N = 6, p = N - 1 = 5. Since X is not known there is no
information about its 2-norm and condition number available. Also, the
bounds for the CARE condition numbers require the knowledge of X and thus
are not returned.


D) The last example given here has the number 19 in [1]. This example has 
four parameters, the order L of an underlying second-order sytem (yielding
n = 2*L), and three real parameters mu, delta, and kappa. Here, L, mu, and
delta are defined by L = 3, mu = 10.0, delta = 3.0, and the remaining
parameter kappa is set by default. One possible calling sequence and the
resulting output are

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

>> [A,G,Q]=carex(19,[3,10.0,3.0]) 

A =

         0         0         0    1.0000         0         0
         0         0         0         0    1.0000         0
         0         0         0         0         0    1.0000
   -0.3000    0.3000         0   -0.4000         0         0
    0.3000   -0.6000    0.3000         0   -0.4000         0
         0    0.3000   -0.3000         0         0   -0.4000


G =

         0         0         0         0         0         0
         0         0         0         0         0         0
         0         0         0         0         0         0
         0         0         0    0.0100         0         0
         0         0         0         0         0         0
         0         0         0         0         0    0.0100


Q =

     1     0     0     0     0     0
     0     1     0     0     0     0
     0     0     1     0     0     0
     0     0     0     1     0     0
     0     0     0     0     1     0
     0     0     0     0     0     1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

---------------------------------------------------------------------------

HELP AND BUGS:
==============

If you have trouble using the functions or if you find any bugs please send
an e-mail message reporting the problem or bug to  

        benner@mathematik.tu-chemnitz.de

We will get in touch with you as soon as possible.

---------------------------------------------------------------------------

REFERENCES:
===========

[1] P. Benner, A.J. Laub, and V.Mehrmann
    'A Collection of Benchmark Examples for the Numerical Solution of
     Algebraic Riccati Equations I: Continuous-Time Case',
    Technical Report SPC 95_22, TU Chemnitz-Zwickau (FRG), 1995.
    (A postscript version is available by anonymous FTP on
     ftp.tu-chemnitz.de in directory /pub/Local/mathematik/Benner)

[2] C. Kenney and G. Hewer 
    'The sensitivity of the algebraic and differential Riccati equations', 
    SIAM J. Control Optim., vol. 28 (1990), pp.50-69.

---------------------------------------------------------------------------
 
CONTRIBUTORS:
=============

Peter Benner and Volker Mehrmann
Fakultaet fuer Mathematik
Technische Universitaet Chemnitz-Zwickau (FRG)
e-mail: benner@mathematik.tu-chemnitz.de
        mehrmann@mathematik.tu-chemnitz.de

Alan J. Laub
Department of Electrical and Computer Engineering
University of California
Santa Barbara, CA 93106-9560 (USA)
e-mail: laub@ece.ucsb.edu

---------------------------------------------------------------------------

Peter Benner, November 13, 1995.
