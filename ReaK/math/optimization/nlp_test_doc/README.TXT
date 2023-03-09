

           TEST EXAMPLES FOR CONSTRAINED NONLINEAR PROGRAMMING     
           ---------------------------------------------------



      Purpose:
      -------

      306 test problems of two previous collections
            
         W. Hock, K. Schittkowski, 'Test Examples for Nonlinear                    
         Programming Codes', Lecture Notes in Economics and Mathematical         
         Systems, Springer, No, 187, 1981                                           

         K. Schittkowski, 'More Test Examples for Nonlinear Programming 
         Codes', Lecture Notes in Economics and Mathematical Systems, 
         Springer, No, 282, 1987                                           
                                  
      are provided in form of Fortran source code together with a test frame.
      A decision is made which of the runs is successful, and performance 
      results are evaluated. With the default tolerances given, all problems 
      can be solved successfully by the code NLPQLP, a new version of the
      SQP implementation NLPQL of the author. Results of NLPQLP are 
      included for comparative studies.




      FILES:
      -----

      1. PROB.FOR: Fortran codes of the test problems of the two collections 
         mentioned above, which contain full documentation of the usage of
         the subroutines. Note that in some cases, gradients are not provided.

         IMPORTANT: The author gives no warranty that the the gradients, as 
         far as included, are correct!


      2. CONV.FOR: Interface between the individual test problem codes and 
         an available optimization routine, to facilitate the calling procedure 
         and to eb able to execute all test problems within a loop. 


      3. TESTP.FOR: Test program that executes 306 problems in a loop. The calling
         sequence for the the SQP code NLPQLP is included to give an example.
         Three different approximation formulae for gradient evaluations are 
         included, and a loop over randomly generated errors added to the problem
         functions is available. The code generates three output files listed below.


      4. TEST.DAT: Output file of the test frame containing numerical results 
         obtained by NLPQLP. Typical contants of TEST.DAT without lines generated 
         by the NLP routine:

         1  2  0  0  0 26 19 178  0.00000000E+00  0.73114619E-10  0.73E-10  0.00E+00 
         2  2  0  0  0 20 15 140  0.50426188E-01  0.50426193E-01  0.11E-06  0.00E+00 
         3  2  0  0  0 10 10  90  0.00000000E+00  0.16103740E-19  0.16E-19  0.00E+00 
         4  2  0  0  0  2  2  18  0.26666667E+01  0.26666667E+01  0.00E+00  0.00E+00 
         5  2  0  0  0  8  6  56 -0.19132230E+01 -0.19132230E+01  0.11E-10  0.00E+00 
          ...............


         The following data are listed, see TESTP.FOR for details:

         NTP      (Test problem number)  
         N        (Number of variables)
         ME       (Number of equality constraints)                                        
         M        (Number of constraints)                                        
         IFAIL    (Convergence criterion)                             
         NF       (Number of objective function evaluations)                     
         NDF      (Number of gradient evaluations of objective function)                                                     
         NEF      (Number of equivalent function evaluations, i.e. NF plus
                   number of function calls needed for gradient approximation)
         FEX      (Exact objective function value)
         F        (Cmputed objective function value)
         DFX      (Relative error in objective function)                   
         DGX      (Sum of constraint violations including bound violations)                                              



      5. TEST.TEX: Same as above, buyt with Latex separators.



      6. RESULT.DAT: The following summary is shown:

         Flag for evaluating gradients  
         Tolerance for gradient approximation
         Termination accuracy for NLP routine
         Randomly generated error added to objective
         Total number of test runs
         Number of successful test runs
         - constraint violation less than squared tolerance and
           either error in objective less than tolerance or
           termination criteria of NLP routine satisfied
         Number of better test runs
         - constraint violation less than squared tolerance and
           error in objective less than -tolerance
         Number of local solutions obtained
         Number of runs without satisfying termination accuracy
         - as indicated by NLP routine, i.e. IFAIL>0
         Tolerance for determining successful return
         Average number of function evaluations
         - NLFUNC calls for successful returns
         Average number of gradient evaluations
         - NLGRAD calls for successful returns
         Average number of equivalent function calls      
         - additionlly counted function calls for gradient 
           approximations
         Total execution time over all test runs (sec)

 

      7. TEMP.DAT: Contains the same data in one row

     


      Compilation and Link:
      --------------------

      To run the code, the following files must be linked:

      TESTP.OBJ     : test program
      PROB.OBJ      : individual test problems
      CONV.OBJ      : interface for executing individual test problems 
      NLPQLP.OBJ    : NLP solver  
      QL.OBJ        : auxiliary routine required by solver NLPQLP 
 



      References:
      ---------
  
      K. Schittkowski. NLPQL: A Fortran subroutine solving
      constrained nonlinear programming problems, Annals of
      Operations Research, Vol.5 (1985/6), 485-500.

      K. Schittkowski. NLPQLP: A Fortran implementation 
      of a sequential quadratic programming algorithm with 
      distributed and non-monotone line search, Report, 
      Department of Computer Science, University of Bayreuth, 
      2006





      Author (C):
      ----------

      Prof. K. Schittkowski       
      Dept. of Computer Science               
      University of Bayreuth
      D-95440 Bayreuth
      Germany

      http:\\www.klaus.schittkowski.de
      klaus.schittkowski@uni-bayreuth.de
      


      Bayreuth, April 2008
