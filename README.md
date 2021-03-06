CURRENT VERSION: Beta 1.6
_______________________________________________________________________________________
Installation Instructions:

Pre-Install) 
Required libraries - only c++ std 11 libraries are required. All other functions were
written in house and are included with the program.

1) Export Install Directory with export SCF_HOME=%INSTALL_DIR%/HFSCFCISProgram/
2) Open configure_install file and set compiler (if other than g++) and 
   necessary flags. 
3) Run the command ./INSTALL configure_install in the src folder

_______________________________________________________________________________________
Example input file:


-----------------STARTS BELOW THIS LINE---------------------
$PARAM
3               !Number of contracted gaussians per orbital "n of STO-nG"
0               !Charge on the molecule
1.0E-30         !SCF Jacobi Convergence
1.0E-7          !SCF Energy Convergence
1               !Output verbosity (0 minimal, 1 all)
$ENDPARAM

Notes can be written outside of the delimiters without causing an error.

$COORD
        1       -1.12798434        1.545547           2.2132132
        6       -0.1224696362       -3.1174832504      -1.5015549412
        8        2.1012819183       4.0397895568      -1.1329875333
        6        3.2785162937       3.3434825848      -0.7159896664
        6       14.1137391171      14.5455261160      -3.6515096380
$ENDCOORD

$END OF FILE$
-------------------ENDS ABOVE THIS LINE-----------------------

_______________________________________________________________________________________
Extra Information:

Program Status: The program keeps ~6 to 7 sig figs. Some molecules are still having 
trouble converging to the correct energy. I.E. CN+ or CN- or similar molecules.
However, molecules as large as nitro benzene have been tested and give correct energy
to ~ 6 sig figs.

Admittedly the program is a memory hog. It explicitely stores and calculates both the 
upper and lower triangular of all symmetric matricies. It also stores all zero values 
from orthogonal atomic orbitals. This can be changed very easily, however, I left it 
this way for testing of the integrals. 

Included methods:
	1) HF

Included Basis Sets:
	1) STO-3G
	2) STO-6G

Following is a run through of how the program works.
_______________________________________________________________________________________
1) Integral Calculations
	1 - Overlap integral calculation
		*Working very well in comparison with other programs
		*Stored in S

	2 - Kinetic energy calculation
		*Working very well in comparison with other programs
		*Stored in Hcore

	3 - Nulcear attraction calculation 
		*Working very well in comparison with other programs
		*Stored in Hcore -- added with values from the kinetic energy and 
		other nuclear centers

	4 - Electron repulsion integrals
		*Working very well in comparison with other programs
		NOTE: The issue with the aux function is now resolved. Problem was 
		overflow of a pow(x,N)/fact(N) factor in the expansion. Fixed by 
		generalizing to: 
			x/n * x/(n-1) * ... * x/2 * x/1
		*Stored in FEI -- Indexing is still slow. Need better way to do this.
		This indexing is by far the slowest part of the program aside from the
		computation of the 2 electron integrals.

_______________________________________________________________________________________
2) Various Linear Algebra Routines
	1 - Matrix inversion and root (all in one)
	2 - Matrix Multiplication with Transpose options
	3 - Matrix Multiplication with constant
	4 - Determinant Calculation - NOT BEING USED
	5 - Null matrix producer
	6 - Cofactor matrix calculator - NOT BEING USED
	7 - Set_EPS functions - Sets values in working matricies to zero if below a
				specific value.
	8 - Matrix Transpose - NOT BEING USED
	9 - Null_Set - Sets a given matrix to a null matrix
	10 - Multi_Elements - Multiplies the elements of a matrix by the elements in
			      another matrix. NOT MATRIX MULTIPLICATION.

_______________________________________________________________________________________
3) Jacobi Algorithm for calculating the eigenvalues and eigenvectors

Works perfect for small matricies of 4x4 or less, The entire thing uses doubles for 
obtaining more precise values.

Basically, the program takes a matrix S and produces V (eigenvector matrix) and 
D (a set of eigenvalues) using the Jacobi iterative algorithm.

	S = V * D * V^(H)

Currently, this function works properly.
_______________________________________________________________________________________
4) SCF Main Implemented

Currently carries out the linear algebra calculations in the following order:

1) Takes S, calculates S^(-1/2) with the following steps:
	1 - Uses the jacobi algorithm and the relation,

		S^(1/2) = V * D^(1/2) * V^T

		Effectively, we calculate D and V with Jacobi, then use the fact that 
		D^(1/2) is simply the square root of the diagonal D matrix obtained from 
		jacobi.

	2 - Calculate the inverse of S^(1/2) by,

		S^(1/2) = [1 / Det(S)] * Adjoint(S)

2) Calculate F' = S^(-1/2) * F * S(-1/2)

3) Calculate E and C' from F'

	This is simply running Jacobi for eigenvalues E and eigenvectors C'

4) Calculate C = S^(-1/2) * C' 

	C is the coefficient matrix and E from (3) is an array of energy values

5) Next we produce the charge density matrix (Rho) from the C matrix using the lowest 
   energy occupied orbitals. These lowest energy occupied orbitals are the eigen 
   vectors that correspond to the N / 2 lowest energy eigen values.

6) Calculate electronic energy from Rho, Hcore, and the 2 electron integrals.

7) Check convergence, repeat if necessary, end if not. Convergence is determined by 
   the change in energy between steps. This can be set in the input file.

_______________________________________________________________________________________
8) Mulliken Population Analysis Program

Currently calculates individual orbital populations, total atom populations, total
atom Mulliken charge, and total molecule population/charge.
 
_______________________________________________________________________________________
9) Huckel Method Initial Density Production

Implemented the Huckel method to produce the initial density matrix "guess" to help 
program convergence. The program produces the valence density matrix, then extends 
that to a full density matrix by adding the core and empty orbitals to the valence 
desity matrix.

_______________________________________________________________________________________
10) CIS Added to the program
