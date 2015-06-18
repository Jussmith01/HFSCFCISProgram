#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iostream>
#include <iomanip>

#include "../utils_cpp/lib_includes.h"
#include "../classes/classes.h"
#include "scf_classes.h"
#include "../funcs_cpp/func_classes.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                   Find max val such that 1 <= i < j <= N
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/29/2014  ----
                   ----Modified By:               ----
*/
void find_max(double *A,int *max, int N)
{
	double testval = 0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			double valij = abs(A[i + j * N]);
			if (valij > testval)
			{
				testval = valij;
				max[0] = i;
				max[1] = j;
			}
		}
	}
}

//________________________________________________________________________//
//      *************************************************************     //
//                 Check Convergence of Jacobi Algorithm
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/29/2014  ----
                   ----Modified By:               ----
*/
double check_convergence(double *A,int N)
{
	double pos = 0;
	double tol = 1.0E-13;

        for (int i = 0; i < N; ++i)
        {
                for (int j = i + 1; j < N; ++j)
                {
			if (i != j)
			{
				double tmpval = (double)A[j + i * N] * (double)A[j + i * N];
				if (tmpval > tol)
					{pos += tmpval;}
                        }
                }
        }

	return pos;
}

//________________________________________________________________________//
//      *************************************************************     //
//                            Hypot Functions
//      *************************************************************     //
/*                 ----Author:                    ----
                   ----Date Modified: 10/29/2014  ----
                   ----Modified By:               ----
*/
double hypot(double x,double y)
{
    double t;
    x = abs(x);
    y = abs(y);
    t = min(x,y);
    x = max(x,y);
    t = t/x;
    return x * sqrt(1 + t * t);
}

//________________________________________________________________________//
//      *************************************************************     //
//                     Produce the Rotation Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/29/2014  ----
                   ----Modified By:               ----
*/
extern void Make_Rot_mat(double &C,double &S,double *A,int i, int j,int N)
{
	//Find submatrix eigen vector for use in rotation matrix
        double a = A[i + i * N];
        double b = A[j + i * N];
        double d = A[j + j * N];	

	double amd =  a - d;

	double mup = ((a + d)/(double) 2.0) + sqrt(b * b + pow((amd) / (double) 2.0,2));

	double r1 = a - mup;
	double r2 = b;
	
	double rp1 = r2;
	double rp2 = r1;

	double modr_inv = 1 / (double)hypot(r1,r2);

	double c = modr_inv * rp1;
	double s = -1 * modr_inv * rp2;

	//Return Rotation Values
	C = c;
	S = s;
}


//________________________________________________________________________//
//      *************************************************************     //
//                  Carry Out the Jacobi Algorithm
//      Matrix ordering is row major format. Calculated eigenvectors
//	and eigenvalues of a given matrix through iterative
//	applications of a rotation matrix.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/

extern void Jacobi_Algorithm(double *Evec_out,double *Eval_out,double *S_in,int N,double *converg_data)
{
	//Allocate memory
	double *Awork;
        Awork = new double [N * N];

	int *max; 
	max = new int [2];

	//Copy S_ini to Awork
	memcpy(Awork,S_in,N * N * sizeof(double));

	/*Awork[0] = 0.58464; Awork[1] = 0.0505; Awork[2] = 0.6289; Awork[3] = 0.2652; Awork[4] = 0.6857;
	Awork[5] = 0.0505; Awork[6] = 0.19659; Awork[7] = 0.2204; Awork[8] = 0.3552; Awork[9] = 0.0088;
	Awork[10] = 0.6289; Awork[11] = 0.2204; Awork[12] = 0.44907; Awork[13] = 0.1831; Awork[14] = 0.5086;
	Awork[15] = 0.2652; Awork[16] = 0.3552; Awork[17] = 0.1831; Awork[18] = 0.21333; Awork[19] = 0.272;
	Awork[20] = 0.6857; Awork[21] = 0.0088; Awork[22] = 0.5086; Awork[23] = 0.272; Awork[24] = 0.49667;
	*/
	double eps = 1;
	int it = 0;
	while(eps > converg_data[0])
	{	
 
		find_max(Awork,max,N);

		int i = max[0]; int j = max[1];
		double c; double s;

		//cout << "----------i: " << i << " j: " << j;

		Make_Rot_mat(c,s,Awork,i,j,N);
		//cout << " c: " << c << " s: " << s << "-----------\n";

                /*cout << "|------------------- (Ain) ------------------|\n";
		cout.setf( std::ios::fixed, std::ios::floatfield );
                for (int i = 0; i < N; ++i)
                {
                        for (int j = 0; j < N; ++j)
                        {
                                cout << setprecision(14) << " " << Awork[j + i * N] << "  ";
                        }
                        cout << "\n";
                }*/
	
		JacobiMult_Special(Awork,c,s,i,j,N);

	        /*cout << "|------------------- (Aout) ------------------|\n";
	        for (int i = 0; i < N; ++i)
	        {
	                for (int j = 0; j < N; ++j)
	                {
	                        cout << " " << Awork[j + i * N] << "  ";
	                }
	                cout << "\n";
	        }
		cout << "\n";*/
		Build_Evec(Evec_out,c,s,i,j,N,it);

		eps = check_convergence(Awork,N);
	
		++it;
		int max_it = converg_data[1];
		if (it == max_it) {break;}
	}

	/*cout << "|------------------- (A) ------------------|\n";
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			cout << " " << Awork[j + i * N] << "  ";
		}
		cout << "\n";
	}*/

 	//cout << "\nJacobi Iterations: " << it << " Final Convergence: " << check_convergence(Awork,N) << " \n ";
	converg_data[2] = it;	

	for (int n = 0; n < N; ++n)
        	{Eval_out[n] = Awork[n + n * N];}

	/*cout << "|------------------- (EvecOut) ------------------|\n";
        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        cout << " " << Evec_out[j + i * N] << "  ";
                }
                cout << "\n";
        }*/

	//Free memory
	delete [] Awork; Awork = NULL;
	delete [] max; max = NULL;
}

