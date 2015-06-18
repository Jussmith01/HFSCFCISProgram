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

//      *************************************************************     //
//                           Zero Machine Vals
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/30/2014  ----
                   ----Modified By:               ----
*/
extern void Multi_Elements(double *RtnMat,double *Matrix1,double *Matrix2,int N)
{
        double eps = 1.0E-14;
	
	Null_Set(RtnMat,N);
        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
			double val = Matrix1[j + i * N] * Matrix2[j + i * N];

                        if (abs(val) < eps)
                        	{val = 0.00;}

			RtnMat[j + i * N] = val;
                }
        }
}

//      *************************************************************     //
//                           Zero Machine Vals
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/30/2014  ----
                   ----Modified By:               ----
*/
extern void Set_EPS_d(double *A, int N)
{
        double eps = 1.0E-12;

        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        if (abs(A[j + i * N]) < eps)
                        {
                                A[j + i * N] = 0;
                        }
                }
        }
}

//________________________________________________________________________//
//      *************************************************************     //
//                           Zero Machine Vals
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/30/2014  ----
                   ----Modified By:               ----
*/
extern void Set_EPS(double *A, int N)
{
        double eps = 1.0E-14;

        for (int i = 0; i < N; ++i)
        {
        	for (int j = 0; j < N; ++j)
		{
                	if (abs(A[j + i * N]) < eps)
			{
				A[j + i * N] = 0;
			}
		}
        }
}

//________________________________________________________________________//
//      *************************************************************     //
//                           Trace a Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/30/2014  ----
                   ----Modified By:               ----
*/
extern double Trace(double *A, int N)
{
        double sum;

        for (int i = 0; i < N; ++i)
        {
        	sum += A[i + i * N];
        }

        return sum;
}

//________________________________________________________________________//
//      *************************************************************     //
//                           Sum a Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/30/2014  ----
                   ----Modified By:               ----
*/
extern double Sum_Mat(double *A, int N)
{
        double sum;

        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        sum += A[j + i * N];
                }
        }

	return sum;
}

//________________________________________________________________________//
//      *************************************************************     //
//                           Transpose a Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/30/2014  ----
                   ----Modified By:               ----
*/
extern void Transpose(double *A, int N)
{
	double *TempA;
        TempA = new double [N * N];

        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        TempA[i + j * N] = A[j + i * N];
                }
        }

	memcpy(A,TempA,N * N * sizeof(double));

        //Free Temp Matrix
        delete [] TempA;
        TempA = NULL;
}

//________________________________________________________________________//
//      *************************************************************     //
//                           Matrix Mult by Const
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/29/2014  ----
                   ----Modified By:               ----
*/
extern void Matrix_Mult_Const(double c,double *A, int N)
{
        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        A[j + i * N] = c * A[j + i * N];
                }
        }
}

//________________________________________________________________________//
//      *************************************************************     //
//                             Create Null Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/29/2014  ----
                   ----Modified By:               ----
*/
extern void Null_Set(double *A, int N)
{
        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        A[j + i * N] = 0;
                }
        }
}

//________________________________________________________________________//
//      *************************************************************     //
//                 Calculate the Determinant of a Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/28/2014  ----
                   ----Modified By:               ----
*/

extern double Determinant(double *A,int N)
{
        double det = 1;
	double rval;
	double dval;

        double *Matrix;
        Matrix = new double [N * N];

        memcpy(Matrix,A,N * N * sizeof(double));
        
	for (int i = 0; i < N - 1; ++i)
        {
		/*cout << "_______________(" << i << ")______________\n ";
                for (int n = 0; n < N; ++n)
                {
                        cout << " | ";
                        for (int m = 0; m < N; ++m)
                        {
                                cout << Matrix[m + n * N] << " ";
                        }
                        cout << " |\n";
                }*/

		rval = Matrix[i + i * N];
		//cout << "rval: " << rval << "\n";
		if (abs(rval) < 1.0E-12)
			{break;}
		for (int j = i + 1; j < N; ++j)
		{
			dval = Matrix[i + j * N];
			//cout << "dval: " << dval << "\n";

			for (int k = 0; k < N - i; ++k)
			{
				Matrix[i + k + j * N] = Matrix[i + k + j * N] + ((-1 * dval) / (double)(rval)) * Matrix[i + k + i * N];
				//cout << "Matrix Index" << i + k + j * N << " rindex: " << i + i * N << " dvalindex: " << i + j * N << "\n ";
			}	
		}
        }

	/*cout << "_______________(FINAL)______________\n ";
        for (int n = 0; n < N; ++n)
        {
                cout << " | ";
                for (int m = 0; m < N; ++m)
                {
                	cout << Matrix[m + n * N] << " ";
                }
        	cout << " |\n";
        }*/


	for (int i = 0; i < N; ++i)
	{
		det = det * Matrix[i + i * N];
	}

        //Free Temp Matrix
        delete [] Matrix;
        Matrix = NULL;

	return det;
}

//________________________________________________________________________//
//      *************************************************************     //
//                 Calculate the Cofactor Matrix of an Input
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/28/2014  ----
                   ----Modified By:               ----
*/

extern void Cofactor_Matrix(double *M_out, double *M_in,int N)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			double *SubMat;
        		SubMat = new double [(N) * (N)];

			int idxr;
			int idxc;

			for (int row = 0; row < N - 1; ++row)
			for (int col = 0; col < N - 1; ++col)
			{{
				if (row >= i)
					{idxr = row + 1;} else {idxr = row;}
                                if (col >= j)
                                        {idxc = col + 1;} else {idxc = col;}
				//cout << "i: " << i << " j: " << j << " row: " << row << " col: " << col << " idxr: " << idxr << " idxc: " << idxc << "\n";
				SubMat[col + row * (N - 1)] = M_in[idxc + idxr * N];
			}}

			M_out[j + i * N] = pow(-1,j + i) * Determinant(SubMat,N - 1);
			
			//cout << "M_out[" << j + i * N << "] = " << M_out[j + i * N] << "\n\n";

			delete [] SubMat;
       			SubMat = NULL;
		}	
	}
}


//________________________________________________________________________//
//      *************************************************************     //
//                 Calculate the Inverse Sqrt of a Matrix
//      Matrix ordering is row major format.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/

extern void Inverse_sqrt(double *Sinv,double *S,int N,dataOutput *optfile)
{
	//********************************************//
        //           Copy S to Temp Matrix            //
        //********************************************//
	double *Matrix;
	Matrix = new double [N * N];
	memcpy(Matrix,S,N * N * sizeof(double));

        double *MatrixSqr;
        MatrixSqr = new double [N * N];

        //********************************************//
        //     Calculate Square Root Temp Matrix      //
        //********************************************//
	double *EigenVecs;
        EigenVecs = new double [N * N];

	double *EigenVals;
        EigenVals = new double [N];

	double *convergence_data;
        convergence_data = new double [3];
	convergence_data[0] = 1.0E-30;
	convergence_data[1] = 8 * N * N;
	//cout << "STARTING S^-1/2 \n";
	Jacobi_Algorithm(EigenVecs,EigenVals,Matrix,N,convergence_data);
	optfile->ofile << "S^(-1/2) Jacobi Convergence: " << (int)convergence_data[2] << "\n";
	//cout << "ENDING S^-1/2: " << convergence_data[2] << " of " << convergence_data[1] << "\n\n";
	delete [] convergence_data; convergence_data = NULL;
	
	
	Null_Set(MatrixSqr,N);
	Null_Set(Matrix,N);

        for (int n = 0; n < N; ++n)
	{
		//MatrixSqr[n + n * N] = EigenVals[n];
		MatrixSqr[n + n * N] = 1 / (double)sqrt(EigenVals[n]);
		//cout << "InvSqr: " << MatrixSqr[n + n * N] << " Evals: " << EigenVals[n] << "\n";
	}

        Matrix_Mult(Matrix,EigenVecs,1,MatrixSqr,0,N);
        Matrix_Mult(MatrixSqr,Matrix,0,EigenVecs,0,N);

       //Free Eigen Matricies
        delete [] EigenVecs;
        EigenVecs = NULL;

        delete [] EigenVals;
        EigenVals = NULL;

	//********************************************//
        //       Calculate Determinant of Matrix      //
        //********************************************//
	//double det = Determinant(MatrixSqr,N);
	//cout << "DETERMINANT: " << det << "\n";

        //********************************************//
        //          Calculate Cofactors               //
        //********************************************//
	//Cofactor_Matrix(Matrix,MatrixSqr,N);
	//Matrix_Mult_Const(1 / (double)det, Matrix,N);
	//Transpose(Matrix, N);

        //********************************************//
        //           Copy Temp Matrix to S           //
        //********************************************//
        Set_EPS(MatrixSqr,N);
	memcpy(Sinv,MatrixSqr,N * N * sizeof(double));
       
	//Free Temp Matrix
        delete [] Matrix;
        Matrix = NULL;

        delete [] MatrixSqr;
        MatrixSqr = NULL;
}

//________________________________________________________________________//
//      *************************************************************     //
//                 Calculate B * A = C with double input
//      Matrix ordering is row major format. Trans = 1 makes transpose
//	of B.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/28/2014  ----
                   ----Modified By:               ----
*/

extern void Matrix_Mult(double *C,double *A,int transA,double *B,int transB,int N)
{
	double tol = 1.0E-13;
        
	//Allocate memory
        double *AD; AD = new double [N * N];
        double *BD; BD = new double [N * N];
        double *CD; CD = new double [N * N];

	for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        AD[j + i * N] = (double)A[j + i * N];
                        BD[j + i * N] = (double)B[j + i * N];
                        CD[j + i * N] = (double)C[j + i * N];
                }
        }

	//A TRANSPOSE CALC
	if (transB == 1)
	{
		for (int i = 0; i < N; ++i)
        	{
        	        for (int j = 0; j < N; ++j)
        	        {
                                double pos = 0;
                                double neg = 0;
                                //cout << "BTRANS ELEMENT i(" << i << ") j(" << j << "): ";

                                for (int k = 0; k < N; ++k)
                                {
                                        double tmpval = B[i + k * N] * A[j + k * N];
                                        if (tmpval > 0 && tmpval > tol)
                                        {
                                                pos += tmpval;
                                        } else if (tmpval < 0 && abs(tmpval) > tol) {
                                                neg += abs(tmpval);
                                        }
                                }
				double diff = pos - neg;
                                CD[j + i * N] = diff;
				//cout << " diff: " << diff << " CD: " << CD[j + i * N] << "\n";
                	}
        	}
	//A TRANSPOSE CALC
	} else if (transA == 1) {
                for (int i = 0; i < N; ++i)
                {
                        for (int j = 0; j < N; ++j)
                        {
                                double pos = 0;
                                double neg = 0;
                                //cout << "ATRANS ELEMENT i(" << i << ") j(" << j << "): ";

                                for (int k = 0; k < N; ++k)
                                {
                                        double tmpval = B[k + i * N] * A[k + j * N];
                                        if (tmpval > 0 && tmpval > tol)
                                        {
                                                pos += tmpval;
                                        } else if (tmpval < 0 && abs(tmpval) > tol) {
                                                neg += abs(tmpval);
                                        }
                                }

				double diff = pos - neg;
                                CD[j + i * N] = diff;
				//cout << " diff: " << diff << " CD: " << CD[j + i * N] << "\n";
                        }
                }	
	//NO TRANSPOSE	
	} else {
		for (int i = 0; i < N; ++i)
                {
                        for (int j = 0; j < N; ++j)
                        {
                                double pos = 0;
                                double neg = 0;
                                //cout << "REG ELEMENT i(" << i << ") j(" << j << ") ";

                                for (int k = 0; k < N; ++k)
                                {
                                        double tmpval = B[k + i * N] * A[j + k * N];
                                        if (tmpval > 0 && tmpval > tol)
                                        {
                                                pos += tmpval;
                                        } else if (tmpval < 0 && abs(tmpval) > tol) {
                                                neg += abs(tmpval);
                                        }
                                }
				
				double diff = pos - neg;
                                CD[j + i * N] = diff;
				//cout << " diff: " << diff << " CD: " << CD[j + i * N] << "\n";
                        }
                }

	}

	//Memcopy CD to C
	//cout << "-----" << "transA: " << transA<< " transB: " << transB << "------\n";
        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
                        C[j + i * N] = (double)CD[j + i * N];
			//cout << "MM C[j + i * N]: " << C[j + i * N] << " CD[j + i * N]: " << CD[j + i * N] << "\n";
                }
        }
	//cout << "-------------------------------------\n";

        //Free memory        
	delete [] AD; AD = NULL;
        delete [] BD; BD = NULL;
        delete [] CD; CD = NULL;
}

//________________________________________________________________________//
//      *************************************************************     //
// 	            Special Jacobi Matrix Multiplication
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/19/2014  ----
                   ----Modified By:               ----
NOTE: i and j are reversed in this, so works for j < i instead of i < j. I'm too lazy to change this.
*/
extern void JacobiMult_Special(double *A,double c,double s,int i,int j,int N)
{
        double *Awork;
        Awork = new double [N * N];

	for (int k = 0;k < N;++k)
	{
		for (int l = 0;l < N;++l)
		{
			double value;

                        if (k == i && l == i)
                                {value = A[i + i * N] * c * c + 2.0 * A[i + j * N] * c * s + A[j + j * N] * s * s;}

                        else if (k == j && l == j)
                                {value = A[i + i * N] * s * s - 2.0 * A[i + j * N] * c * s + A[j + j * N] * c * c;}

                        else if ((k == i && l == j) || (k == j && l == i))
                                {value = A[i + j * N] * (c * c - s * s) + s * c * (A[j + j * N] - A[i + i * N]);}

			else if (k == i && (l != j || l != i))
				{value = A[l + i * N] * c + A[l + j * N] * s;}
			
			else if ((k != j || k != i) && l == i)
				{value = A[k + i * N] * c + A[k + j * N] * s;}

			else if (k == j && (l != j || l != i))
				{value = A[l + j * N] * c - A[l + i * N] * s;}
			
			else if ((k != j || k != i) && l == j)
				{/*cout << " 7 ";*/value = A[k + j * N] * c - A[k + i * N] * s;}
                        
			//else if ((k != i && l != j) || (k != j && l != i))
                         else {value = A[l + k * N];}


			Awork[l + k * N] = value;
		}
		//cout << "\n";
	}
        memcpy(A,Awork,N * N * sizeof(double));
        delete [] Awork; Awork = NULL;
}

//________________________________________________________________________//
//      *************************************************************     //
//                    Special Jacobi Evec Builder
//	Same as the special jacobi matrix multiplication but left 
//	sided only.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/19/2014  ----
                   ----Modified By:               ----
*/
extern void Build_Evec(double *A,double c,double s,int i,int j,int N,int it)
{
        double *Awork;
        Awork = new double [N * N];

	if (it == 0)
	{
                for (int k = 0;k < N;++k)
                {
                        for (int l = 0;l < N;++l)
                        {
                                double value;

                                if (k == i && l == j)
                                        {value = -1 * s;}
                                else if (k == j && l == i)
                                        {value = s;}	
                                else if ((k == i && l == i) || (k == j && l == j))
                                        {value = c;}
                                else if (k == l)
                                        {value = 1.00;}
                                else
					{value = 0.00;}

                                Awork[l + k * N] = value;
                        }
                }
	} else {
        	for (int k = 0;k < N;++k)
        	{
                	for (int l = 0;l < N;++l)
                	{
        	                double value;
	
                        	if (l == i)
                        	        {value = A[i + k * N] * c + A[j + k * N] * s;}
                        	else if (l == j)
                        	        {value = A[j + k * N] * c - A[i + k * N] * s;}
                                else	{value = A[l + k * N];}

                        	Awork[l + k * N] = value;
               		}
        	}
	}
	memcpy(A,Awork,N * N * sizeof(double));
        delete [] Awork; Awork = NULL;
}

