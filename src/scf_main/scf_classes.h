#ifndef scf_includes
#define scf_includes
//C/C++ includes

//Custom includes
#include "../utils_cpp/lib_includes.h"
#include "../classes/classes.h"

// ********************************************************************* //
// ***************************DEFINE CLASSES**************************** //
// ********************************************************************* //

//________________________________________________________________________//
/*                 ----	      SCF Main Class 	  ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/14/2014  ----
                   ----Modified By:               ----			  */
//      *************************************************************     //
//                      Program Memory Handler Class
//      This class hold all needed memory handling declarations and
//	actions, such as memory allocations, transfers, to and from 
//	device and host, in the form of class member functions.
//      *************************************************************     //
class scf_data
{
	public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
		int orbitals; //Number of Orbitals
		int td_flag;
		long long int memory_req; //Stores the required system memory

		//Overlap Matrix
		double *S;
		
		//S inverse sqrt Matrix
		double *Sinv;

		//Hcore Matrix
		double *Hcore;
		
		//Fock Matrix
		double *F;

                //C constants Matrix
                double *C;

		//Eigen Values
                double *E;

		//Density Matrix
		double *Rho;
                
		//Working Matrix Matrix
                double *WM;

		//Transition Dipole Integerals
		double *TDx;
                double *TDy;
                double *TDz;
		
		//Electron Repulsion Integrals
		Four_Elec_Index FEI;
		//double *ERI;
		//string *ERI_idx;
		int num_deg;


        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
		
		//Allocated scf Memory
		void mem_alloc (dataOutput* optfile);

		//Free All Memory
		void mem_free(void);

		//Returns the main degenerate state
		//string ERI_deg_rtn(int i,int j,int k,int l);

		//Finds degenerate match and returns value
		//double ERI_Match(int i,int j,int k,int l);
};

// ********************************************************************* //
// ***************************DEFINE FUNCTIONS************************** //
// ********************************************************************* //
string scf_iteration_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);
void Produce_Density(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);
void Produce_Electronic_Terms(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);
double Calc_Electronic_Energy(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);
void Produce_Initial_Density(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);

void huckel_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);
void find_mins(int *min,double *E,int N,int NO);
void Density_Matrix(double *Rho,double *C,int *occ_orb,int N,int NO);

// ********************************************************************* //
// ***************************DEFINE LA FUNCTIONS*********************** //
// ********************************************************************* //
double Determinant(double *A,int N);
double Sum_Mat(double *A, int N);
void Inverse_sqrt(double *MatrixInv,double *Matrix,int N,dataOutput *optfile);
void Null_Set(double *A, int N);
void Cofactor_Matrix(double *M_out, double *M_in,int N);
void Matrix_Mult(double *C,double *A,int transA,double *B,int transB,int N);
void JacobiMult_Special(double *A,double c,double s,int i,int j,int N);
void Build_Evec(double *A,double c,double s,int i,int j,int N,int it);
void Matrix_Mult_Const(double c,double *A, int N);
void Transpose(double *A, int N);
void Make_Rot_mat(double &C,double &S,double *A,int i, int j,int N);
void Jacobi_Algorithm(double *Evec_out,double *Eval_out,double *S_in,int N,double *converg_data);
double Trace(double *A, int N);
void Set_EPS(double *A, int N);
void Set_EPS_d(double *A, int N);
void Multi_Elements(double *RtnMat,double *Matrix1,double *Matrix2,int N);

#endif
