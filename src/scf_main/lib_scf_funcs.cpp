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
//                 Produce Density Matrix Given Coefficents
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/

extern void Density_Matrix(double *Rho,double *C,int *occ_orb,int N,int NO)
{
	double *Evec_in = new double [N * NO];

	double tol = 1.0E-14;

	for (int o = 0; o < NO; ++o)
	{
		int num_neg = 0;
		for (int i = 0; i < N; ++i)
		{
			Evec_in[o + i * NO] = C[occ_orb[o] + i * N];
		}
	}

        for (int i = 0; i < N; ++i)
        {
                for (int j = 0; j < N; ++j)
                {
			//cout << "------------(" << i << "," << j << ")--------------\n";
			double pos = 0;
			double neg = 0;

                	for (int o = 0; o < NO; ++o)
			{
				double val = Evec_in[o + i * NO] * Evec_in[o + j * NO];

				if (val > 0 && val > tol)
					{pos += val;}
				else if (val < 0 && abs(val) > tol)
					{neg += abs(val);}
				//cout << "C[" << i << "," << o <<  "]*C[" << j << ","  << o << "] pos = " << pos << " neg = " << neg <<  " | ";
			}

			Rho[j + i * N] = (double) 2.0 * (pos - neg);
		}
	}

	delete [] Evec_in; Evec_in = NULL;	
}

//________________________________________________________________________//
//      *************************************************************     //
//                 Find Minimum Value from a Vector
//      Returns memory index
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/

extern void find_mins(int *min,double *E,int N,int NO)
{
	for (int o = 0; o < NO; ++o)
	{
		double min_val = 1.0E10;
		for (int i = 0; i < N; ++i)
		{
			if (o == 0)
			{
				if (E[i] < min_val)
				{
					min[o] = i;
					min_val = E[i];
				}
			} else {
				int chk = 0;
				for (int o2 = 0; o2 < o; ++o2)
				{
					if (min[o2] == i)
						{chk = 1;}
				}
				if (chk == 0 && E[i] < min_val)
					{
						min[o] = i;
						min_val = E[i];
					}
			}
		}
	}
}

//________________________________________________________________________//
//      *************************************************************     //
//                 Produce Density Matrix Occupied Orbitals
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/
extern void Produce_Density(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	int N = scf->orbitals;
	int NO = data_mem->num_elec / 2;
	//********************************************//
        //   Find Minimum Energy Orbital Occupancy    //
        //********************************************//
	
	int *min = new int [NO];
	find_mins(min,scf->E,N,NO);
	
	optfile->ofile << "Minimum Energy Eigenvalues: ";
	for (int i = 0; i < NO; ++i)
	{
		optfile->ofile << min[i];
	
		if (i != NO - 1)
		{
			optfile->ofile << ", ";
		}
	}
	optfile->ofile << "\n";

        //********************************************//
        //  Produce Density Matrix from Coefficients  //
        //********************************************//
	Density_Matrix(scf->Rho,scf->C,min,N,NO);
	delete [] min; min = NULL;
}

//________________________________________________________________________//
//      *************************************************************     //
//               Produce Electronic Terms of the Fock Matrix
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/
extern void Produce_Electronic_Terms(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
        timer ET_timer;
        ET_timer.set_timer();

	double tol = 1.0E-14;
        int N = scf->orbitals;
	memcpy(scf->F,scf->Hcore,N * N * sizeof(double));

        //********************************************//
        //   Find Electronic Terms then Add to FOCK   //
        //********************************************//
        for (int i = 0; i < N; ++i) {
        	for (int j = 0; j < N; ++j) {//i and j gives the fock posistions
			
			double pos = 0;
			double neg = 0;
			int count = 0;
			//cout << "i: " << i << " j: " << j << "\n";
        		
			for (int k = 0; k < N; ++k) {
        			for (int l = 0; l < N; ++l) {//k and l run through all other orbitals

					double density = scf->Rho[l + k * N];

					if (abs(density) > tol)
					{
		        			//Orbital Index
        					//int basis_A = data_mem->orb_idx[i]; //orbital type 1,2,3,4,5 ...
        					//int basis_B = data_mem->orb_idx[j];
        					//int basis_C = data_mem->orb_idx[k];
        					//int basis_D = data_mem->orb_idx[l];

        					//Atmoic Index
				  	        //int atomic_data_A = data_mem->atom_idx[i]; //Gives index of the atom in question
 				 	        //int atomic_data_B = data_mem->atom_idx[j];
					        //int atomic_data_C = data_mem->atom_idx[k];
					        //int atomic_data_D = data_mem->atom_idx[l];

						//int ortho_chk = 0;

						double E1; double E2;

						//This statement decides if the basis functions are orthogonal P functions, if they are set rtnval to 0.
        					/*if ((basis_A >= 3 && basis_B >= 3 && atomic_data_A == atomic_data_B && basis_A != basis_B) || (basis_C >= 3 && basis_D >= 3 && atomic_data_C == atomic_data_D && basis_C != basis_D))
        						{E1 = 0.00;++ortho_chk;}
						else*/
							{E1 = scf->FEI.get_val(i,j,k,l);}

						/*if ((basis_A >= 3 && basis_C >= 3 && atomic_data_A == atomic_data_C && basis_A != basis_C) || (basis_B >= 3 && basis_D >= 3 && atomic_data_B == atomic_data_D && basis_B != basis_D))
                                                        {E2 = 0.00;++ortho_chk;}
                                                else*/
                                                        {E2 = scf->FEI.get_val(i,k,j,l);}

						double val = density * (E1 - 0.5 * E2);

						//if (ortho_chk != 2)
						//{
							if (val > 0 && abs(val) > tol)
								{pos += val;} 
							else if (val < 0 && abs(val) > tol) 
								{neg += abs(val);}
						//}

					} else {++count;}

					//cout << "A: " << i << " A: " << j << " A: " << k << " A: " << l << " Unique-IDENT: " << scf->FEI.produce_unique_ident(i,j,k,l) << " E1: " << E1 << "\n";
					//cout << "B: " << i << " B: " << k << " B: " << j << " B: " << l << " Unique-IDENT: " << scf->FEI.produce_unique_ident(i,k,j,l) << " E2: " << E2 << "\n";
				}
			}
			//cout << "COUNT SKIPPED: " << count << "\n";

			scf->F[j + i * N] += pos - neg;
			//scf->F[j + i * N] += 0;
		}
	}
        ET_timer.end_timer();
        string message = "Electronic Term Clock Time: ";
        ET_timer.print_clock_time(message,optfile->ofile);
}

//      *************************************************************     //
//                  Calculate the Electronic Energy
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/
extern double Calc_Electronic_Energy(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
        int N = scf->orbitals;
	double rtnval = 0;

	//cout << "TEST\n";
        //********************************************//
        //         Carryout Energy Summation          //
        //********************************************//
        for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {//i and j gives the fock positions

                        rtnval += scf->Rho[j + i * N] * (scf->Hcore[j + i * N] + scf->F[j + i * N]);
			//cout << "rtnval: " << rtnval << "\n";
			//cout << "Rho: " << scf->Rho[j + i * N] << " Hcore: " << scf->Hcore[j + i * N] << " F: " << scf->F[j + i * N] << "\n"; 

                }
        }

	rtnval = rtnval * 0.5;	

	return rtnval;
}

//      *************************************************************     //
//                 Produce Basic Initial Density Matrix
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/3/2014  ----
                   ----Modified By:               ----
*/
extern void Produce_Initial_Density(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
        timer den_timer;
        den_timer.set_timer();

        int N = scf->orbitals;
        int NE = data_mem->num_elec;

	Null_Set(scf->Rho,N);

	//cout << "NE: " << NE << "\n";

        //********************************************//
        // 		Fill Shells     	      //
        //********************************************//
	for (int i = 0; i < 3; ++i)
	{
		int NUMELEC = 2;

		if (data_mem->orb_idx[i] == 3)
		{
			NUMELEC = 6;
		}

		for (int j = 0; j < NUMELEC * N; ++j)
		{
			int k = j % N;
			//cout << "K: " << k << "\n";
			if (data_mem->orb_idx[k] == data_mem->orb_idx[i] && data_mem->orb_idx[i] <= 2)
			{
				scf->Rho[k + k * N] += 1.0;
				NE = NE - 1;
				//cout << "Filling Orbital: " << data_mem->orb_idx[i] << " NE: " << NE << "\n";
				if (NE == 0)
					{break;}
			} else if (data_mem->orb_idx[k] >= data_mem->orb_idx[i])
			{
                                scf->Rho[k + k * N] += 1.0;
                                NE = NE - 1;
                                //cout << "Filling Orbital: " << data_mem->orb_idx[i] << " NE: " << NE << "\n";
                                if (NE == 0)
                                        {break;}
			}
		}

		if (NE == 0)
			{break;}
	}
	//cout << "NE FINAL: " << NE << "\n";
	den_timer.end_timer();
        string message = "SCF Cycle Clock Time: ";
        den_timer.print_clock_time(message,optfile->ofile);
}

