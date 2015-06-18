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
//                        Obtain Valence Orbitals
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/21/2014  ----
                   ----Modified By:               ----
*/

void Rebuild_Density(double *Rho,double *Rho_val,int *val_orbs,int *core_orbs,int Nval,int Ncore,int N)
{
	for (int i = 0; i < Nval;++i)
        {
                for (int j = 0; j < Nval;++j)
                {
				Rho[val_orbs[j] + val_orbs[i] * N] = Rho_val[j + i * Nval];
                }
        }

        for (int i = 0; i < Ncore;++i)
        {
                                Rho[core_orbs[i] + core_orbs[i] * N] = 2.0000;
        }
}

//________________________________________________________________________//
//      *************************************************************     //
//                        Obtain Valence Orbitals
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/21/2014  ----
                   ----Modified By:               ----
*/

void Build_Heff(double *Heff,double *Sval,double *val_orbs,int Nval)
{
	double k = 1.75/(double)2;

        for (int i = 0; i < Nval;++i)
        {
                for (int j = 0; j < Nval;++j)
                {
			if (i == j)
				{Heff[j + i * Nval] = -val_orbs[i];}
			else
				{Heff[j + i * Nval] = -k * (val_orbs[i] + val_orbs[j]) * Sval[j + i * Nval];}
                }
        }
}

//________________________________________________________________________//
//      *************************************************************     //
//                        Overlap of Valence Orbitals
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/21/2014  ----
                   ----Modified By:               ----
*/

void Build_Sval(double *S,double *Sval,int *val_orbs_idx,int Nval,int N)
{
	for (int i = 0; i < Nval;++i)
	{
		for (int j = 0; j < Nval;++j)
		{
			Sval[j + i * Nval] = S[val_orbs_idx[j] + val_orbs_idx[i] * N];
		}
	}
}

//________________________________________________________________________//
//      *************************************************************     //
//                        Obtain Valence Orbitals
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/21/2014  ----
                   ----Modified By:               ----
*/

void obtain_valence_orbitals(double *&val_orbs,int *&val_orbs_idx,int *&core_orbs,int &Nval,int &Ncore,MemHandler *data_mem)
{
	int Na = data_mem->ipt_parms.num_atoms;

	Nval = 0;
	Ncore = 0;

	for (int i = 0; i < Na;++i)
	{
		int atomic_num = data_mem->atom_data[i].atomic_num;
		switch (atomic_num)
		{
			case 1: {Nval += 1;break;}
			case 2: {Nval += 1;break;}
			case 3: {Ncore += 1;Nval += 1;break;}
			case 4: {Ncore += 1;Nval += 1;break;}
			case 5: {Ncore += 1;Nval += 4;break;}
			case 6: {Ncore += 1;Nval += 4;break;}
			case 7: {Ncore += 1;Nval += 4;break;}
			case 8: {Ncore += 1;Nval += 4;break;}
			case 9: {Ncore += 1;Nval += 4;break;}
			case 10: {Ncore += 1;Nval += 4;break;}
		}
	}

	//cout << "test 1\n";

	val_orbs = new double [Nval];
	val_orbs_idx = new int [Nval];
	core_orbs = new int [Ncore];

	//cout << "test 2\n";
	int j = 0;
	int k = 0;
	int ei = 0;

	for (int i = 0; i < Na;++i)
        {
                int atomic_num = data_mem->atom_data[i].atomic_num;
                
		//cout << "I: " << i << " J: " << j << " K: " << k << " EI: " << ei << " Nval: " << Nval << " Ncore: " << Ncore << "\n";
		switch (atomic_num)
                {
                        case 1: {val_orbs[j] = 13.6;val_orbs_idx[j] = j + k + ei;j += 1;break;}
                        case 2: {val_orbs[j] = 24.5;val_orbs_idx[j] = j + k + ei;j += 1;break;}
                        case 3: {val_orbs[j] = 5.45;val_orbs_idx[j] = j + k + ei + 1;core_orbs[k] = j + k;j += 1;k += 1;ei += 3;break;}
                        case 4: {val_orbs[j] = 9.30;val_orbs_idx[j] = j + k + ei + 1;core_orbs[k] = j + k;j += 1;k += 1;ei += 3;break;}
                        case 5: {
                                        core_orbs[k] = j + k;
                                        val_orbs[j] = 14.0;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs[j + l] = 8.30;}

                                        core_orbs[k] = j + k;

					k += 1;
                                        val_orbs_idx[j] = j + k + ei;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs_idx[j + l] = j + k + l + ei;}

					j += 4;
                                        break;
                                }
                        case 6: {
                                        core_orbs[k] = j + k;
                                        val_orbs[j] = 19.5;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs[j + l] = 10.7;}

                                        core_orbs[k] = j + k;

                                        k += 1;
                                        val_orbs_idx[j] = j + k + ei;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs_idx[j + l] = j + k + l + ei;}

                                        j += 4;
                                        break;
                                }
                        case 7: {
                                        core_orbs[k] = j + k;
                                        val_orbs[j] = 25.5;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs[j + l] = 13.1;}

                                        core_orbs[k] = j + k;

                                        k += 1;
                                        val_orbs_idx[j] = j + k + ei;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs_idx[j + l] = j + k + l +ei;}

                                        j += 4;
                                        break;
                                }
                        case 8: {
                                        core_orbs[k] = j + k;
                                        val_orbs[j] = 32.3;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs[j + l] = 15.9;}

                                        core_orbs[k] = j + k;

                                        k += 1;
                                        val_orbs_idx[j] = j + k + ei;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs_idx[j + l] = j + k + l + ei;}

                                        j += 4;
                                        break;
                                }
                        case 9: {
                                        core_orbs[k] = j + k;
                                        val_orbs[j] = 46.4;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs[j + l] = 18.7;}

                                        core_orbs[k] = j + k;

                                        k += 1;
                                        val_orbs_idx[j] = j + k + ei;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs_idx[j + l] = j + k + l + ei;}

                                        j += 4;
                                        break;
                                 }
                        case 10: {
                                        core_orbs[k] = j + k;
                                        val_orbs[j] = 48.5;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs[j + l] = 21.5;}

                                        core_orbs[k] = j + k;

                                        k += 1;
                                        val_orbs_idx[j] = j + k + ei;
                                        for (int l = 1; l <= 3; ++l)
                                                {val_orbs_idx[j + l] = j + k + l + ei;}

                                        j += 4;
                                        break;
                                 }
                }
        }
	//cout << "TEST 3\n";
}

//________________________________________________________________________//
//      *************************************************************     //
//                     		Huckel Main
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/21/2014  ----
                   ----Modified By:               ----
*/

extern void huckel_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	optfile->ofile << "\n|===============Obtaining Huckel Orbitals================|\n";

	//cout << "BEGINNING HUCKEL TEST\n";
	int N = data_mem->orbs;
	double *val_orbs; int *val_orbs_idx; int Nval;
	int *core_orbs; int Ncore;

	obtain_valence_orbitals(val_orbs,val_orbs_idx,core_orbs,Nval,Ncore,data_mem);

	/*for (int i = 0; i < Nval; ++i)
		{cout << "(" << val_orbs_idx[i] << "): " << val_orbs[i] << " ";}
	cout << "\n";

	for (int i = 0; i < Ncore; ++i)
		{cout << "(" << i << "): " << core_orbs[i] << " ";}
	cout << "\n";*/

	double *Heff; Heff = new double [Nval * Nval];
        double *Sval; Sval = new double [Nval * Nval];

	Build_Sval(scf->S,Sval,val_orbs_idx,Nval,N);
	Build_Heff(Heff,Sval,val_orbs,Nval);

        double *Sval_inv; Sval_inv = new double [Nval * Nval];
        Inverse_sqrt(Sval_inv,Sval,Nval,optfile);

        Matrix_Mult(scf->WM,Heff,0,Sval_inv,0,Nval);//Calculate F = S^(-1/2) * F' * S(-1/2)
        Matrix_Mult(Heff,Sval_inv,0,scf->WM,0,Nval);

	double *Evecs; Evecs = new double [Nval * Nval];
	double *Evals; Evals = new double [Nval];

	double *converg_data;
        converg_data = new double [3];
        converg_data[0] = 1.0E-30;
        converg_data[1] = 8 * Nval * Nval;	
	Jacobi_Algorithm(Evecs,Evals,Heff,Nval,converg_data);
	optfile->ofile << "Heff' Jacobi Convergence: " << (int)converg_data[2] << " of " << (int)converg_data[1] << "\n";
        delete [] converg_data; converg_data = NULL;
	delete [] Heff; Heff = NULL;

	memcpy(Sval,Evecs,Nval * Nval * sizeof(double));
	Matrix_Mult(Evecs,Sval,0,Sval_inv,0,Nval);//Calculate C = S^(-1/2) * C'

	optfile->Matrix_printer("Huckel Eigen Vectors",Evecs,Nval);
        optfile->Array_printer("Huckel Eigen Values",Evals,Nval,4);

	int Val_Elec = data_mem->num_elec - 2 * Ncore;
	
	int NO = Val_Elec / 2;
	optfile->ofile << "Number of Valence Electrons: " << Val_Elec << " in " << NO << " orbitals.\n";
        //********************************************//
        //   Find Minimum Energy Orbital Occupancy    //
        //********************************************//

        int *min = new int [NO];
        find_mins(min,Evals,Nval,NO);

        optfile->ofile << "Huckel Minimum Energy Eigenvalues: ";
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
	double *Rho_val; Rho_val = new double [Nval * Nval];
        Density_Matrix(Rho_val,Evecs,min,Nval,NO);
        delete [] min; min = NULL;

	optfile->Matrix_printer("Huckel Density Matrix",Rho_val,Nval);

	Null_Set(scf->Rho,N);
	Rebuild_Density(scf->Rho,Rho_val,val_orbs_idx,core_orbs,Nval,Ncore,N);

	optfile->ofile << "\n|============End Huckel Orbitals Calculation===============|\n";

	delete [] Evecs; Evecs = NULL;
	delete [] Evals; Evals = NULL;
	delete [] Sval; Sval = NULL;
	delete [] Sval_inv; Sval_inv = NULL;
	delete [] val_orbs; val_orbs = NULL;
	delete [] val_orbs_idx; val_orbs_idx = NULL;
	delete [] core_orbs; core_orbs = NULL;
}
