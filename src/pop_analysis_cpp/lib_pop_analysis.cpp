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
#include "../scf_main/scf_classes.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                  Returns Working Info for Pop Analysis
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/16/2014  ----
                   ----Modified By:               ----
*/

int rtn_num_orbs(int AtomicNum)
{
	int orbs;
	switch (AtomicNum)
	{
		case 1: {orbs = 1;break;}
		case 2: {orbs = 1;break;}
		case 3: {orbs = 5;break;}
		case 4: {orbs = 5;break;}
		case 5: {orbs = 5;break;}
		case 6: {orbs = 5;break;}
		case 7: {orbs = 5;break;}
		case 8: {orbs = 5;break;}
		case 9: {orbs = 5;break;}
		case 10: {orbs = 5;break;}
	}
	return orbs;
}

string rtn_atom_lett(int AtomicNum)
{
	string AtomLett;
        switch (AtomicNum)
        {
                case 1: {AtomLett = "H";break;}
                case 2: {AtomLett = "He";break;}
                case 3: {AtomLett = "Li";break;}
                case 4: {AtomLett = "Be";break;}
                case 5: {AtomLett = "B";break;}
                case 6: {AtomLett = "C";break;}
                case 7: {AtomLett = "N";break;}
                case 8: {AtomLett = "O";break;}
                case 9: {AtomLett = "F";break;}
                case 10: {AtomLett = "Ne";break;}
        }
	return AtomLett;
}

int get_num_orbsB4(int idx,MemHandler *data_mem)
{
	int rtnval = 0;
	for (int k = 0; k < idx; ++k)
	{
		int Aidk = data_mem->atom_idx[k];
		int Anumk = data_mem->atom_data[k].atomic_num;
		int orbsk = rtn_num_orbs(Anumk);
		rtnval += orbsk;
	}
	return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                       Population Analysis Main
//      Main Function for Carrying out population analysis of the 
//	calculated density matrix form SCF.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 11/16/2014  ----
                   ----Modified By:               ----
*/

extern void population_analysis(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
        optfile->ofile << "|------------------------------------------------------|\n";
	optfile->ofile << "|---------------Begin Population Analysis--------------|\n";
        optfile->ofile << "|------------------------------------------------------|\n";

	int N = scf->orbitals;        
	int AN = data_mem->ipt_parms.num_atoms;

	optfile->Matrix_printer("Final Density Matrix",scf->Rho,N);
        double *PopMat; PopMat = new double [N * N];
	Multi_Elements(PopMat,scf->Rho,scf->S,N);
	optfile->Matrix_printer("Mulliken Population Matrix",PopMat,N);
	optfile->ofile << "\nOrbital Occupation Data:\n\n";

	int orbs_comp = 0;
	double MolPop = 0;
	double MolCharge = 0;
	for (int i = 0; i < AN; ++i)
	{
		double AtomPop = 0;

		int AtomID = data_mem->atom_idx[i + orbs_comp];
		int AtomicNum = data_mem->atom_data[AtomID].atomic_num;
		string AtomLett = rtn_atom_lett(AtomicNum);
		int orbs = rtn_num_orbs(AtomicNum);
		//cout << " ATOMID: " << AtomID << " AtomicNum: " << AtomicNum << " AtomLett: " << AtomLett << " orbs: " << orbs << "\n";

		optfile->ofile << "|==== Atom(" << i << "): " << AtomicNum << "-" << AtomLett << " Number of Orbitals: " << orbs << " ====|\n\n";
		for (int j = 0; j < orbs; ++j)
		{
			int col_idx = j + orbs_comp + i;
			//cout << " Colidx: " << col_idx << "\n";
			string orbtype = data_mem->basis_type(data_mem->orb_idx[col_idx]);
			double pop_val = 0;
			
			for (int m = 0; m < N; ++m)
				{pop_val += PopMat[col_idx + m * N];}
			
			optfile->ofile << "   Orbital(" <<  col_idx << "): " << orbtype << " Population: " << pop_val << "\n";
			AtomPop += pop_val;
		}

		optfile->ofile << "\n   Total Atom Population: " <<  AtomPop << "\n";
		optfile->ofile << "   Atom Mulliken Charge: " <<  AtomicNum - AtomPop << "\n\n";
		
		orbs_comp += orbs - 1;
		MolCharge += AtomicNum - AtomPop;
		MolPop += AtomPop;
	}

	optfile->ofile << "Total System Population: " << MolPop << "\n";
	optfile->ofile << "Total System Charge: " << MolCharge <<  "\n\n";

	optfile->ofile << "|------------------------------------------------------|\n";
	optfile->ofile << "|----------------End Population Analysis---------------|\n";
        optfile->ofile << "|------------------------------------------------------|\n\n";

	optfile->ofile << "Producing Mayer Bond Order Matrix... \n";
        
        double *BO; BO = new double [AN * AN];
        double *PS; PS = new double [N * N];
        Matrix_Mult(PS,scf->Rho,0,scf->S,0,N);

        optfile->Matrix_printer("PS Matrix",PS,N);

	Null_Set(BO,AN);

	optfile->ofile << "\nBond Order Matrix (Aprroximate to one decimal): ";
	optfile->ofile << "\n     [ ";
	for (int i = 1; i < AN; ++i)
		{optfile->ofile << "(" << i+1 << ") | ";}
	optfile->ofile << "\n";

	for (int i = 0; i < AN - 1; ++i)
	{
		optfile->ofile << " (" << i+1 << ") [ ";

		for (int j = i + 1; j --> 1; )
                	{optfile->ofile << "    | ";}

		for (int j = i + 1; j < AN; ++j)
		{
			int AidA = data_mem->atom_idx[i];
			int AidB = data_mem->atom_idx[j];

                	int AnumA = data_mem->atom_data[i].atomic_num;
                	int AnumB = data_mem->atom_data[j].atomic_num;

			int orbsA = rtn_num_orbs(AnumA);
			int orbsB = rtn_num_orbs(AnumB);

			int startA = get_num_orbsB4(i,data_mem);
			int startB = get_num_orbsB4(j,data_mem);

			//cout << i << "(" << orbsA << "-" << AnumA << "):" << startA << " -> " << j << "(" << orbsB << "-" << AnumB << "):" << startB << "  |  ";

			for (int k = startA; k < startA + orbsA; ++k)
			{
				for (int l = startB; l < startB + orbsB; ++l)
				{
					//cout << k << "*" << l << "|";
					BO[i + j * AN] += PS[k + l * N] * PS[l + k * N];
					//BO[i + j * AN] += scf->Rho[k + l * N] * scf->S[k + l * N] * scf->Rho[l + k * N] * scf->S[l + k * N];
				}
			}
			optfile->Prec_Printer(BO[i + j * AN],1);
			optfile->ofile << " | ";
			//cout << "\n";
		} 
		optfile->ofile << "\n";
		//cout << "\n";
	}
        
	//optfile->Matrix_printer("BO Matrix",BO,AN);

	optfile->ofile << "\n";
	delete [] PopMat; PopMat = NULL;
	delete [] PS; PS = NULL;
	delete [] BO; BO = NULL;
}
