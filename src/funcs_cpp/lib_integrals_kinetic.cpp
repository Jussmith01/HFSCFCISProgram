#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <math.h>
#include <fstream>
#include <time.h>
#include <iostream>

#include "../utils_cpp/lib_includes.h"
#include "../classes/classes.h"
#include "../scf_main/scf_classes.h"
#include "func_classes.h"
#include "lib_integrals.h"

using namespace std;


//________________________________________________________________________//
//      *************************************************************     //
//                     Calculate the Kinetic Integrals
//      Function that calculates the kinetic contribution to the fock
//	matrix.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/20/2014  ----
                   ----Modified By:               ----
Math functions nCk(), fact(), d_fact()
*/

double calculate_kinetic(int basis1,int basis2,MemHandler *data_mem,scf_data *scf,STOnG_handler *STO)
{
        double rtnval = 0;
	int N = STO->N;
	
	//Orbital Index
	int basis_row = data_mem->orb_idx[basis1]; //orbital type 1,2,3,4,5 ...
	int basis_col = data_mem->orb_idx[basis2];

	//Atmoic Index
	int atomic_data_row = data_mem->atom_idx[basis1]; //Gives index of the atom in question
	int atomic_data_col = data_mem->atom_idx[basis2];

	//Atomic Number
	int Anum = data_mem->atom_data[atomic_data_row].atomic_num;
	int Bnum = data_mem->atom_data[atomic_data_col].atomic_num;

	//Row Position	
	double Ax = data_mem->atom_data[atomic_data_row].pos_xyz[0];
	double Ay = data_mem->atom_data[atomic_data_row].pos_xyz[1];
	double Az = data_mem->atom_data[atomic_data_row].pos_xyz[2];

	//Column Position
	double Bx = data_mem->atom_data[atomic_data_col].pos_xyz[0];
	double By = data_mem->atom_data[atomic_data_col].pos_xyz[1];
	double Bz = data_mem->atom_data[atomic_data_col].pos_xyz[2];
	
	//|A-B|^2
	double AmB = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
	
	//Used for indexing the constants array and obtaining the correct constants
	int basis_idx_row = STO[Anum].rtn_bas_idx(basis_row);
	int basis_idx_col = STO[Bnum].rtn_bas_idx(basis_col);

	//The following sets the angular numbers for the given orbital
	int ax,ay,az;
	int bx,by,bz;

	switch (basis_row)
	{
		case 1:{ax=0;ay=0;az=0;break;}
		case 2:{ax=0;ay=0;az=0;break;}
		case 3:{ax=1;ay=0;az=0;break;}
		case 4:{ax=0;ay=1;az=0;break;}
		case 5:{ax=0;ay=0;az=1;break;}
	} 
	
	switch (basis_col)
        {
                case 1:{bx=0;by=0;bz=0;break;}
                case 2:{bx=0;by=0;bz=0;break;}
                case 3:{bx=1;by=0;bz=0;break;}
                case 4:{bx=0;by=1;bz=0;break;}
                case 5:{bx=0;by=0;bz=1;break;}
        }

	//The following calculates the kinetic contribution of the given orbitals
	for (int n = 0; n < N ;++n)
	{	
		for (int m = 0; m < N ;++m)
		{
			double an = STO[Anum - 1].a[n + (N * basis_idx_row - N)];
			double bm = STO[Bnum - 1].a[m + (N * basis_idx_col - N)];
			double dn = STO[Anum - 1].d[n + (N * basis_idx_row - N)];
			double dm = STO[Bnum - 1].d[m + (N * basis_idx_col - N)];

			double Term1 = bm * (4.0 * (bx + by + bz) + 6.0) * S_factor(ax,bx,an,bm,Ax,Bx) * S_factor(ay,by,an,bm,Ay,By) * S_factor(az,bz,an,bm,Az,Bz);
			
			double Term2 = -4.0 * bm * bm * (S_factor(ax,bx + 2,an,bm,Ax,Bx) * S_factor(ay,by,an,bm,Ay,By) * S_factor(az,bz,an,bm,Az,Bz)
						       + S_factor(ax,bx,an,bm,Ax,Bx) * S_factor(ay,by + 2,an,bm,Ay,By) * S_factor(az,bz,an,bm,Az,Bz)
						       + S_factor(ax,bx,an,bm,Ax,Bx) * S_factor(ay,by,an,bm,Ay,By) * S_factor(az,bz + 2,an,bm,Az,Bz));

			double Term3 = -1.0 * (bx * (bx - 1) * S_factor(ax,bx - 2,an,bm,Ax,Bx) * S_factor(ay,by,an,bm,Ay,By) * S_factor(az,bz,an,bm,Az,Bz)
				             + by * (by - 1) * S_factor(ax,bx,an,bm,Ax,Bx) * S_factor(ay,by - 2,an,bm,Ay,By) * S_factor(az,bz,an,bm,Az,Bz)
				             + bz * (bz - 1) * S_factor(ax,bx,an,bm,Ax,Bx) * S_factor(ay,by,an,bm,Ay,By) * S_factor(az,bz - 2,an,bm,Az,Bz));

			double Term_sum = Term1 + Term2 + Term3;

			rtnval += 0.5 * norm(ax,ay,az,an) * norm(bx,by,bz,bm) * dn * dm * Term_sum * E_factor(an,bm,AmB) * pow(M_PI/(double)(an + bm),1.5);
		}
	}
        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                        Kinetic Main Function
//      Initial fill the fock (F) matrix in the scf class.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/20/2014  ----
                   ----Modified By:               ----
*/

extern void kinetic_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO)
{
	//Set some variables
	timer int_timer;
        optfile->ofile << "\nBeginning Kinetic Calculations...\n";
	int_timer.set_timer();
	array_idx_funcs aif;//load 2d array manipulation class
	aif.width = scf->orbitals;

	

	//**************************
	//Start Kinetic Computations	
	//*************************

	for (int row = 0; row < scf->orbitals; ++row) 
	{
		for (int col = 0; col < scf->orbitals; ++col)
		{
			//cout << "ROW: " << row << " COL: " << col << "\n";
			double zero = 0;
			aif.rpl_2dval(row,col,zero,scf->Hcore);
			
			double Tij;
        		Tij = calculate_kinetic(row,col,data_mem,scf,STO);
        		
			aif.rpl_2dval(row,col,Tij,scf->Hcore);
		}
	}
	
	//*************************
	//End Kinetic Computations
	//*************************

	//Som ending functions
	int_timer.end_timer();
        string message = "Kinetic Calculation Clock Time: ";
        int_timer.print_clock_time(message,optfile->ofile);
        optfile->ofile << "\n";
        //optfile->ofile << "Ending Overlap Calculations.\n\n";
}

