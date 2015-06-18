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
//                      Overlap Analytical Solution
//      Analytical Solution to the STOnG basis set for the overlap
//	values.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/17/2014  ----
                   ----Modified By:               ----
*/

/*extern double S_factor(int axyz,int bxyz,double an,double bm,double Axyz,double Bxyz)
{
	double rtnval = 0;

	double Pxyz =  (an * Axyz + bm * Bxyz) / (double)(an + bm);
	
	for (int i = 0; i <= axyz ;++i)
        {
                for (int j = 0; j <= bxyz; ++j)
                {
			if ((i + j) == 0 || (i + j) % 2 == 0)
			{
				rtnval += nCk(axyz,i) * nCk(bxyz,j) * pow(Pxyz - Axyz,(double)(axyz - i)) * pow(Pxyz - Bxyz,(double)(bxyz - j)) * (d_fact(i + j - 1)/(double)pow(2 * (an + bm),(i + j)/(double)(2)));
			}
		}
	}
	rtnval = rtnval * sqrt(M_PI / (double)(an + bm));

	return rtnval;
}*/

//________________________________________________________________________//
//      *************************************************************     //
//                     Overlap Analytical Solution 2
//      Analytical Solution to the STOnG basis set for the overlap
//      values.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/17/2014  ----
                   ----Modified By:               ----
*/
extern double S_factor(int axyz,int bxyz,double an,double bm,double Axyz,double Bxyz)
{
        double rtnval = 0;

	for (int i1 = 0; i1 <= floor(0.5 * axyz) ;++i1)
        {
                for (int i2 = 0; i2 <= floor(0.5 * bxyz) ;++i2)
                {
        		double omega = axyz + bxyz - 2 * (i1 + i2);
                        for (int o = 0; o <= floor(0.5 * omega); ++o)
			{
                                rtnval += ((pow(-1,o) * fact(omega) * pow(an,bxyz - i1 - 2 * i2 - o) * pow(bm,axyz - 2 * i1 - i2 - o)) / (double)(pow(4,i1 + i2 + o) * fact(i1) * fact(i2) * fact(o))) * ((pow(an + bm,2 * (i1 + i2) + o) * pow(Axyz - Bxyz,omega - 2 * o)) / (double)(fact(axyz - 2 * i1) * fact(bxyz - 2 * i2) * fact(omega -  2 * o)));
                        }
                }
        }
        rtnval = ((pow(-1,axyz) * fact(axyz) * fact(bxyz)) / (double)(pow(an + bm,axyz + bxyz))) * rtnval;

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                        Calculate the Overlap
//      Function that calculates the overlap value with given inputs.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/17/2014  ----
                   ----Modified By:               ----
Math functions nCk(), fact(), d_fact()
*/

double calculate_overlap(int basis1,int basis2,MemHandler *data_mem,scf_data *scf,STOnG_handler *STO)
{
        double rtnval = 0;
	int N = STO->N;
	
	int basis_row = data_mem->orb_idx[basis1]; //orbital type 1,2,3,4,5 ...
	int basis_col = data_mem->orb_idx[basis2];

	int atomic_data_row = data_mem->atom_idx[basis1]; //Gives index of the atom in question
	int atomic_data_col = data_mem->atom_idx[basis2];

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
	
	//(B-A)(B-A)
	double AmB = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
	
	//Used for indexing the constants array and obtaining the correct constants
	int basis_idx_row = STO[Anum].rtn_bas_idx(basis_row);
	int basis_idx_col = STO[Bnum].rtn_bas_idx(basis_col);

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

	//cout << "basis1: " << basis1 << " basis2: " << basis2 << " Anum: " << Anum << " Bnum: " << Bnum <<"\n";
	for (int n = 0; n < N ;++n)
	{	
		for (int m = 0; m < N ;++m)
		{
			double an = STO[Anum - 1].a[n + (N * basis_idx_row - N)];
			double bm = STO[Bnum - 1].a[m + (N * basis_idx_col - N)];
			double dn = STO[Anum - 1].d[n + (N * basis_idx_row - N)];
			double dm = STO[Bnum - 1].d[m + (N * basis_idx_col - N)];
			//cout << "   n: " << n << " m: " << m << "\n";
			//cout << "      Ax: " << Ax << " Ay: " << Ay << " Az: " << Az << "\n";
			//cout << "      ax: " << ax << " ay: " << ay << " az: " << az << "\n";
			//cout << "      Bx: " << Bx << " By: " << By << " Bz: " << Bz << "\n";
			//cout << "      bx: " << bx << " by: " << by << " bz: " << bz << "\n";
			//cout << "      an: " << an << " dn: " << dn << " bm: " << bm << " dm: " << dm << "\n";
        		double val = norm2(ax,ay,az,an) * norm2(bx,by,bz,bm) * dn * dm * S_factor(ax,bx,an,bm,Ax,Bx) * S_factor(ay,by,an,bm,Ay,By) * S_factor(az,bz,an,bm,Az,Bz) * E_factor(an,bm,AmB) * pow(M_PI/(double)(an + bm),1.5);
			rtnval += val;
			//cout << "      rtnval: " << rtnval << " val: " << val << " norma: " << norm2(ax,ay,az,an) << " normb: " << norm2(bx,by,bz,bm) << " dn: " << dn << " dm: " << dm << " Sx: " << S_factor(ax,bx,an,bm,Ax,Bx) << " Sy: " << S_factor(ay,by,an,bm,Ay,By) << " Sz: " << S_factor(az,bz,an,bm,Az,Bz) << " Efactor: " << E_factor(an,bm,BmA) << " lastterm: " << pow(M_PI/(double)(an + bm),1.5) << "\n";
		}
	}
        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                        Overlap Main Function
//      Fills the overlap (S) matrix in the scf class.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/15/2014  ----
                   ----Modified By:               ----
*/

extern void overlap_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO)
{
	//Set some variables
	timer int_timer;
        optfile->ofile << "\nBeginning Overlap Calculations...\n";
	int_timer.set_timer();
	array_idx_funcs aif;//load 2d array manipulation class
	aif.width = scf->orbitals;

	//**************************
	//Start Overlap Computations	
	//*************************

	for (int row = 0; row < scf->orbitals; ++row) 
	{
		for (int col = 0; col < scf->orbitals; ++col)
		{
			//cout << "ROW: " << row << " COL: " << col << "\n";
			double Sij;
        		//Sij = calculate_overlap(row,col,data_mem,scf,STO);
        		
			if (row == col)
               			{Sij = 1.0;}
        		else
        			{Sij = calculate_overlap(row,col,data_mem,scf,STO);}

			aif.rpl_2dval(row,col,Sij,scf->S);
		}
	}
	
	//*************************
	//End Overlap Computations
	//*************************

	//Some ending functions
	int_timer.end_timer();
	string message = "Overlap Calculation Clock Time: ";
	int_timer.print_clock_time(message,optfile->ofile);
        optfile->ofile << "\n";
        //optfile->ofile << "Ending Overlap Calculations.\n\n";
}

