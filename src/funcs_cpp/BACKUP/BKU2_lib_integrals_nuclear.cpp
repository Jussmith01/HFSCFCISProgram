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
//            Calculate the A factor from the Nuc Att Solution
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/
double A(int i,int r,int u,double C,double P,double ep)
{
        double rtnval;

	rtnval = (pow(-1,u) * fact(i) * pow(C - P,i - 2 * r - 2 * u) * pow(ep,r + u)) / (double)(fact(r) * fact(u) * fact(i - 2 * r - 2 * u));

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//      			Calc_A Terms 
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/

double calc_sum_A(int ax,int ay,int az,int bx,int by,int bz,double Ax,double Ay,double Az, double Bx,double By,double Bz,double Cx,double Cy,double Cz,double an,double bm)
{
        double rtnval;
	double tol = 1.0E-13;
	double pos = 0;
	double neg = 0;

        double Px = (an * Ax + bm * Bx) / (double)(an + bm);
        double Py = (an * Ay + bm * By) / (double)(an + bm);
        double Pz = (an * Az + bm * Bz) / (double)(an + bm);
        
        double PmC = (Px - Cx) * (Px - Cx) + (Py - Cy) * (Py - Cy) + (Pz - Cz) * (Pz - Cz);

	double ep = 1 / (double)(4 * (an + bm));

        //cout << "PmC: " << PmC << " an: " << an << " bm:" << bm << "\n";
        //cout << "Cx: " << Cx << " Cy: " << Cy << " Cz: " << Cz << " Px: " << Px  << " Py: " << Py << " Pz: " << Pz << "\n";
        //cout << "Ax: " << Ax << " Ay: " << Ay << " Az: " << Az << " Bx: " << Bx  << " By: " << By << " Bz: " << Bz << "\n\n";
	cout << "|--------------ax: " << ax << " bx: " << bx << "-----------|\n";
        //Indicies for X calculations
        for (int ix = 0; ix <= ax + bx; ++ix) {
        for (int rx = 0; rx <= floor(ix/(double)2); ++rx) {
        for (int ux = 0; ux <= floor((ix - 2 * rx)/(double)2); ++ux) {
	
	double Avalx = pow(-1,ix) * fi(ix,ax,bx,Px - Ax,Px - Bx) * A(ix,rx,ux,Cx,Px,ep);
	cout << "I: " << ix - 2 * rx - ux << " iru: " << ix << rx << ux <<"\n";
        	//Indicies for Y calculations
        	for (int iy = 0; iy <= ay + by; ++iy) {
        	for (int ry = 0; ry <= floor(iy/(double)2); ++ry) {
        	for (int uy = 0; uy <= floor((iy - 2 * ry)/(double)2); ++uy) {

        	double Avaly = pow(-1,iy) * fi(iy,ay,by,Py - Ay,Py - By) * A(iy,ry,uy,Cy,Py,ep);
		
                	//Indicies for Z calculations
                	for (int iz = 0; iz <= az + bz; ++iz) {
                	for (int rz = 0; rz <= floor(iz/(double)2); ++rz) {
                	for (int uz = 0; uz <= floor((iz - 2 * rz)/(double)2); ++uz) {

                	double Avalz = pow(-1,iz) * fi(iz,az,bz,Pz - Az,Pz - Bz) * A(iz,rz,uz,Cz,Pz,ep);
                
                                double nu = (ix + iy + iz) - 2.0 * (rx + ry + rz) - (ux + uy + uz);
				//cout << "nu: " << nu << "\n";
                                double val = Avalx * Avaly * Avalz * Fnu(nu,PmC / (double)(4 * ep));
                                if (val > 0 && val > tol)
                                        {pos += val;}
                                else if (val < 0 && abs(val) > tol)
                                        {neg += abs(val);}

                	}}}
		}}}
	}}}

	rtnval = pos - neg;

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                Calculate the Nuclear Attraction Integrals
//      Function that calculates the kinetic contribution to the fock
//	matrix.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
Math functions nCk(), fact(), d_fact()
*/

double calculate_nuclear(int basis1,int basis2,int ac,MemHandler *data_mem,scf_data *scf,STOnG_handler *STO)
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
	//Column Position
     	double Cx = data_mem->atom_data[ac].pos_xyz[0];
        double Cy = data_mem->atom_data[ac].pos_xyz[1];
        double Cz = data_mem->atom_data[ac].pos_xyz[2];

	double NucCharge = data_mem->atom_data[ac].atomic_num;
	
	//cout << "Cx: " << Cx << " Cy: " << Cy << " Cz: " << Cz << " NucCharge: " << NucCharge << "\n";
	//cout << "Ax: " << Ax << " Ay: " << Ay << " Az: " << Az << " Bx: " << Bx  << " By: " << By << " Bz: " << Bz << "\n\n";

	//cout << "BASIS SET -- basis1: " << basis1 << " basis2: " << basis2 << "\n";

	for (int n = 0; n < N ;++n)
	{	
		for (int m = 0; m < N ;++m)
		{
			//cout << "  PRIM GAUSS -- n: " << n << " m: " << m << "\n";
			
			double an = STO[Anum - 1].a[n + (N * basis_idx_row - N)];
			double bm = STO[Bnum - 1].a[m + (N * basis_idx_col - N)];
			double dn = STO[Anum - 1].d[n + (N * basis_idx_row - N)];
			double dm = STO[Bnum - 1].d[m + (N * basis_idx_col - N)];

			//cout << "BEFORE: an = " << an << " bm = " << bm << "\n";

			double Ifact = calc_sum_A(ax,ay,az,bx,by,bz,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,an,bm);

			rtnval += norm(ax,ay,az,an) * norm(bx,by,bz,bm) * dn * dm * Ifact * C_fact(an,bm) * E_factor(an,bm,AmB);
			//cout << "    Charge: " << NucCharge << " Norma: " << norm(ax,ay,az,an) << " Normb: " << norm(bx,by,bz,bm) << " Ifact: " << Ifact << " C_fact: " << C_fact(an,bm) << " E_factor: " << E_factor(an,bm,AmB) << "\n";
		}
	}
	rtnval = -1 * NucCharge * rtnval;
	//cout << "    RtnVal: " << rtnval << "\n";

       	return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                    Nuclear Attraction Main Function
//      Add the nuclear attaction contribution to the fock (F) matrix 
//	in the scf class.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/

extern void nuclear_main(int ac,MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO)
{
	//Set some variables
	timer int_timer;
        optfile->ofile << "\nBeginning Nuclear Attraction Calculations...\n";
	int_timer.set_timer();
	array_idx_funcs aif;//load 2d array manipulation class
	aif.width = scf->orbitals;

	//**************************************
	//Start Nuclear Attraction Computations	
	//**************************************

	for (int row = 0; row < scf->orbitals; ++row) 
	{
		for (int col = 0; col < scf->orbitals; ++col)
		{
			//cout << "ROW: " << row << " COL: " << col << "\n";
        		double Vij;
			Vij = calculate_nuclear(row,col,ac,data_mem,scf,STO);
			aif.rpl_2dval(row,col,Vij,scf->WM);

			double Hnewij = scf->Hcore[col + row * scf->orbitals] + Vij; 
			aif.rpl_2dval(row,col,Hnewij,scf->Hcore);
			//cout << "ROW: " << row << " COL: " << col << " Vij: " << Vij << "\n";	
		}
	}

	//************************************
	//End Nuclear Attraction Computations
	//************************************

	//Some ending functions
	int_timer.end_timer();
        string message = "Nuclear Attraction Calculation Clock Time: ";
        int_timer.print_clock_time(message,optfile->ofile);
        optfile->ofile << "\n";
}

