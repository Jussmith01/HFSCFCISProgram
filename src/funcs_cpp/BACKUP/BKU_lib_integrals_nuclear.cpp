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
//            Calculate the MuTerms factor from the Nuc Att Solution
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/

double MuTerms(int ox1,int ox2,int rx,int ux,int mux,int ax,int bx, double Ax,double Bx,double Cx,double an, double bm)
{
        double rtnval;

        double Px = (an * Ax + bm * Bx) / (double)(an + bm);

        rtnval = ((pow(-1,ux) * fact(mux) * pow(Px - Cx,mux - 2 * ux)) / (double)(pow(4,ux) * fact(ux) * fact(mux - 2 * ux) * pow(an + bm,ox1 + ox2 - rx + ux)));

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//            Calculate the A factor from the Nuc Att Solution
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/
double A(int ix1,int ix2,int ox1,int ox2,int rx,int ax,int bx, double Ax,double Bx,double an, double bm)
{
        double rtnval;

	double fact1;
	double fact2;

        fact1 = ((pow(-1,ox2 + rx) * fact(ox1 + ox2)) / (double)(pow(4,ix1 + ix2 + rx) * fact(ix1) * fact(ix2) * fact(ox1) * fact(ox2) * fact(rx)));
	fact2 = ((pow(an,ox2 - ix1 - rx) * pow(bm,ox1 - ix2 - rx) * pow(Ax - Bx,ox1 + ox2 - 2 * rx)) / (double)(fact(ax - 2 * ix1 - ox1) * fact(bx - 2 * ix2 - ox2) * fact(ox1 + ox2 - 2 * rx)));

	//cout << "          Fact1: " << fact1 << " Fact2: " << fact2 << "\n";

	rtnval = fact1 * fact2;

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
	double tol = 1.0E-14;
	double pos = 0;
	double neg = 0;

        double Px = (an * Ax + bm * Bx) / (double)(an + bm);
        double Py = (an * Ay + bm * By) / (double)(an + bm);
        double Pz = (an * Az + bm * Bz) / (double)(an + bm);
        
        double PmC = (Px - Cx) * (Px - Cx) + (Py - Cy) * (Py - Cy) + (Pz - Cz) * (Pz - Cz);

        //cout << "PmC: " << PmC << " an: " << an << " bm:" << bm << "\n";
        //cout << "Cx: " << Cx << " Cy: " << Cy << " Cz: " << Cz << " Px: " << Px  << " Py: " << Py << " Pz: " << Pz << "\n";
        //cout << "Ax: " << Ax << " Ay: " << Ay << " Az: " << Az << " Bx: " << Bx  << " By: " << By << " Bz: " << Bz << "\n\n";

        //Indicies for X calculations
	double Aconstx = pow(-1,ax + bx) * fact(ax) * fact(bx);
        for (int ix1 = 0; ix1 <= floor(ax / (double)2.0); ++ix1) {
        for (int ix2 = 0; ix2 <= floor(bx / (double)2.0); ++ix2) {
        for (int ox1 = 0; ox1 <= ax - 2.0 * ix1; ++ox1) {
        for (int ox2 = 0; ox2 <= bx - 2.0 * ix2; ++ox2) {
        for (int rx = 0; rx <= floor((ox1 + ox2) / (double)2.0); ++rx) {

	double Avalx = A(ix1,ix2,ox1,ox2,rx,ax,bx,Ax,Bx,an,bm);

	int mux = ax + bx - 2 * (ix1 + ix2) - (ox1 + ox2);
        for (int ux = 0; ux <= floor(mux / (double)2.0); ++ux) {

	double Uvalx = MuTerms(ox1,ox2,rx,ux,mux,ax,bx,Ax,Bx,Cx,an,bm);

                //Indicies for Y calculations
		double Aconsty = pow(-1,ay + by) * fact(ay) * fact(by);
        	for (int iy1 = 0; iy1 <= floor(ay / (double)2.0); ++iy1) {
        	for (int iy2 = 0; iy2 <= floor(by / (double)2.0); ++iy2) {
        	for (int oy1 = 0; oy1 <= ay - 2.0 * iy1; ++oy1) {
        	for (int oy2 = 0; oy2 <= by - 2.0 * iy2; ++oy2) {
        	for (int ry = 0; ry <= floor((oy1 + oy2) / (double)2.0); ++ry) {

		double Avaly = A(iy1,iy2,oy1,oy2,ry,ay,by,Ay,By,an,bm);
        	
		int muy = ay + by - 2 * (iy1 + iy2) - (oy1 + oy2); 
        	for (int uy = 0; uy <= floor(muy / (double)2.0); ++uy) {
		
		double Uvaly = MuTerms(oy1,oy2,ry,uy,muy,ay,by,Ay,By,Cy,an,bm);

                        //Indicies for Z calculations
			double Aconstz = pow(-1,az + bz) * fact(az) * fact(bz);
                	for (int iz1 = 0; iz1 <= floor(az / (double)2.0); ++iz1) {
                	for (int iz2 = 0; iz2 <= floor(bz / (double)2.0); ++iz2) {
                	for (int oz1 = 0; oz1 <= az - 2.0 * iz1; ++oz1) {
                	for (int oz2 = 0; oz2 <= bz - 2.0 * iz2; ++oz2) {
                	for (int rz = 0; rz <= floor((oz1 + oz2) / (double)2.0); ++rz) {

			double Avalz = A(iz1,iz2,oz1,oz2,rz,az,bz,Az,Bz,an,bm);
                	
			int muz = az + bz - 2 * (iz1 + iz2) - (oz1 + oz2);
                	for (int uz = 0; uz <= floor(muz / (double)2.0); ++uz) {

			double Uvalz = MuTerms(oz1,oz2,rz,uz,muz,az,bz,Az,Bz,Cz,an,bm);

                                double Atotx = Aconstx * Avalx * Uvalx; 
                                double Atoty = Aconsty * Avaly * Uvaly; 
                                double Atotz = Aconstz * Avalz * Uvalz; 
               
				//cout << "      powz: " << pow(-1,az + bz) << " factaz: " << fact(az) << " factbz: " << fact(bz) << " A: " << A(iz1,iz2,oz1,oz2,rz,uz,muz,az,bz,Az,Bz,Cz,an,bm) << "\n"; 
	
                                int nu = mux + muy + muz - (ux + uy + uz);
				double val = Atotx * Atoty * Atotz * 2.0 * Fnu(nu,PmC * (an + bm));
                                //cout << "    Atotx: " << Atotx << " Atoty: " << Atoty << " Atotz: " << Atotz << " Avalx: " << Avalx << " Avaly: " << Avaly << " Atotz: " << Avalz << " Uvalx: " << Uvalx << " Uvaly: " << Uvaly << " Utotz: " << Uvalz << "\n";

				if (val > 0 && val > tol)
				{
					pos += val;
				} else if (val < 0 && abs(val) > tol) {
					neg += abs(val);
				}

                                //rtnval += Avalx * Avaly * Avalz * 2;
                                //cout << "      VAL: " << val << " Avalx: " << Avalx << " Avaly: " << Avaly << " Avalz: " << Avalz << " Fnu: " << Fnu(nu,PmC * (an + bm)) << "\n";

        		}}}}}}
        	}}}}}}
        }}}}}}
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

			rtnval += norm2(ax,ay,az,an) * norm2(bx,by,bz,bm) * dn * dm * Ifact * C_fact(an,bm) * E_factor(an,bm,AmB);
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

