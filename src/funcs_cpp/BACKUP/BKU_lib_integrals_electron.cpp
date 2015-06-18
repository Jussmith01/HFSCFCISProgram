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
//           Calculate P0 portion of the Electron Integrals
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/
double J0_func(int axI,int axJ,double an1,double an2)
{
        double rtnval;

	rtnval = pow(-1,axI + axJ) * fact(axI) * fact(axJ) / (double)(pow(an1 + an2,axI + axJ));

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//           Calculate J1 portion of the Electron Integrals
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/
double J1_func(int i1,int i2,int o1,int o2,int r1,int aA,int aB,int aC,int aD,double an1,double an2,double an3,double an4,double Axyz,double Bxyz)
{
        double rtnval;

        rtnval = (pow(-1, o2 + r1) * fact(o1 + o2) / (double)(pow(4,i1 + i2 + r1) * fact(i1) * fact(i2) * fact(o1) * fact(o2) * fact(r1))) * ((pow(an1,o2 - i1 - r1) * pow(an2,o1 - i2 - r1) * pow(an1 + an2,2 * (i1 + i2) + r1) * pow(Axyz - Bxyz,o1 + o2 - 2 * r1)) / (double)(fact(aA - 2 * i1 - o1) * fact(aB - 2 * i2 - o2) * fact(o1 + o2 - 2 * r1))) * ((fact(aC) * fact(aD)) / (double)(pow(an3 + an4,aC + aD)));

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//           Calculate J2 portion of the Electron Integrals
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/
double J2_func(int i3,int i4,int o3,int o4,int r2,int aC,int aD,double an3,double an4,double Cxyz,double Dxyz)
{
        double rtnval;

        rtnval = ((pow(-1,o3 + r2) * fact(o3 + o4)) / (double)(pow(4,i3 + i4 + r2) * fact(i3) * fact(i4) * fact(o3) * fact(o4) * fact(r2))) * ((pow(an3,o4 - i3 - r2) * pow(an4,o3 - i4 - r2) * pow(an3 + an4,2 * (i3 + i4) + r2) * pow(Cxyz - Dxyz,o3 + o4 - 2 * r2)) / (double)(fact(aC - 2 * i3 - o3) * fact(aD - 2 * i4 - o4) * fact(o3 + o4 - 2 * r2)));

	//rtnval = 1;
        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//           Calculate J3 portion of the Electron Integrals
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/
double J3_func(int u,int mu,double n,double Pxyz,double Qxyz)
{
        double rtnval;

        rtnval = ((pow(-1,u) * fact(mu) * pow(n,mu - u) * pow(Pxyz - Qxyz,mu - 2 * u)) / (double)(pow(4,u) * fact(u) * fact(mu - 2 * u)));

	//rtnval = 1;
        return rtnval;
}


//________________________________________________________________________//
//      *************************************************************     //
//      			Calc_sum_J Terms 
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/

double calc_sum_J(int test,int axA,int ayA,int azA,int axB,int ayB,int azB,int axC,int ayC,int azC,int axD,int ayD,int azD,double Ax,double Ay,double Az,double Bx,double By,double Bz,double Cx,double Cy,double Cz,double Dx,double Dy,double Dz,double an1,double an2,double an3,double an4)
{
	double tol = 1.0E-13;
        double rtnval = 0;
	double pos = 0;
	double neg = 0;

        double Px = (an1 * Ax + an2 * Bx) / (double)(an1 + an2);
        double Py = (an1 * Ay + an2 * By) / (double)(an1 + an2);
        double Pz = (an1 * Az + an2 * Bz) / (double)(an1 + an2);
        
        double Qx = (an3 * Cx + an4 * Dx) / (double)(an3 + an4);
        double Qy = (an3 * Cy + an4 * Dy) / (double)(an3 + an4);
        double Qz = (an3 * Cz + an4 * Dz) / (double)(an3 + an4);
        
        double PmQ = (Px - Qx) * (Px - Qx) + (Py - Qy) * (Py - Qy) + (Pz - Qz) * (Pz - Qz);
	//double PmQ = ((an1 * pow(pow(Ax,2.0) + pow(Ay,2.0) + pow(Az,2.0),0.5) - (an2 * pow(pow(Bx,2.0) + pow(By,2.0) + pow(Bz,2.0),0.5)) / (double)(an1 + an2))) - ((an3 * pow(pow(Cx,2.0) + pow(Cy,2.0) + pow(Cz,2.0),0.5) - (an4 * pow(pow(Dx,2.0) + pow(Dy,2.0) + pow(Dz,2.0),0.5)) / (double)(an3 + an4)));
        //cout << "CmP: " << CmP << " an: " << an << " bm:" << bm << "\n";
        //cout << "Cx: " << Cx << " Cy: " << Cy << " Cz: " << Cz << " Px: " << Px  << " Py: " << Py << " Pz: " << Pz << "\n";
        cout << "axA: " << axA << " ayA: " << ayA << " azA: " << azA << " | " ;
        cout << " axB: " << axB << " ayB: " << ayB << " azB: " << azB << " | " ;
        cout << " axC: " << axC << " ayC: " << ayC << " azC: " << azC << " | " ;
        cout << " axD: " << axD << " ayD: " << ayD << " azD: " << azD << " | "  << " PmQ: " << PmQ << "\n";


	double nFACT = ((an1 + an2) * (an3 + an4)) / (double)((an1 + an2) + (an3 + an4));

        //Indicies for X calculations
	//double J0x = J0_func(axA,axB,an1,an2) * J0_func(axC,axD,an3,an4);
	//double J0y = J0_func(ayA,ayB,an1,an2) * J0_func(ayC,ayD,an3,an4);
	//double J0z = J0_func(azA,azB,an1,an2) * J0_func(azC,azD,an3,an4);
	double J0x = J0_func(axA,axB,an1,an2);
	double J0y = J0_func(ayA,ayB,an1,an2);
	double J0z = J0_func(azA,azB,an1,an2);
        for (int ix1 = 0; ix1 <= floor(axA / (double)2.0); ++ix1) {
        for (int ix2 = 0; ix2 <= floor(axB / (double)2.0); ++ix2) {
        for (int ox1 = 0; ox1 <= axA - 2.0 * ix1; ++ox1) {
        for (int ox2 = 0; ox2 <= axB - 2.0 * ix2; ++ox2) {
        for (int rx1 = 0; rx1 <= floor((ox1 + ox2) / (double)2.0); ++rx1) {

	double J1x = J1_func(ix1,ix2,ox2,ox2,rx1,axA,axB,axC,axD,an1,an2,an3,an4,Ax,Bx);
	
        for (int ix3 = 0; ix3 <= floor(axC / (double)2.0); ++ix3) {
        for (int ix4 = 0; ix4 <= floor(axD / (double)2.0); ++ix4) {
        for (int ox3 = 0; ox3 <= axC - 2.0 * ix3; ++ox3) {
        for (int ox4 = 0; ox4 <= axD - 2.0 * ix4; ++ox4) {
        for (int rx2 = 0; rx2 <= floor((ox3 + ox4) / (double)2.0); ++rx2) {

	double J2x = J2_func(ix3,ix4,ox3,ox4,rx2,axC,axD,an3,an4,Cx,Dx);

	int mux = axA + axB + axC + axD - 2.0 * (ix1 + ix2 + ix3 + ix4) - (ox1 + ox2 + ox3 + ox4);
        for (int ux = 0; ux <= floor(mux / (double)2.0); ++ux) {

	double J3x = J3_func(ux,mux,nFACT,Px,Qx);

        	//Indicies for Y calculations
        	//double J0y = J0_func(ayA,ayB,an1,an2);
        	for (int iy1 = 0; iy1 <= floor(ayA / (double)2.0); ++iy1) {
        	for (int iy2 = 0; iy2 <= floor(ayB / (double)2.0); ++iy2) {
        	for (int oy1 = 0; oy1 <= ayA - 2.0 * iy1; ++oy1) {
        	for (int oy2 = 0; oy2 <= ayB - 2.0 * iy2; ++oy2) {
        	for (int ry1 = 0; ry1 <= floor((oy1 + oy2) / (double)2.0); ++ry1) {

		double J1y = J1_func(iy1,iy2,oy2,oy2,ry1,ayA,ayB,ayC,ayD,an1,an2,an3,an4,Ay,By);

        	for (int iy3 = 0; iy3 <= floor(ayC / (double)2.0); ++iy3) {
        	for (int iy4 = 0; iy4 <= floor(ayD / (double)2.0); ++iy4) {
        	for (int oy3 = 0; oy3 <= ayC - 2.0 * iy3; ++oy3) {
        	for (int oy4 = 0; oy4 <= ayD - 2.0 * iy4; ++oy4) {
        	for (int ry2 = 0; ry2 <= floor((oy3 + oy4) / (double)2.0); ++ry2) {

		double J2y = J2_func(iy3,iy4,oy3,oy4,ry2,ayC,ayD,an3,an4,Cy,Dy);
        	
		int muy = ayA + ayB + ayC + ayD - 2.0 * (iy1 + iy2 + iy3 + iy4) - (oy1 + oy2 + oy3 + oy4);
        	for (int uy = 0; uy <= floor(muy / (double)2.0); ++uy) {

		double J3y = J3_func(uy,muy,nFACT,Py,Qy);
                
			//Indicies for Z calculations
                	//double J0z = J0_func(azA,azB,an1,an2);
                	for (int iz1 = 0; iz1 <= floor(azA / (double)2.0); ++iz1) {
                	for (int iz2 = 0; iz2 <= floor(azB / (double)2.0); ++iz2) {
                	for (int oz1 = 0; oz1 <= azA - 2.0 * iz1; ++oz1) {
                	for (int oz2 = 0; oz2 <= azB - 2.0 * iz2; ++oz2) {
                	for (int rz1 = 0; rz1 <= floor((oz1 + oz2) / (double)2.0); ++rz1) {

			double J1z = J1_func(iz1,iz2,oz2,oz2,rz1,azA,azB,azC,azD,an1,an2,an3,an4,Az,Bz);
                	
			for (int iz3 = 0; iz3 <= floor(azC / (double)2.0); ++iz3) {
                	for (int iz4 = 0; iz4 <= floor(azD / (double)2.0); ++iz4) {
                	for (int oz3 = 0; oz3 <= azC - 2.0 * iz3; ++oz3) {
                	for (int oz4 = 0; oz4 <= azD - 2.0 * iz4; ++oz4) {
                	for (int rz2 = 0; rz2 <= floor((oz3 + oz4) / (double)2.0); ++rz2) {

			double J2z = J2_func(iz3,iz4,oz3,oz4,rz2,azC,azD,an3,an4,Cz,Dz);
                	
			int muz = azA + azB + azC + azD - 2.0 * (iz1 + iz2 + iz3 + iz4) - (oz1 + oz2 + oz3 + oz4);
                	for (int uz = 0; uz <= floor(muz / (double)2.0); ++uz) {

			double J3z = J3_func(uz,muz,nFACT,Pz,Qz);
			double nu = mux + muy + muz - (ux + uy + uz);

				//Putting it all together
				double val = J3x * J2x * J1x * J0x * J3y * J2y * J1y * J0y * J3z * J2z * J1z * J0z * 2.0 * Fnu(nu,PmQ * nFACT);
				if (val > 0 && val > tol)
					{pos += val;}
				else if (val < 0 && abs(val) > tol)
					{neg += abs(val);}
				cout.setf( std::ios::fixed, std::ios::floatfield );
				cout << "   J3x: " << J3x << " J2x: " << J2x << " J1x: " << J1x << " J0x: " << J0x << " J3y: " << J3y << " J2y: " << J2y << " J1y: " << J1y << " J0y: " << J0y << " J3z: " << J3z << " J2z: " << J2z << " J1z: " << J1z << " J0z: " << J0z << " Fnu: " << Fnu(nu,PmQ * nFACT) << "\n";
        		}}}}}}}}}}}
        	}}}}}}}}}}}
        }}}}}}}}}}}

	rtnval = pos - neg;
	cout << "   RTNVAL: " << rtnval << "\n";
        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                Calculate the Nuclear Attraction Integrals
//      Function that calculates the kinetic contribution to the fock
//	matrix.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
Math functions nCk(), fact(), d_fact()
*/

double calculate_electron_repulsion(int basisA,int basisB, int basisC,int basisD,MemHandler *data_mem,scf_data *scf,STOnG_handler *STO,int assu,int test)
{
        double rtnval = 0;
	int N = STO->N;
	
	//Orbital Index
	int basis_A = data_mem->orb_idx[basisA]; //orbital type 1,2,3,4,5 ...
	int basis_B = data_mem->orb_idx[basisB];
	int basis_C = data_mem->orb_idx[basisC];
	int basis_D = data_mem->orb_idx[basisD];

	//Atmoic Index
	int atomic_data_A = data_mem->atom_idx[basisA]; //Gives index of the atom in question
	int atomic_data_B = data_mem->atom_idx[basisB];
	int atomic_data_C = data_mem->atom_idx[basisC];
	int atomic_data_D = data_mem->atom_idx[basisD];

	//Atomic Number
	int Anum = data_mem->atom_data[atomic_data_A].atomic_num;
	int Bnum = data_mem->atom_data[atomic_data_B].atomic_num;
	int Cnum = data_mem->atom_data[atomic_data_C].atomic_num;
	int Dnum = data_mem->atom_data[atomic_data_D].atomic_num;

	//A Position	
	double Ax = data_mem->atom_data[atomic_data_A].pos_xyz[0];
	double Ay = data_mem->atom_data[atomic_data_A].pos_xyz[1];
	double Az = data_mem->atom_data[atomic_data_A].pos_xyz[2];

	//B Position
	double Bx = data_mem->atom_data[atomic_data_B].pos_xyz[0];
	double By = data_mem->atom_data[atomic_data_B].pos_xyz[1];
	double Bz = data_mem->atom_data[atomic_data_B].pos_xyz[2];
	
        //C Position
        double Cx = data_mem->atom_data[atomic_data_C].pos_xyz[0];
        double Cy = data_mem->atom_data[atomic_data_C].pos_xyz[1];
        double Cz = data_mem->atom_data[atomic_data_C].pos_xyz[2];

        //D Position
        double Dx = data_mem->atom_data[atomic_data_D].pos_xyz[0];
        double Dy = data_mem->atom_data[atomic_data_D].pos_xyz[1];
        double Dz = data_mem->atom_data[atomic_data_D].pos_xyz[2];

	//|A-B|^2
	double AmB = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);

        //|C-D|^2
        double CmD = (Cx - Dx) * (Cx - Dx) + (Cy - Dy) * (Cy - Dy) + (Cz - Dz) * (Cz - Dz);
	
	//Used for indexing the constants array and obtaining the correct constants
	int basis_idx_A = STO[Anum].rtn_bas_idx(basis_A);
	int basis_idx_B = STO[Bnum].rtn_bas_idx(basis_B);
	int basis_idx_C = STO[Bnum].rtn_bas_idx(basis_C);
	int basis_idx_D = STO[Bnum].rtn_bas_idx(basis_D);

	//The following sets the angular numbers for the given orbital
	int axA,ayA,azA;
	int axB,ayB,azB;
	int axC,ayC,azC;
	int axD,ayD,azD;

	switch (basis_A)
	{
		case 1:{axA=0;ayA=0;azA=0;break;}
		case 2:{axA=0;ayA=0;azA=0;break;}
		case 3:{axA=1;ayA=0;azA=0;break;}
		case 4:{axA=0;ayA=1;azA=0;break;}
		case 5:{axA=0;ayA=0;azA=1;break;}
	} 

        switch (basis_B)
        {
                case 1:{axB=0;ayB=0;azB=0;break;}
                case 2:{axB=0;ayB=0;azB=0;break;}
                case 3:{axB=1;ayB=0;azB=0;break;}
                case 4:{axB=0;ayB=1;azB=0;break;}
                case 5:{axB=0;ayB=0;azB=1;break;}
        }

        switch (basis_C)
        {
                case 1:{axC=0;ayC=0;azC=0;break;}
                case 2:{axC=0;ayC=0;azC=0;break;}
                case 3:{axC=1;ayC=0;azC=0;break;}
                case 4:{axC=0;ayC=1;azC=0;break;}
                case 5:{axC=0;ayC=0;azC=1;break;}
        }

        switch (basis_D)
        {
                case 1:{axD=0;ayD=0;azD=0;break;}
                case 2:{axD=0;ayD=0;azD=0;break;}
                case 3:{axD=1;ayD=0;azD=0;break;}
                case 4:{axD=0;ayD=1;azD=0;break;}
                case 5:{axD=0;ayD=0;azD=1;break;}
        }	

	//This statement decides if the basis functions are orthogonal P functions, if they are set rtnval to 0.
	if ((basis_A >= 3 && basis_B >= 3 && basis_A != basis_B) || (basis_C >= 3 && basis_D >= 3 && basis_C != basis_D) && assu == 0)
	{
		rtnval = 0.00;
	//The following else statement calculates the rest of the integrals which are not orthogonal P functions.
	} else {
	for (int n1 = 0; n1 < N ;++n1) {
	for (int n2 = 0; n2 < N ;++n2) {
	for (int n3 = 0; n3 < N ;++n3) {	
	for (int n4 = 0; n4 < N ;++n4) {

			//Load primitive gaussian exponent constant	
			double an1 = STO[Anum - 1].a[n1 + (N * basis_idx_A - N)];
			double an2 = STO[Bnum - 1].a[n2 + (N * basis_idx_B - N)];
			double an3 = STO[Cnum - 1].a[n3 + (N * basis_idx_C - N)];
			double an4 = STO[Dnum - 1].a[n4 + (N * basis_idx_D - N)];

			//Load primitive gaussian contraction coefficient			
			double dn1 = STO[Anum - 1].d[n1 + (N * basis_idx_A - N)];
			double dn2 = STO[Bnum - 1].d[n2 + (N * basis_idx_B - N)];
			double dn3 = STO[Cnum - 1].d[n3 + (N * basis_idx_C - N)];
			double dn4 = STO[Dnum - 1].d[n4 + (N * basis_idx_D - N)];

			//Define GammaP and GammaQ
			double gammaP = an1 + an2;
			double gammaQ = an3 + an4;

			//Define constant factor
			double NRI_fact = pow(M_PI,5 / (double)2) / (double)(gammaP * gammaQ * sqrt(gammaP + gammaQ)); 
			
			cout << " n1:  " << n1 << " n2: " << n2 << " n3: " << n3 << " n4: " << n4 << "\n";

			double Ifact = calc_sum_J(test,axA,ayA,azA,axB,ayB,azB,axC,ayC,azC,axD,ayD,azD,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,an1,an2,an3,an4);

			rtnval += norm(axA,ayA,azA,an1) * norm(axB,ayB,azB,an2) * norm(axC,ayC,azC,an3) * norm(axD,ayD,azD,an4) * dn1 * dn2 * dn3 * dn4 * NRI_fact * Ifact * E_factor(an1,an2,AmB) * E_factor(an3,an4,CmD);
			//cout << "NORMA: " << norm(axA,ayA,azA,an1) << " NORMB: " << norm(axB,ayB,azB,an2) << " NORMC: " << norm(axC,ayC,azC,an3) << " NORMD: " << norm(axD,ayD,azD,an4) <<  " dn1: " << dn1 << " dn2: " << dn2 << " dn3: " << dn3 << " dn4: " << dn4 << " NRI_fact: " << NRI_fact << " IFACT: " << Ifact << " E_factor1: " << E_factor(an1,an2,AmB) << " E_factor2: " << E_factor(an3,an4,CmD) << "\n";
	}}}}
	}

       	return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                   Electron Repulsion Main Function
//      Create initial index and calculate the integrals for all
//	degenerate electron repulsion integrals. 
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/

extern void electron_repulsion_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO)
{
	//Set some variables
	timer int_timer;
        optfile->ofile << "\nBeginning Electron Repulsion Calculations...\n";
	int_timer.set_timer();

	//**************************************
	//Start Nuclear Attraction Computations	
	//**************************************

        int count = 0;
        scf->num_deg = 0;
	
	int tot1 = scf->orbitals * (scf->orbitals + 1) / 2;
	int total_calcs = tot1 * (tot1 + 1) / 2;
	int check = 0;

        int strider;
        for (int i = 0; i < scf->orbitals; ++i) {
        for (int j = i; j < scf->orbitals; ++j) {
        for (int k = i; k < scf->orbitals; ++k) {
        for (int l = k; l < scf->orbitals; ++l) {

                if ((long int)scf->FEI.produce_unique_ident(i,j,k,l) == (long int)scf->FEI.produce_ident(i,j,k,l))
                {
			//int testval = scf->FEI.produce_ident(i,j,k,l);
			//ofstream myfile;
                        //myfile.open ("edata.bin", ios::app);
			//myfile << "TEST: " << i << "-" << j << "|" << k << "-" << l << "\n";
                        //myfile.close();
			//scf->ERI_idx[scf->num_deg] = uint_val;

        		double ERIij;
			//ERIij = calculate_electron_repulsion(i,j,k,l,data_mem,scf,STO,1,0);
			ERIij = 0.0;
			//scf->ERI[scf->num_deg] = ERIij;
			
			scf->FEI.set_val(ERIij,i,j,k,l);
			//if (abs(ERIij - ERIij2) > 1.0E-10)
			//{
				//cout << i << "," << j << "|" << k << "," << l << " = " << ERIij << " = " << ERIij2 << "\n";
			//}
			//cout << scf->FEI.get_val(i,j,k,l) << "\n";
			double perc = (scf->num_deg / (double)total_calcs) * 100;
			
			if (perc >= check)
			{
				//cout << perc << "\%\n";
				check += 10;
			}

                        ++scf->num_deg;
                }
                ++count;
        }}}}

	cout << "|---------------------------------FIRST: \n\n";
	double Ifact1 = calculate_electron_repulsion(2,5,5,5,data_mem,scf,STO,1,0);
	cout << "|---------------------------------SECOND: \n\n";
	double Ifact2 = calculate_electron_repulsion(5,2,5,5,data_mem,scf,STO,1,0);
	cout.setf( std::ios::fixed, std::ios::floatfield );
	cout << "IFACT1: " << Ifact1 << " IFACT2: " << "" << Ifact2  << "\n";
	/*cout.setf( std::ios::fixed, std::ios::floatfield );
        for (int n = 0;n < scf->num_deg;++n)
              {cout << scf->ERI_idx[n] << " = " << scf->ERI[n] << " | ";if (n % 6 == 0) {cout << "\n";}}
        */
	//cout << "\n"; 
        cout << " count: " << count << "\n";
        cout << "FINISHED ARRAY: CALC_SIZE(" << scf->FEI.num_val << ")\n";
	
	//************************************
	//End Nuclear Attraction Computations
	//************************************

	//Some ending functions
	int_timer.end_timer();
	string message = "Electron Repulsion Calculationi Clock Time: ";
        int_timer.print_clock_time(message,optfile->ofile);
        optfile->ofile << "\n";
}

