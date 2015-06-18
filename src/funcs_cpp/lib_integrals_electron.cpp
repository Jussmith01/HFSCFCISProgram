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

//      *************************************************************     //
//                              Calc_J_main Terms 
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/

double J_main(int i1,int i2,int r1,int r2,int u,double gamma1,double gamma2,double delta,double PmQxyz)
{
        double rtnval;

	double mu = i1 + i2 - 2.0 * (r1 + r2);

        double fact1 = ((fact(i1) * pow(4.0 * gamma1,r1 - i1)) / (double)(fact(r1) * fact(i1 - 2.0 * r1)));
        double fact2 = pow(-1,i2) * ((fact(i2) * pow(4.0 * gamma2,r2 - i2)) / (double)(fact(r2) * fact(i2 - 2.0 * r2)));
        //double fact1 = ((fact(i1) * fact(i2)) / (double)(pow(4 * gamma1, i1) * pow(4 * gamma2, i2) * pow(delta, i1 + i2)));
        //double fact2 = pow(-1,i2) * ((pow(4 * gamma1, r1) * pow(4 * gamma2, r2) * pow(delta, r1 + r2)) / (double)(fact(r1) * fact(r2) * fact(i1 - 2 * r1) * fact(i2 - 2 * r2)));
	double fact3 = (pow(-1,u) * fact(mu) * pow(PmQxyz,mu - 2.0 * u)) / (double)(fact(u) * fact(mu - 2.0 * u) * pow(delta,mu - u));

	rtnval = fact1 * fact2 * fact3;

	//cout << "FACT1: " <<  fact1 << " FACT2: " << fact2 << " FACT3: " << fact3 << "\n";
	//cout << "num: " <<  fact(i1) * fact(i2) << " Denom: " << pow(4 * gamma1, i1) * pow(4 * gamma2, i2) * pow(delta,i1 + i2) << "\n";

        return rtnval;
}

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

	double gamma1 = an1 + an2;
	double gamma2 = an3 + an4;

	double delta = (1.0/(double)(4.0 * gamma1)) + (1.0/(double)(4.0 * gamma2));

        double Px = (an1 * Ax + an2 * Bx) / (double)(an1 + an2);
        double Py = (an1 * Ay + an2 * By) / (double)(an1 + an2);
        double Pz = (an1 * Az + an2 * Bz) / (double)(an1 + an2);
        
        double Qx = (an3 * Cx + an4 * Dx) / (double)(an3 + an4);
        double Qy = (an3 * Cy + an4 * Dy) / (double)(an3 + an4);
        double Qz = (an3 * Cz + an4 * Dz) / (double)(an3 + an4);
        
        double PmQx = Qx - Px;
        double PmQy = Qy - Py;
        double PmQz = Qz - Pz;
        double PmQ = PmQx * PmQx + PmQy * PmQy + PmQz * PmQz;
	//double PmQ = ((an1 * pow(pow(Ax,2.0) + pow(Ay,2.0) + pow(Az,2.0),0.5) - (an2 * pow(pow(Bx,2.0) + pow(By,2.0) + pow(Bz,2.0),0.5)) / (double)(an1 + an2))) - ((an3 * pow(pow(Cx,2.0) + pow(Cy,2.0) + pow(Cz,2.0),0.5) - (an4 * pow(pow(Dx,2.0) + pow(Dy,2.0) + pow(Dz,2.0),0.5)) / (double)(an3 + an4)));
        //cout << "CmP: " << CmP << " an: " << an << " bm:" << bm << "\n";
        //cout << "Cx: " << Cx << " Cy: " << Cy << " Cz: " << Cz << " Px: " << Px  << " Py: " << Py << " Pz: " << Pz << "\n";
        //cout << "axA: " << axA << " ayA: " << ayA << " azA: " << azA << " | " ;
        //cout << " axB: " << axB << " ayB: " << ayB << " azB: " << azB << " | " ;
        //cout << " axC: " << axC << " ayC: " << ayC << " azC: " << azC << " | " ;
        //cout << " axD: " << axD << " ayD: " << ayD << " azD: " << azD << " | "  << " PmQ: " << PmQ << "\n";


	double nFACT = ((an1 + an2) * (an3 + an4)) / (double)((an1 + an2) + (an3 + an4));

        //Indicies for X calculations
        for (int ix1 = 0; ix1 <= axA + axB; ++ix1) {
        for (int ix2 = 0; ix2 <= axC + axD; ++ix2) {
        for (int rx1 = 0; rx1 <= floor(ix1/(float)2.0); ++rx1) {
        for (int rx2 = 0; rx2 <= floor(ix2/(float)2.0); ++rx2) {
        for (int ux = 0; ux <= floor(((ix1 + ix2)/(float)2.0) - rx1 - rx2); ++ux) {

	double Jx = fi(ix1,axA,axB,Px - Ax,Px - Bx) * fi(ix2,axC,axD,Qx - Cx,Qx - Dx) * J_main(ix1,ix2,rx1,rx2,ux,gamma1,gamma2,delta,PmQx);

        	//Indicies for Y calculations
        	for (int iy1 = 0; iy1 <= ayA + ayB; ++iy1) {
        	for (int iy2 = 0; iy2 <= ayC + ayD; ++iy2) {
        	for (int ry1 = 0; ry1 <= floor(iy1/(float)2.0); ++ry1) {
        	for (int ry2 = 0; ry2 <= floor(iy2/(float)2.0); ++ry2) {
        	for (int uy = 0; uy <= floor(((iy1 + iy2)/(float)2.0) - ry1 - ry2); ++uy) {

        	double Jy = fi(iy1,ayA,ayB,Py - Ay,Py - By) * fi(iy2,ayC,ayD,Qy - Cy,Qy - Dy) * J_main(iy1,iy2,ry1,ry2,uy,gamma1,gamma2,delta,PmQy);

                	//Indicies for Y calculations
                	for (int iz1 = 0; iz1 <= azA + azB; ++iz1) {
                	for (int iz2 = 0; iz2 <= azC + azD; ++iz2) {
                	for (int rz1 = 0; rz1 <= floor(iz1/(float)2.0); ++rz1) {
                	for (int rz2 = 0; rz2 <= floor(iz2/(float)2.0); ++rz2) {
                	for (int uz = 0; uz <= floor(((iz1 + iz2)/(float)2.0) - rz1 - rz2); ++uz) {

                	double Jz = fi(iz1,azA,azB,Pz - Az,Pz - Bz) * fi(iz2,azC,azD,Qz - Cz,Qz - Dz) * J_main(iz1,iz2,rz1,rz2,uz,gamma1,gamma2,delta,PmQz);
			//cout << "   fi1: " << fi(iz1,azA,azB,Pz - Az,Pz - Bz) << " fi2: " << fi(iz2,azC,azD,Qz - Cz,Qz - Dz) << " J_main: " << J_main(iz1,iz2,rz1,rz2,uz,gamma1,gamma2,delta,PmQz) << "\n";
				
				double nu = (ix1 + ix2 + iy1 + iy2 + iz1 + iz2) - 2.0 * (rx1 + rx2 + ry1 + ry2 + rz1 + rz2) - (ux + uy + uz);
				double val = Jx * Jy * Jz * Fnu(nu, PmQ / (double)(4 * delta));
                                if (val > 0 && val > tol)
                                        {pos += val;}
                                else if (val < 0 && abs(val) > tol)
                                        {neg += abs(val);}
				//cout << "      Jx: " << Jx << " Jy: " << Jy << " Jz: " << Jz << " Fnu: " << Fnu(nu,PmQ * nFACT) << "\n";
                	}}}}}
        	}}}}}
        }}}}}

	rtnval = pos - neg;
	//cout << "   RTNVAL: " << rtnval << "\n";
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

double calculate_electron_repulsion(int basisA,int basisB, int basisC,int basisD,MemHandler *data_mem,scf_data *scf,STOnG_handler *STO,int assu,int test,dataOutput *optfile)
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
	/*if ((basis_A >= 3 && basis_B >= 3 && atomic_data_A == atomic_data_B && basis_A != basis_B) || (basis_C >= 3 && basis_D >= 3 && atomic_data_C == atomic_data_D && basis_C != basis_D))
	{
		rtnval = 0.00;
		//cout << "SKIPPING...\n";
	//The following else statement calculates the rest of the integrals which are not orthogonal P functions.
	} else {*/
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
			double NRI_fact = ((2.0 * M_PI * M_PI) / (double)(gammaP * gammaQ)) * sqrt((M_PI) / (double)(gammaP + gammaQ)); 
			
			//cout << " n1:  " << n1 << " n2: " << n2 << " n3: " << n3 << " n4: " << n4 << "\n";

			double Ifact = calc_sum_J(test,axA,ayA,azA,axB,ayB,azB,axC,ayC,azC,axD,ayD,azD,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,an1,an2,an3,an4);

			rtnval += norm(axA,ayA,azA,an1) * norm(axB,ayB,azB,an2) * norm(axC,ayC,azC,an3) * norm(axD,ayD,azD,an4) * dn1 * dn2 * dn3 * dn4 * NRI_fact * Ifact * E_factor(an1,an2,AmB) * E_factor(an3,an4,CmD);
			//cout << "NORMA: " << norm(axA,ayA,azA,an1) << " NORMB: " << norm(axB,ayB,azB,an2) << " NORMC: " << norm(axC,ayC,azC,an3) << " NORMD: " << norm(axD,ayD,azD,an4) <<  " dn1: " << dn1 << " dn2: " << dn2 << " dn3: " << dn3 << " dn4: " << dn4 << " NRI_fact: " << NRI_fact << " IFACT: " << Ifact << " E_factor1: " << E_factor(an1,an2,AmB) << " E_factor2: " << E_factor(an3,an4,CmD) << "\n";
	}}}}

        if (optfile->verbose == 1)
        	{optfile->ofile << "I: " << basisA + 1 << " J: " << basisB + 1 << " K: " << basisC + 1 << " L: " << basisD + 1 << " Value =  " << rtnval << "\n";}

	//}

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

        if (optfile->verbose == 1)
                {optfile->ofile << "\n|-------TWO ELECTRON INTEGRALS--------|\n";}

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
			ERIij = calculate_electron_repulsion(i,j,k,l,data_mem,scf,STO,1,0,optfile);
			//ERIij = 0.0;
			//scf->ERI[scf->num_deg] = ERIij;
			
			scf->FEI.set_val(ERIij,i,j,k,l);
			//cout << scf->FEI.get_val(i,j,k,l) << "\n";
			//double perc = (scf->num_deg / (double)total_calcs) * 100;
			
			//if (perc >= check)
			//{
			//	cout << perc << "\%\n";
			//	check += 10;
			//}

                        ++scf->num_deg;
                }
                ++count;
        }}}}

        if (optfile->verbose == 1)
                {optfile->ofile << "\nUnique Integrals: " << scf->num_deg << "\n\n";}

	/*int unmatched = 0;
        for (int i = 0; i < scf->orbitals; ++i) {
        for (int j = 0; j < scf->orbitals; ++j) {
        for (int k = 0; k < scf->orbitals; ++k) {
        for (int l = 0; l < scf->orbitals; ++l) {

                        double ERI = calculate_electron_repulsion(i,j,k,l,data_mem,scf,STO,1,0);
                        double ERI2 = scf->FEI.get_val(i,j,k,l);

                        if (abs(ERI - ERI2) > 1.0E-8)
                        {
                                cout << i << "," << j << "|" << k << "," << l << " = " << ERI << " /= " << ERI2 << " IDENT: " << scf->FEI.produce_unique_ident(i,j,k,l) << " - " << scf->FEI.produce_ident(i,j,k,l) << "\n";
                        	++unmatched;
			}
        }}}}
	*/
	/*cout << "|---------------------------------FIRST: \n\n";
	double Ifact1 = calculate_electron_repulsion(0,5,2,2,data_mem,scf,STO,1,0);
	cout << "|---------------------------------SECOND: \n\n";
	double Ifact2 = calculate_electron_repulsion(2,2,5,0,data_mem,scf,STO,1,0);
	cout.setf( std::ios::fixed, std::ios::floatfield );
	cout << "IFACT1: " << Ifact1 << " IFACT2: " << "" << Ifact2  << "\n";
	*/
	/*cout.setf( std::ios::fixed, std::ios::floatfield );
        for (int n = 0;n < scf->FEI.num_val;++n)
             {cout << scf->FEI.Val_Table[n].key << " = " << scf->FEI.Val_Table[n].value  << " | ";if (n % 6 == 0) {cout << "\n";}}
        */
	//cout << "\n";
	//int tot = pow(2,scf->orbitals + scf->orbitals);
        //cout << " count: " << count << " Unmatched: " << unmatched << " TotPossible: " << tot <<"\n";
        //cout << "FINISHED ARRAY: CALC_SIZE(" << scf->FEI.num_val << ")\n";
	
	//************************************
	//End Nuclear Attraction Computations
	//************************************

	//Some ending functions
	int_timer.end_timer();
	string message = "Electron Repulsion Calculationi Clock Time: ";
        int_timer.print_clock_time(message,optfile->ofile);
        optfile->ofile << "\n";
}

