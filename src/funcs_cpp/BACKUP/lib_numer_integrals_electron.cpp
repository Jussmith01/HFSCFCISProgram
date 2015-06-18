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
#include <random>
#include <omp.h>

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

double function(double x1,double y1,double z1,double x2,double y2,double z2)
{
	double rtnval;

	rtnval = x1 * x1 * y1 * y1 * z1 * z1 * x2 * x2 * y2 * y2 * z2 * z2;
	
	return rtnval;
}

double w(double dx,int k,int M)
{
	double rtnval;

	if (k == 0 || k == M)
		{rtnval = 0.333333333 * dx;}
	else 
	{
		switch (k % 2)
		{
			case 0: {rtnval = 1.333333333 * dx;break;}
			case 1: {rtnval = 0.666666666 * dx;break;}
		}
	}

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

//double calc_numer_elec(int axA,int ayA,int azA,int axB,int ayB,int azB,int axC,int ayC,int azC,int axD,int ayD,int azD,double Ax,double Ay,double Az,double Bx,double By,double Bz,double Cx,double Cy,double Cz,double Dx,double Dy,double Dz,double an1,double an2,double an3,double an4)
extern double calc_numer_elec()
{
	double tol = 1.0E-13;
        double I = 0;
	
	int M = 100;

	double x10 = 0.0;
	double x1f = 2.0;
	double y10 = 0.0;
	double y1f = 2.0;
	double z10 = 0.0;
	double z1f = 2.0;

        double x20 = 0.0;
        double x2f = 2.0;
        double y20 = 0.0;
        double y2f = 2.0;
        double z20 = 0.0;
        double z2f = 2.0;

	double d1x = (x1f - x10) / (double)M;
	double d1y = (y1f - y10) / (double)M;
	double d1z = (z1f - z10) / (double)M;

        double d2x = (x2f - x20) / (double)M;
        double d2y = (y2f - y20) / (double)M;
        double d2z = (z2f - z20) / (double)M;

	#pragma omp parallel for shared(x10,x20,y10,y20,z10,z20,d1x,d1y,d1z,d2x,d2y,d2z,M,I)
	for (int k1 = 0; k1 <= M; ++k1)
	{
		double K1sum = 0;
		double x1val = x10 + k1 * d1x;
		double w1 = w(d1x,k1,M);

		for (int k2 = 0; k2 <= M; ++k2){
		double y1val = y10 + k2 * d1y;
		double w2 = w(d1y,k2,M);
		double w2tot = w1 * w2;		

		for (int k3 = 0; k3 <= M; ++k3){
		double z1val = z10 + k3 * d1z;
		double w3 = w(d1z,k3,M);
		double w3tot = w3 * w2tot;		

		for (int k4 = 0; k4 <= M; ++k4){
		double x2val = x20 + k4 * d2x;
		double w4 = w(d2x,k4,M);
		double w4tot = w4 * w3tot;		
	
		for (int k5 = 0; k5 <= M; ++k5){
		double y2val = y20 + k5 * d2y;
		double w5 = w(d2y,k5,M);
		double w5tot = w5 * w4tot;		
		
		for (int k6 = 0; k6 <= M; ++k6){
		double z2val = z20 + k6 * d2z;
		double w6 = w(d2z,k6,M);
		double w6tot = w6 * w5tot;		

			double fxyz = function(x1val,y1val,z1val,x2val,y2val,z2val);

			K1sum += w6tot * fxyz;
		}}}}}

		#pragma omp critical
			{I += K1sum;cout << "K1 = " << k1 << "\n";}
	}

        return I;
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

/*double calculate_electron_repulsion(int basisA,int basisB, int basisC,int basisD,MemHandler *data_mem,scf_data *scf,STOnG_handler *STO,int assu,int test)
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
	//if ((basis_A >= 3 && basis_B >= 3 && basis_A != basis_B) || (basis_C >= 3 && basis_D >= 3 && basis_C != basis_D) && assu == 0)
	//{
	//	rtnval = 0.00;
	//The following else statement calculates the rest of the integrals which are not orthogonal P functions.
	//} else {
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

			//cout << " n1:  " << n1 << " n2: " << n2 << " n3: " << n3 << " n4: " << n4 << "\n";

			double NumerIfact = calc_numer_elec(axA,ayA,azA,axB,ayB,azB,axC,ayC,azC,axD,ayD,azD,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,an1,an2,an3,an4);

			rtnval += norm(axA,ayA,azA,an1) * norm(axB,ayB,azB,an2) * norm(axC,ayC,azC,an3) * norm(axD,ayD,azD,an4) * dn1 * dn2 * dn3 * dn4 * NumerIfact;
			//cout << "NORMA: " << norm(axA,ayA,azA,an1) << " NORMB: " << norm(axB,ayB,azB,an2) << " NORMC: " << norm(axC,ayC,azC,an3) << " NORMD: " << norm(axD,ayD,azD,an4) <<  " dn1: " << dn1 << " dn2: " << dn2 << " dn3: " << dn3 << " dn4: " << dn4 << " NRI_fact: " << NRI_fact << " IFACT: " << Ifact << " E_factor1: " << E_factor(an1,an2,AmB) << " E_factor2: " << E_factor(an3,an4,CmD) << "\n";
	}}}}
	//}

       	return rtnval;
}*/

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

/*extern void electron_repulsion_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO)
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
			ERIij = calculate_electron_repulsion(i,j,k,l,data_mem,scf,STO,1,0);
			//ERIij = 0.0;
			//scf->ERI[scf->num_deg] = ERIij;
			
			scf->FEI.set_val(ERIij,i,j,k,l);
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

	cout.setf( std::ios::fixed, std::ios::floatfield );
        for (int n = 0;n < scf->FEI.num_val;++n)
             {cout << scf->FEI.Val_Table[n].key << " = " << scf->FEI.Val_Table[n].value  << " | ";if (n % 6 == 0) {cout << "\n";}}
        
	cout << "FINISHED ARRAY: CALC_SIZE(" << scf->FEI.num_val << ")\n";
	
	//************************************
	//End Nuclear Attraction Computations
	//************************************

	//Some ending functions
	int_timer.end_timer();
	string message = "Electron Repulsion Calculationi Clock Time: ";
        int_timer.print_clock_time(message,optfile->ofile);
        optfile->ofile << "\n";
}*/

