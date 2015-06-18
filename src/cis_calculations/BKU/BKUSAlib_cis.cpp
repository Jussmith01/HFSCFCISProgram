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
#include <vector>

//#include "../utils_cpp/lib_includes.h"
//#include "../classes/classes.h"
//#include "../scf_main/scf_classes.h"
#include "lib_cis.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//              Produce Pseudo Density Matrix Given Coefficents
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/


//________________________________________________________________________//
//      *************************************************************     //
//                    CIS MAIN CALC SINGLET SPIN ADAPTED
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 3/26/2015  ----
                   ----Modified By:               ----
*/

extern void cis_energy_calc(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	if (data_mem->ipt_parms.cis_flag==1)
	{
        	optfile->ofile << "|------------------------------------------------------|\n";
		optfile->ofile << "|-----------------Begin CIS Calculation----------------|\n";
        	optfile->ofile << "|------------------------------------------------------|\n";

		int N = scf->orbitals;
		optfile->Matrix_printer("SCF Fock Matrix: ",scf->F,N);

		optfile->ofile << "\nBuilding CIS Singlet Spin Adapted Hamitonian Matrix...\n";

		//Note, Excluding the ground state since it does not couple with the single excited states.
	        int NO = data_mem->num_elec / 2;
		int NV = N-NO;//Number of virtual orbitals

	        //********************************************//
	        //   Find Minimum Energy Orbital Occupancy    //
	        //********************************************//
	        int *min = new int [NO];
		find_mins(min,scf->E,N,NO);


                //********************************************//
                //         Determine Occ. Orb. Index          //
                //********************************************//
		optfile->ofile << "Sum of Occupied Eigenvectors: \n";
		vector<float> OccVecSum;
		for (int i=0;i<N;++i)
		{
			OccVecSum.push_back(0.0f);
			for (int j=0;j<NO;++j)
                	{
				OccVecSum[i] += scf->C[min[j] + i*N];
			}
			OccVecSum[i] = abs(OccVecSum[i]*2.0f);
			optfile->ofile << " (" << i << "): " << OccVecSum[i] << endl;
		}

		vector<int> Occupied;
		vector<int> Unoccupied;

		Occupied.push_back(0);
	//	Occupied.push_back(1);
		//Occupied.push_back(2);
		//Occupied.push_back(3);
		//Occupied.push_back(4);
		//Occupied.push_back(6);

                Unoccupied.push_back(0);
          //      Unoccupied.push_back(0);
                //Unoccupied.push_back(5);

                //********************************************//
                //        Build CIS Hamiltonian Matrix        //
                //********************************************//
		//cout << "MinSB Size: " << (int)minsb.size() << " MaxSB Size: " << (int)maxsb.size() << endl;
		double *Hcis = new double [2*NO * 2*NO];
		int Nuocc = (int)Unoccupied.size();
		//double Escf = data_mem->Eelec;
		for (int i=0;i<NO;++i)
                for (int a=0;a<Nuocc;++a)
                {
		   	for (int j=0;j<NO;++j)
		   	for (int b=0;b<Nuocc;++b)
		   	{
				int iidx = Occupied[i];
				int aidx = Unoccupied[a];
				int jidx = Occupied[j];
				int bidx = Unoccupied[b];

				cout << "i: " << iidx << " a: " << aidx << " j: " << jidx << " b: " << bidx; 

				double val = 0;
				//val += (((iidx == jidx) && (aidx == bidx)) ? (scf->E[aidx] - scf->E[iidx]) : 0);
				val += ((iidx == jidx) ? scf->F[aidx + bidx * N] : 0);
				val -= ((aidx == bidx) ? scf->F[iidx + jidx * N] : 0);
				val += 2.0 * scf->FEI.get_val(aidx,jidx,iidx,bidx) - scf->FEI.get_val(aidx,jidx,bidx,iidx);

                           	Hcis[(j*Nuocc+b) + (i*Nuocc+a) * NO*2] = val;

				cout << " " << (j*Nuocc+b) + (i*Nuocc+a) * NO*2 << " \n";
		   	}
			//cout << endl;
                }

		optfile->Matrix_printer("CIS Hamiltonian",Hcis,2*NO);
		
                //********************************************//
                //    Diagonalize CIS Hamiltonian Matrix      //
                //********************************************//
		cout << "TEST6\n";
		double *Evec = new double [2*NO * 2*NO];
		double *Eval = new double [2*NO];

                //Find eigenvectors and eigenvalues
                optfile->ofile << "Running Jacobi Algorithm... \n";

                double *convergence_data;
                convergence_data = new double [3];

		cout << "TEST7\n";
                convergence_data[0] = data_mem->ipt_parms.j_conv;
                convergence_data[1] = 32.0 * NO * NO;
                Jacobi_Algorithm(Evec,Eval,Hcis,2*NO,convergence_data);//Diagonlize H CIS and get eigenvectors and values
                optfile->ofile << "Jacobi Iterations: " << int_to_str(convergence_data[2]) <<" of " << int_to_str(convergence_data[1]) <<  "\n";
                delete [] convergence_data; convergence_data = NULL;

		optfile->Matrix_printer("CIS Eigenvectors",Evec,2*NO);
                optfile->Array_printer("CIS Eigenvalues (Hartree)",Eval,2*NO,4);

		//Convert to eV
		for (int i=0;i<2*NO;++i) {Eval[i]=Eval[i]*27.2107f;}

                optfile->Array_printer("CIS Eigenvalues (eV)",Eval,2*NO,4);
                //optfile->ofile << "Producing Density... \n";


		int *coeffidx = new int [2*NO];
		double *Rho = new double [2*NO * 2*NO];

		for (int i=0;i<2*NO;++i) {coeffidx[i] = i;}

		Density_Matrix(Rho,Evec,coeffidx,2*NO,NO);
                //optfile->Matrix_printer("CIS Density Matrix",Rho,2*NO);
		
		//********************************************//
                //    		   Clean Up                   //
                //********************************************//
		cout << "TEST8\n";
		delete [] Evec; Evec = NULL;
		delete [] Eval; Eval = NULL;
		delete [] Hcis; Hcis = NULL;
		delete [] Rho; Rho = NULL;
		delete [] coeffidx; coeffidx = NULL;

		cout << "TEST9\n";
                optfile->ofile << "|------------------------------------------------------|\n";
                optfile->ofile << "|------------------End CIS Calculation-----------------|\n";
                optfile->ofile << "|------------------------------------------------------|\n";
	} else {
        	optfile->ofile << "CIS Calculation Disabled... moving on.\n\n";
	}
}
