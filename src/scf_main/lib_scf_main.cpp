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
#include "../pop_analysis_cpp/lib_pop_analysis.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                     Self Consistant Field Function
//      Carrys out the SCF iteration calculations and obtains energy.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/31/2014  ----
                   ----Modified By:               ----
*/

extern string scf_iteration_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	optfile->ofile << "\n|------------------------------------------------------|\n";
	optfile->ofile << "|----------------Begin SCF Calculations----------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n";
	
	string errstr = "";

	//********************************************//
        //              Calc S^(-1/2)                 //
        //********************************************//
	optfile->ofile << "Calculate S^(-1/2)... \n";
	int N = scf->orbitals;

	Inverse_sqrt(scf->Sinv,scf->S,N,optfile);//Produce S^(1/2)

	optfile->Matrix_printer("S^(-1/2)",scf->Sinv,N);

	int it = 0;
		
	//memcpy(scf->F,scf->Hcore,N * N * sizeof(double));
	
	double energy_old = 0;
	double scf_converg = 0;

	Null_Set(scf->Rho,N);
	Null_Set(scf->F,N);
	Null_Set(scf->WM,N);

	huckel_main(data_mem,scf,optfile);

	//Produce_Initial_Density(data_mem,scf,optfile);
	Produce_Electronic_Terms(data_mem,scf,optfile);

        double EnergyStart = Calc_Electronic_Energy(data_mem,scf,optfile);
        optfile->ofile << "\nCalculated Starting Electronic Energy: " << EnergyStart << "\n";

	optfile->Matrix_printer("Rho Initial",scf->Rho,N);
	optfile->ofile << "Trace of Rho: " << Trace(scf->Rho,N) << "\n";

	while (abs(scf_converg) > data_mem->ipt_parms.e_conv || it < 2)
	{
        	optfile->ofile << "\n|========================= Iteration " << it + 1 << " =========================|\n";
		
		timer scf_timer;
	        scf_timer.set_timer();

		Null_Set(scf->C,N);
		Null_Set(scf->WM,N);

		//Produce_Electronic_Terms(data_mem,scf,optfile);
		
		optfile->Matrix_printer("F",scf->F,N);

		Matrix_Mult(scf->WM,scf->F,0,scf->Sinv,0,N);//Calculate F = S^(-1/2) * F' * S(-1/2)
        	Matrix_Mult(scf->F,scf->Sinv,0,scf->WM,0,N);

		Set_EPS_d(scf->F,N);

		optfile->Matrix_printer("F'",scf->F,N);

		//Find eigen vectors and values
		optfile->ofile << "Running Jacobi Algorithm... \n";

		double *convergence_data;
		convergence_data = new double [3];

        	convergence_data[0] = data_mem->ipt_parms.j_conv;
	        convergence_data[1] = 8.0 * N * N;
		Jacobi_Algorithm(scf->WM,scf->E,scf->F,N,convergence_data);//Calculate E and C' from F'
		optfile->ofile << "Jacobi Iterations: " << int_to_str(convergence_data[2]) <<" of " << int_to_str(convergence_data[1]) <<  "\n";
		delete [] convergence_data; convergence_data = NULL;

		optfile->Matrix_printer("Eigen Vectors'",scf->WM,N);
        	Matrix_Mult(scf->C,scf->WM,0,scf->Sinv,0,N);//Calculate C = S^(-1/2) * C'

		Set_EPS_d(scf->C,N);
		optfile->Matrix_printer("Eigen Vectors",scf->C,N);
		optfile->Array_printer("Eigen Values",scf->E,N,4);

		Null_Set(scf->Rho,N);
		Produce_Density(data_mem,scf,optfile);
		Null_Set(scf->F,N);
		
		Produce_Electronic_Terms(data_mem,scf,optfile);
		double Energy = Calc_Electronic_Energy(data_mem,scf,optfile);

		scf_converg = Energy - energy_old;
		energy_old = Energy;

        	optfile->ofile << "\nCalculated Electronic Energy of this Iteration: " << Energy << "\n";
        	//cout << "STEP(" << it << "): " << Energy << "\n";
        	//cout << "\nIteration (" << it << ") Energy: " << Energy << " Convergence: " << abs(scf_converg) << "\n";

		optfile->Matrix_printer("Rho",scf->Rho,N);
		optfile->ofile << "Trace of Rho: " << Trace(scf->Rho,N) << "\n";
        	
		scf_timer.end_timer();
        	string message = "SCF Cycle Clock Time: ";
        	scf_timer.print_clock_time(message,optfile->ofile);

		if (it == 128)
		{
			optfile->ofile << "\n!!!!!! SCF FAILED TO CONVERGE !!!!!!\n";
			errstr = "!! SCF FAILED TO CONVERGE !!";
			break;
		}
		++it;
	}
	data_mem->Eelec = energy_old;
        //cout << "Iteration (" << it << ") Energy: " << energy_old << "\n";
        optfile->ofile << "\n|------------------------------------------------------|\n";
        optfile->ofile << "|-----------------End SCF Calculations-----------------|\n";
        optfile->ofile << "|------------------------------------------------------|\n\n";

	return errstr;
}
