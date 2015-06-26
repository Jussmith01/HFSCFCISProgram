#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "../utils_cpp/lib_includes.h"
#include "../classes/classes.h"
#include "scf_classes.h"

using namespace std;

// ********************************************************************* //
// ******************memhandler CLASS MEMBER FUNCTIONS****************** //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Allocates Necessary Memory ----
		   ----for SCF data.		  ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/14/2014  ----
                   ----Modified By:               ----
This program allocates all necessary memory for the self consistant field data.
*/
void scf_data::mem_alloc (dataOutput* optfile)
{
	//Allocate Overlap Matrix
	S = new double [orbitals * orbitals];

        //Allocate Overlap Matrix
        Sinv = new double [orbitals * orbitals];

        //Allocate Hcore Matrix
        Hcore = new double [orbitals * orbitals];
	
	//Allocate Fock Matrix
	F = new double [orbitals * orbitals];

        //Allocate C constants Matrix
        C = new double [orbitals * orbitals];

        //Allocate E Eigen Value Matrix
        E = new double [orbitals];
        
	//Allocate Density Matrix
        Rho = new double [orbitals * orbitals];
	for (int i = 0; i < orbitals * orbitals; ++i) {Rho[i] = 0;}//Inititialize Rho will NULL matrix

	if (td_flag==1)
	{
		TDx = new double [orbitals * orbitals];
		TDy = new double [orbitals * orbitals];
		TDz = new double [orbitals * orbitals];
	}

	//Allocate Workin  Matrix
	WM = new double [orbitals * orbitals];

	//Set Number of Degenerate Electron Repulsion integrals and the index
	num_deg = FEI.alloc_val_table(orbitals);


        int fset = orbitals * (orbitals + 1) / (int) (2);
        int degenerate = fset * (fset + 1) / (int) (2);

	memory_req = (((2 * orbitals) + (13 * orbitals * orbitals))* sizeof(double));//Num Matricies * orbitals^2 * numbytes
	optfile->ofile << "Total Number of Basis Functions: " << orbitals << "\n";

	//Set Unit Size for Memory Required output
	string unit = "B";
	double div = 1;
	if (memory_req >= 1024 && memory_req < 1024 * 1024)
		{div = 1024;unit = "KB";} 
	else if (memory_req >= 1024 * 1024 && memory_req < 1024 * 1024 *1024)
                {div = 1024 * 1024;unit = "MB";}
        else if (memory_req >= 1024 * 1024 * 1024)
                {div = 1024 * 1024 * 1024;unit = "GB";}

	optfile->ofile << "Memory Required for Computation: " << memory_req / (double)div << unit <<"\n";

}

/*____________________________________________________________________________
                   ----Frees All Memory Allocations ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/14/2014    ----
                   ----Modified By:                 ----
Frees any memory allocated by this class.
*/
void scf_data::mem_free ()
{
	//Free Overlap Matrix
	delete [] S;
        S = NULL;

        //Free Overlap Matrix
        delete [] Sinv;
        Sinv = NULL;

        //Free Fock Matrix
        delete [] Hcore;
        Hcore = NULL;

        //Free Fock Matrix
        delete [] C;
        C = NULL;

        //Free Fock Matrix
        delete [] E;
        E = NULL;

        //Free Fock Matrix
        delete [] Rho;
        Rho = NULL;

        //Free Fock Matrix
        delete [] F;
        F = NULL;

        //Free Working Matrix
        delete [] WM;
        WM = NULL;

        //Free TD Matricies
	if (td_flag==1)
	{
        	delete [] TDx;TDx = NULL;
        	delete [] TDy;TDy = NULL;
        	delete [] TDz;TDz = NULL;
	}

	//Free two electron integral values
	FEI.free_val_table();
}

