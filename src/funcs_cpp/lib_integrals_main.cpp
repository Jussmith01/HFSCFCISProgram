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
#include "../scf_main/scf_classes.h"
#include "func_classes.h"
#include "lib_integrals.h"

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                        Integral Calculator Function
//      Calculates all integrals required for the SCF iterations.
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/15/2014  ----
                   ----Modified By:               ----
*/

extern void calculate_integrals(MemHandler *data_mem,scf_data *scf,dataOutput *optfile)
{
	optfile->ofile << "\n|------------------------------------------------------|\n";
	optfile->ofile << "|--------------Begin Integral Calculations-------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n";
	//********************************************//
        //          Load Basis Set Constants          //
        //********************************************//
	STOnG_handler *STO;
	STO = new STOnG_handler [10]; //Data will be shifted left by one i.e. Atomic number 1 = 0...
	optfile->ofile << "Reading Gaussian Type Orbital Constants...\n";
	cout << "Reading Gaussian Type Orbital Constants...\n";
	for (int i = 0; i < 10; ++i)
		{STO[i].read_input(data_mem->ipt_parms.nCGO,i+1,optfile);}
	
        //********************************************//
        //        Calculate Overlap Integrals         //
        //********************************************//
	overlap_main(data_mem,scf,optfile,STO);//Calculates the overlap matrix elements and stores them.
	optfile->Matrix_printer("Overlap Matrix",scf->S,scf->orbitals);

	//********************************************//
        //        Calculate Kinetic Integrals         //
        //********************************************//
        kinetic_main(data_mem,scf,optfile,STO);//Calculates the kinetic elements of F and stores them.
	optfile->Matrix_printer("Kinetix Matrix",scf->Hcore,scf->orbitals);

	//********************************************//
        //  Calculate Transition Density Integrals    //
        //********************************************//
	if (scf->td_flag==1)
	{
        	TD_main(data_mem,scf,optfile,STO);//Calculates the kinetic elements of F and stores them.
        	optfile->Matrix_printer("TDx Matrix",scf->TDx,scf->orbitals);
        	optfile->Matrix_printer("TDy Matrix",scf->TDy,scf->orbitals);
        	optfile->Matrix_printer("TDz Matrix",scf->TDz,scf->orbitals);
	}
 
	//********************************************//
        //     Calculate Nuc Attraction Integrals     //
        //********************************************//
	Null_Set(scf->WM,scf->orbitals);
        for (int ac = 0; ac < data_mem->ipt_parms.num_atoms; ++ac)
        {
        	nuclear_main(ac,data_mem,scf,optfile,STO);//Calculates the nuclear attraction elements of F and adds them.
		string title = "Nuclear Matrix ATOM(" + int_to_str(ac) + ")";
		optfile->Matrix_printer(title,scf->WM,scf->orbitals);
	}

        //********************************************//
        //           	   PRINT Hcore                //
        //********************************************//
	optfile->Matrix_printer("Hcore Matrix",scf->Hcore,scf->orbitals);

        //********************************************//
        //        Calculate Exchange Integrals        //
        //********************************************//
       	electron_repulsion_main(data_mem,scf,optfile,STO);//Calculates the electron repulsion integrals and stores them.

        for (int i = 0; i < 10; ++i)
                {STO[i].mem_free();}

        delete [] STO;
        STO = NULL;

	optfile->ofile << "\n|------------------------------------------------------|\n";
	optfile->ofile << "|---------------End Integral Calculations--------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n";
}
