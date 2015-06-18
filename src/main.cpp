//**C/C++ included libraries**
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>

//**Included CUDA libraries**

//**Included custom headers**
#include "utils_cpp/lib_includes.h"
#include "classes/classes.h"
#include "scf_main/scf_classes.h"
#include "funcs_cpp/func_classes.h"
#include "pop_analysis_cpp/lib_pop_analysis.h"
#include "cis_calculations/lib_cis.h"

using namespace std;

int main (int argc, char *argv[])
{
        //********************************************//
        //               Declarations                 //
        //********************************************//
	char *inputfile = argv[1];
	char *opfilename = argv[2];
	MemHandler data_mem;//Define Memory Handler
	dataOutput optfile;//Define Output Handler
	timer wall_timer;//Define Timer Handler

        //********************************************//
        //              Preamble Programs             //
        //********************************************//
	wall_timer.set_timer();//Set the program wall timer
	optfile.set_output(opfilename);//Set output file

        //********************************************//
        //              Read Inputs                   //
        //****************************************** **//
	read_input(inputfile,&data_mem,&optfile);//Read the input file
	data_mem.ConvertDistances();//Convert distances to AU
	optfile.initial_data_printer(data_mem);//Print input data to output file

        //********************************************//
        //        Program Memory Allocations          //
        //********************************************//
	scf_data scf;//Define SCF data class
	data_mem.ipt_parms.define_parms(&optfile);//Defines the parameters from the input
	data_mem.num_basis_funcs();//Set number of basis functions
	data_mem.set_orb_idx(&optfile);//Indexes the type of orbitals
	scf.orbitals = data_mem.orbs;// set orbitals in scf class
	scf.td_flag = data_mem.ipt_parms.td_flag;// Set td flag
	scf.mem_alloc(&optfile);//Allocate Matricies used during the scf

	//********************************************//
        //          Main Program Execution            //
        //********************************************//
	calculate_integrals(&data_mem,&scf,&optfile);//Calculate all integrals for SCF
	string errchk = scf_iteration_main(&data_mem,&scf,&optfile);//Carry out SCF iterations and obtain E
	calc_nuclear_repulsion(&data_mem,&optfile);//Calculate the nuclear repulsion term under the BO appoximation
	population_analysis(&data_mem,&scf,&optfile);//Carryout Density Population Analysis
	cout << " " << data_mem.Eelec + data_mem.Enuc << "   " << data_mem.Eelec << "    " << data_mem.Enuc << " " << errchk << "\n";
	cis_energy_calc(&data_mem,&scf,&optfile);//Carryout CIS Calculation

        //********************************************//
        //           Print Energy Results             //
        //********************************************//
	optfile.ofile << "HF Energies: \n";
	optfile.ofile << "Etot = ";
	optfile.Prec_Printer(data_mem.Eelec + data_mem.Enuc,10); 
	optfile.ofile << " Eelec = ";
	optfile.Prec_Printer(data_mem.Eelec,10); 
	optfile.ofile << " Enuc = "; 
	optfile.Prec_Printer(data_mem.Enuc,10); 
	optfile.ofile << " " << errchk <<  "\n\n";	


        //********************************************//
        //            Free Allocated Memory           //
        //********************************************//
	data_mem.mem_free(1);//Free memory from data_mem
	scf.mem_free();// Free SCF memory

        //********************************************//
        //              Post Run Programs             //
        //********************************************//
        wall_timer.end_timer();//End the program wall timer
	string message = "Computation Wall Time: ";
	string message2 = "Computation Clock Time: ";
        wall_timer.print_clock_time(message2,optfile.ofile);//Print the wall time to output
        wall_timer.print_time(message,optfile.ofile);//Print the wall time to output
        optfile.close_output(0);//Close the output file

	return 0;
}
