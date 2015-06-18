#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <iomanip>

//**Included custom headers**
#include "lib_includes.h"
#include "../classes/classes.h"

using namespace std;

// ********************************************************************* //
// ****************DATA OUTPUT CLASS MEMBER FUNCTIONS******************* //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Precision String Return    ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
*/
void dataOutput::Prec_Printer(double A,int N)
{
        ofile.setf( std::ios::fixed, std::ios::floatfield );
        double val = A;
        int size;

        if (abs(val) >= 10)
                {size = N-1;}
        else if (abs(val) >= 100)
                {size = N-2;}
        else if (abs(val) >= 1000)
                {size = N-3;}
        else if (abs(val) >= 10000)
                {size = N-4;}
        else {size = N;}

        if (val < 0)
        {
                size -= 1;
        }

        ofile << setprecision(size) << val;
};

/*____________________________________________________________________________
                   ----Print an Array             ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
Row major formatting.
*/
void dataOutput::Array_printer(string Title,double *A,int N,int numPcol)
{
        ofile.setf( std::ios::fixed, std::ios::floatfield );
        if (verbose <= 1)
	{
		ofile << "\n_______________" << Title << "______________\n";
        	for (int i = 0; i < N; ++i)
        	{
        		double val = A[i];
        	        int size;
	
        	        if (abs(val) >= 10)
        	                {size = 7;}
        	        else if (abs(val) >= 100)
        	                {size = 6;}
        	        else if (abs(val) >= 1000)
        	                {size = 5;}
        	        else if (abs(val) >= 10000)
        	                {size = 4;}
        	        else {size = 8;}

        	        if (val < 0)
        	        {
        	                size -= 1;
        	        }

        		ofile << "(" << i << ") " << val << "  ";
			if (i % numPcol == numPcol - 1)
				{ofile << "\n";}
        	}
		ofile << "\n\n";
	}
};

/*____________________________________________________________________________
                   ----Print an Array (Hartree->eV)----
                   ----Author: Justin Smith        ----
                   ----Date Modified: 4/22/2015    ----
                   ----Modified By:                ----
Row major formatting.
*/
void dataOutput::Array_printer_harTOeV(string Title,double *A,int N,int numPcol)
{
	double convo = 27.21138505;
        ofile.setf( std::ios::fixed, std::ios::floatfield );
        if (verbose <= 1)
        {
                ofile << "\n_______________" << Title << "______________\n";
                for (int i = 0; i < N; ++i)
                {
                        double val = A[i] * convo;
                        int size;

                        if (abs(val) >= 10)
                                {size = 7;}
                        else if (abs(val) >= 100)
                                {size = 6;}
                        else if (abs(val) >= 1000)
                                {size = 5;}
                        else if (abs(val) >= 10000)
                                {size = 4;}
                        else {size = 8;}

                        if (val < 0)
                        {
                                size -= 1;
                        }

                        ofile << "(" << i << ") " << val << "  ";
                        if (i % numPcol == numPcol - 1)
                                {ofile << "\n";}
                }
                ofile << "\n\n";
        }
};

/*____________________________________________________________________________
                   ----Print a Matrix 		  ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
Row major formatting.
*/
void dataOutput::Matrix_printer(string Title,double *A,int N)
{
	ofile.setf( std::ios::fixed, std::ios::floatfield );
        if (verbose == 1)
	{
		ofile << "\n_______________" << Title << "______________\n";
        	for (int i = 0; i < N; ++i)
        	{
        		ofile << " | ";
        	        for (int j = 0; j < N; ++j)
	                {
				double val = A[j + i * N];
				int size;
	
				if (abs(val) >= 10)
					{size = 7;}
				else if (abs(val) >= 100)
					{size = 6;}
				else if (abs(val) >= 1000)
					{size = 5;}
				else if (abs(val) >= 10000)
					{size = 4;}
				else {size = 8;}

				if (val < 0)
				{
					size -= 1;
				}
                	
				ofile << setprecision(size) << A[j + i * N] << "  ";
        	        }
        	        ofile << " |\n";
        	}
	}
	/*ofile << "mat = {";
	for (int i = 0; i < N; ++i)
        {
                ofile << "{";
                for (int j = 0; j < N; ++j)
                {
                        double val = A[j + i * N];
                        int size;

                        if (abs(val) >= 10)
                                {size = 7;}
                        else if (abs(val) >= 100)
                                {size = 6;}
                        else if (abs(val) >= 1000)
                                {size = 5;}
                        else if (abs(val) >= 10000)
                                {size = 4;}
                        else {size = 8;}

                        if (val < 0)
                        {
                                size -= 1;
                        }
			if (j != N - 1)
			{
                        	ofile << setprecision(size) << A[j + i * N] << ",";
                	} else {
                        	ofile << setprecision(size) << A[j + i * N];
			}
		}
		if (i != N - 1)
		{
                	ofile << "},";
		} else {
                	ofile << "}";
		}
        }
	ofile << "}\n";*/
};

/*____________________________________________________________________________
                   ----Create and Open Output File----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
This opens the output file by the name of 'filename'. If the output filename 
was not included as the second argument upon command execution then the
output will be sent to a file named DEFAULT_OUTPUT.opt.
*/
void dataOutput::set_output (char *filename)
{
	//Define output name string
	string OutputFile;

	//Check if name is not included as the argument, then set output name.
	if (filename == NULL)
	{
		OutputFile = "DEFAULT_OUTPUT.opt";
	} else {
		OutputFile = filename;	
	}

	//Open the output file by given filename.
	ofile.open(OutputFile.c_str());
};

/*____________________________________________________________________________
                   ----Close Output File          ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/21/2014   ----
                   ----Modified By:               ----
This saves the data and closes the output file. Do this before ending the 
program. When input int close_type = 0 closes normally, when 1, abnormal
termination.
*/
void dataOutput::close_output (int close_type)
{
	switch(close_type)
	{
		case 0: {ofile << "Normal Termination\n";break;}
		case 1: {ofile << "Abnormal Termination\n";break;}
	}

        ofile.close();
};

/*____________________________________________________________________________
                   ----Prints Initial Data to Output----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 7/21/2014     ----
                   ----Modified By:                 ----
This prints the initial data to the output file.
*/
void dataOutput::initial_data_printer (MemHandler data_mem)
{
	ofile << "Molecular Data Loaded -- Number of Atoms = " << data_mem.ipt_parms.num_atoms << "\n\n";
	ofile << "|-------------------ATOMIC INPUT DATA----------------------|\n";
	for (int i = 0;i < data_mem.ipt_parms.num_atoms; ++i)
        {
                ofile << "ATOM(" << i << ") ";
		ofile << "Atomic Class: " << data_mem.atom_data[i].AtomLetter().c_str() << "\n";
		if (data_mem.atom_data[i].AtomLetter().compare("error") == 0)
		{
			ostringstream error;
			error << "Error - Atomic number " << data_mem.atom_data[i].atomic_num << " is unsupported by this software.";
			Prog_err_quit(error.str(),data_mem);
		}
                ofile << "Cartesian Coords: X(" << data_mem.atom_data[i].pos_xyz[0] << ") ";
                ofile << "Y(" << data_mem.atom_data[i].pos_xyz[1] << ") ";
                ofile << "Z(" << data_mem.atom_data[i].pos_xyz[2] << ") \n";

                //data_mem.atom_data[i].CartToSphere();

                //ofile << "Spherical Coords: R(" << data_mem.atom_data[i].pos_sph[0] << ") ";
                //ofile << "phi(" << data_mem.atom_data[i].pos_sph[1] << ") ";
                //ofile << "theta(" << data_mem.atom_data[i].pos_sph[2] << ") \n";
       	}
	ofile << "|-------------------ATOMIC INPUT DATA----------------------|\n\n";
};

/*____________________________________________________________________________
                   ----Error Printer and Ender      ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 7/21/2014     ----
                   ----Modified By:                 ----
This prints a closing error message and ends the program. 
*/
void dataOutput::Prog_err_quit(string Message,MemHandler data_mem)
{
	ofile << "\n" << Message << "\n";
	ofile << "Freeing Memory Allocations...\n";
	data_mem.mem_free(0);
	close_output(1);
	exit(1);
};
