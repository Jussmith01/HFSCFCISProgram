#ifndef class_header
#define class_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

// ********************************************************************* //
// ***************************DEFINE CLASSES**************************** //
// ********************************************************************* //
//Predeclarations to prevent co-dependency compilation errors
class atoms; 
class parameters;
class output_file_set;
class MemHandler;

//________________________________________________________________________//
//      *************************************************************     //
//				Data Output Class
//      Holds everything needed to output the data too the output files.
//      *************************************************************     //
class dataOutput
{
	//-------------------------
	//Public Class Declarations
	//-------------------------
        public:
		// Declare the output stream:
		ofstream ofile; 

		//Verbosity of output
		int verbose;	
	//-----------------------------
	//Public Member Class Functions
	//-----------------------------

		//Creates/Opens the output file:
		void set_output(char *filename);

		//Closes the output file:
		void close_output(int close_type);

		//Prints the initial input data to the output:
		void initial_data_printer(MemHandler data_mem); //Prints the initial input data from the input

                //Program to print a matrix "A" of size "N" with title "Title"
                void Matrix_printer(string Title,double *A,int N);

                //Program to print an array "A" of size "N" with title "Title" with numPcol per row
		void Array_printer(string Title,double *A,int N,int numPcol);

		void Array_printer_harTOeV(string Title,double *A,int N,int numPcol);
	
		//Program to help with error handling printing. Will cut down on lines of code later.
		void Prog_err_quit(string Message,MemHandler data_mem);

		void Prec_Printer(double A,int N);
};

//________________________________________________________________________//
//      *************************************************************     //
//                              Atomic Data Class
//      Holds the data for atoms along with conversion and transfor-
//	mation functions. An array of these give molecular data.
//      *************************************************************     //
class atoms
{

        //-------------------------
        //Public Class Declarations
        //-------------------------
        public:
		//Holds the atomic number:
        	int atomic_num;

		//Holds the cartesian position vector of the atom:
        	double pos_xyz[3]; // 0 = x; 1 = y; 2 = z;
		
		//Holds the spherical position vector of the atom:
        	double pos_sph[3]; //0 = r; 1 = theta; 2 = phi
	
        //-----------------------------
        //Public Member Class Functions
        //-----------------------------

		//Converts the cartesian vector in pos_xyz to a spherical vector in pos_sph:
		void CartToSphere(void);
		
		//Converts the sperical vector in pos_sph to a cartesian vector in pos_xyz:
		void SphereToCart(void);
		
		//Convert angstrom units of pos_xyz to a.u.
		void AngToAU (void);

		//Returns the periodic element letter code for the given value in atomic_num:
		string AtomLetter(void);
};

//________________________________________________________________________//
//      *************************************************************     //
//                            Parameter Data Class
//      	    Holds any parameters from the input file.
//      *************************************************************     //
class parameters
{
        //-------------------------
        //Public Class Declarations
        //-------------------------
        public:
		//This array holds the input parameters as strings
		//Numerical values must be converted before use
		//Specificed below:
		/*
		  -----------------------------------------------
		  |Array Index   |      Description		|
		  -----------------------------------------------
		  | parmchars[0] |n from STO-nG to nCGO		|
		  |     [1]      |System Charge	to sys_charge	|
		  |     [2]      |Jacobi Convergence to j_conv	|
		  |	[3]      |Energy Convergence to e_conv	|
		  |     [4]      |Output Verbosity to opt_verb	|
		  |     [5]      |CIS Calculation Flag   	|
		  |     [6]      |Transition Dipole Flags	|
		  -----------------------------------------------
		*/
		string parmchars[7];

		string dataout;

		//Holds the total number of atoms in the molecule(declared during file read):
		int num_atoms;

		//Number of contracted gaussian type orbitals in the basis function
		int nCGO;
		int sys_charge;
		int opt_verb;
		int cis_flag;
		int td_flag;
		double j_conv;
		double e_conv;

		//-----------------------------
        	//Public Member Class Functions
	        //-----------------------------
		void define_parms(dataOutput *optfile);
};

//________________________________________________________________________//
//      *************************************************************     //
//                               Timer Class
//                  Holds timer variables and class functions
//      *************************************************************     //
class timer
{
        //--------------------------
        //Private Class Declarations
        //--------------------------

        time_t start_time; //Holds start time
        time_t end_time; //Holds end time
        double run_time; //Holds run time = end_time - start_time
	clock_t t;
	double CLOCK_TIME;

        //------------------------------
        //Private Member Class Functions
        //------------------------------

	//Intakes a value of time in seconds and returns a string formmatted as:
	//Computation Wall Time: d days h hours m minutes s seconds
        string mk_time_string(double time_val); 

	public:
        //-----------------------------
        //Public Member Class Functions
        //-----------------------------

		//Sets the start time for the timer
        	void set_timer(void);

		//Sets the end time for the timer
        	void end_timer(void);

		//Prints run_timer = end_time - start_time 
        	void print_time(string message,ofstream &file1);
		
		//Prints the clock time
		void print_clock_time (string message,ofstream &file1);

                //Prints run_timer = end_time - start_time 
                void print_time_s(string message,ofstream &file1);

		//Resets the timer if needed
        	void reset(void);
};
// ********************************************************************* //
// **********************DEFINE EXTERNAL FUNCTIONS********************** //
// ********************************************************************* //

//This function reads the input file. Call with:
// read_input([Name of input file], &[parameters class declaration], &[a NULL atoms class declaration])
void read_input (char *inputfile,MemHandler *data_mem,dataOutput *optfile);
string int_to_str(int a_value);

#endif

