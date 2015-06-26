#ifndef cuda_includes
#define cuda_includes
//C/C++ includes

//Custom includes
#include "../utils_cpp/lib_includes.h"
#include <unordered_map>

// ********************************************************************* //
// ***************************DEFINE CLASSES**************************** //
// ********************************************************************* //

//________________________________________________________________________//
//      *************************************************************     //
//                      Program Memory Handler Class
//      This class hold all needed memory handling declarations and
//	actions, such as memory allocations, transfers, to and from 
//	device and host, in the form of class member functions.
//      *************************************************************     //
class MemHandler
{
	public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
		//Parameters handler
		parameters ipt_parms;
	
		//Molecular data handler
		size_t atomsMemSize;
		atoms *atom_data;

		//Number of Orbitals
		unsigned int orbs;

		//Orbital index array -- gives order to the matrices defined later
		unsigned int *orb_idx;//Indexes the type of orbitals
		unsigned int *atom_idx;//Indexes which orbitals belong to which atom_data[] value

		int num_elec;

                //Electronic Energy Storage
                double Eelec;
		double Enuc;

        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
		
		//Allocated Memory
		void mem_alloc(void);

		//Free All Memory
		void mem_free(int inc);

		//Convert input units to au
		void ConvertDistances (void);
	
		//Determine Number of Basis Functions
		void num_basis_funcs(void);

		//Set the orb_idx array
		void set_orb_idx(dataOutput *optfile);

		//Returns the orbital type given value
		string basis_type(int value);
};

/*____________________________________________________________________________
                   ----Store STOnG Constant Data    ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/15/2014    ----
                   ----Modified By:                 ----
*/
//      Class designed to store the orbital coefficient and contraction
//	coefficient values for the given value of N primative gaussians.
//	This data is stored in an external file named STOnG_constants 
//	in a directory called bs_constant where the executable is
//	located. Declare STOnG_handler as an array, where the index
//	is the atomic number. If atom only has 1x1S then norbs = 1,
//	if 1x1S, 1x2S then norbs = 2 ... if it also contains all
//	three 2P orbtials then only set norbs = 3. This is since all
//	P orbitals have the same coefficients.
//      *************************************************************     //
class STOnG_handler
{
        public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
		//Number of contracted gaussians in basis
		int N;

		//Number of orbitals for this atom type
		int norbs;

		//Orbital Coefficients
		double *a;

		//Contraction Coefficients
		double *d;

		//Normalization constants
                double norm;

        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
		//Function allocates a and d
		void mem_alloc(void);

		//Function frees a and d
		void mem_free(void);

		//Function to fix string for finding given constants
		string AtomLetter (int an);

		//Function that converts ints to strings
		string int_to_str (int a_value);

		//Function to trim whitespace
		string trim (string line);

		//Get double value from line
		double get_value(string line,int i);

                //Set the orb_idx array
                int rtn_bas_idx(int orb);

		//Function reads the
		void read_input(int NCG,int atomic_number,dataOutput *optfile);
};

//________________________________________________________________________//
//      *************************************************************     //
//                              Hash Class
//      *************************************************************     //
/*class Val_Store
{
        public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
	long int key;
	double value;

	//-----------------------------
        //Public Member Class Functions
        //-----------------------------
	int compare_key(long int key_in)
	{
		int rtnval;

		if (key_in == key)
			{rtnval = 0;} 
		else 
			{rtnval = 1;}

		return rtnval;
	};
};*/

//________________________________________________________________________//
//      *************************************************************     //
//                            Key Builder Class
//      *************************************************************     //
/*class Four_Elec_Index
{
        public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
	Val_Store *Val_Table;
	long int num_val;

        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
	int alloc_val_table(int orbitals);
	long int produce_unique_ident(int i,int j,int k,int l);
	long int produce_ident(int i,int j,int k,int l);
	int get_mem_loc(int i,int j,int k,int l);
	void set_val(double val,int i,int j,int k,int l);
	double get_val(int i,int j,int k,int l);
	void free_val_table(void);
};*/

//________________________________________________________________________//
//      *************************************************************     //
//                            Key Builder Class
//      *************************************************************     //
class Four_Elec_Index
{
        public:
        //-------------------------
        //Public Class Declarations
        //-------------------------
        unordered_map<std::string,double> values;
        long int num_val;

        //-----------------------------
        //Public Member Class Functions
        //-----------------------------
        int alloc_val_table(int orbitals);
        std::string produce_unique_ident(int i,int j,int k,int l);
        long int produce_ident(int i,int j,int k,int l);
        //int get_mem_loc(int i,int j,int k,int l);
        void set_val(double val,int i,int j,int k,int l);
        double get_val(int i,int j,int k,int l);
        void free_val_table(void);
};


#endif
