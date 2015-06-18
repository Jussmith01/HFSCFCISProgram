#ifndef funcs_includes
#define funcs_includes
//C/C++ includes

//Custom includes
#include "../utils_cpp/lib_includes.h"
#include "../scf_main/scf_classes.h"
#include "../classes/classes.h"

// ********************************************************************* //
// ***************************DEFINE CLASSES**************************** //
// ********************************************************************* //

//________________________________________________________________________//
	/**************************************************************
                     Array Indexing Structure Class
	This struct hold variables and functions used to index arrays
	of various types.
	**************************************************************/
struct array_idx_funcs
{
	public:
        //-------------------------
        //Public Struct Declarations
        //-------------------------
		//Width of working arrays - set externally
		int width;

        //-----------------------------
        //Public Member Struct Functions
        //-----------------------------
		
		//Function that returns the value of the index given
		double rtn_2dval (int row,int col,double *array);

		//Function that returns the value of the index given
		void rpl_2dval (int row,int col,double value,double *array);

		//Function that adds the value of the index given
		void add_2dval (int row,int col,double value,double *array);
};

// ********************************************************************* //
// ***************************DEFINE FUNCTIONS**************************** //
// ********************************************************************* //

void calculate_integrals(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);
void calc_nuclear_repulsion(MemHandler *data_mem,dataOutput *optfile);

#endif
