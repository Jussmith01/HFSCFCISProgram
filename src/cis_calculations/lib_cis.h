#ifndef cis_includes
#define cis_includes
//C/C++ includes

//Custom includes
#include "../utils_cpp/lib_includes.h"
#include "../scf_main/scf_classes.h"
#include "../classes/classes.h"

// --------------------------------------------------------------------- //
// 			DEFINE CIS RELATED FUNCTIONS			 //
// --------------------------------------------------------------------- //

//Carry out a CIS calculation
void cis_energy_calc(MemHandler *data_mem,scf_data *scf,dataOutput *optfile);

#endif

