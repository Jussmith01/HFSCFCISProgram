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
//             Calculates the Coefficient of the NucAtt Factor
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/

extern void calc_nuclear_repulsion(MemHandler *data_mem,dataOutput *optfile)
{
	double rtnval = 0;
	for (int i = 0; i < data_mem->ipt_parms.num_atoms; ++i)	
	{
		for (int j = 0; j < i; ++j)	
		{
			//cout << "i: " << i << " j: " << j << "\n";

			double Ax = data_mem->atom_data[i].pos_xyz[0];
			double Ay = data_mem->atom_data[i].pos_xyz[1];
			double Az = data_mem->atom_data[i].pos_xyz[2];

			double Bx = data_mem->atom_data[j].pos_xyz[0];
			double By = data_mem->atom_data[j].pos_xyz[1];
			double Bz = data_mem->atom_data[j].pos_xyz[2];

			double Za = data_mem->atom_data[i].atomic_num;
			double Zb = data_mem->atom_data[j].atomic_num;

			double r = sqrt((Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz));

			rtnval += (Za * Zb) / (double)r;
		}
	}
	data_mem->Enuc = rtnval;
}
