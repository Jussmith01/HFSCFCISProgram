#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../utils_cpp/lib_includes.h"
#include "classes.h"

// ********************************************************************* //
// ******************memhandler CLASS MEMBER FUNCTIONS****************** //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Convert Distances in Atoms ----
                   ----for input data.            ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/24/2014   ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
void MemHandler::ConvertDistances ()
{
	for (int i = 0; i < ipt_parms.num_atoms;++i)
	{
		atom_data[i].AngToAU();
	}
}

/*____________________________________________________________________________
                   ----Allocates Necessary Memory ----
		   ----for input data.		  ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 7/24/2014   ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
void MemHandler::mem_alloc ()
{
	atomsMemSize = ipt_parms.num_atoms * sizeof(atoms);
	//cout << "ATOMS MEM SIZE: " << atomsMemSize << "\n";
	atom_data = (atoms *) malloc (atomsMemSize);
}

/*____________________________________________________________________________
                   ----Frees All Memory Allocations ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/15/2014     ----
                   ----Modified By:                 ----
Frees any memory allocated by this class.
*/
void MemHandler::mem_free (int inc)
{
	free(atom_data);

	if (inc == 1)
	{
        	delete [] orb_idx;
        	orb_idx = NULL;

                delete [] atom_idx;
                atom_idx = NULL;
	}
}

/*____________________________________________________________________________
                   ----Determination of Basis       ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/15/2014     ----
                   ----Modified By:                 ----
Determines the number, but not type, of basis functions required for the memory
allocation in the scf class.
*/
void MemHandler::num_basis_funcs ()
{
	orbs = 0;
	for (int i = 0; i < ipt_parms.num_atoms; ++i)
	{
	        switch (atom_data[i].atomic_num)
        	{
        	        case 1: {orbs += 1;break;}// 1x1S
        	        case 2: {orbs += 1;break;}// 1x1S
        	        case 3: {orbs += 5;break;}// 1x1S,1x2S
        	        case 4: {orbs += 5;break;}// 1x1S,1x2S
        	        case 5: {orbs += 5;break;}// 1x1S,1x2S,3x2P
        	        case 6: {orbs += 5;break;}// 1x1S,1x2S,3x2P
        	        case 7: {orbs += 5;break;}// 1x1S,1x2S,3x2P
        	        case 8: {orbs += 5;break;}// 1x1S,1x2S,3x2P
        	        case 9: {orbs += 5;break;}// 1x1S,1x2S,3x2P
        	        case 10: {orbs += 5;break;}// 1x1S,1x2S,3x2P
        	}	
	}
}

/*____________________________________________________________________________
                   ----Define Basis Ident Array     ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/15/2014     ----
                   ----Modified By:                 ----
Defines the basis function identification array called orb_idx.
*/
void MemHandler::set_orb_idx (dataOutput *optfile)
{
        orb_idx = new unsigned int [orbs];
        atom_idx = new unsigned int [orbs];
	int j = 0;
	num_elec = 0;
	optfile->ofile << "|---------------------Orbitals-------------------|\n";
	for (int i = 0; i < ipt_parms.num_atoms; ++i)
        {
		num_elec += atom_data[i].atomic_num;
                switch (atom_data[i].atomic_num)
                {
                        case 1: 
			{
				orb_idx[j] = 1;
				atom_idx[j] = i;
				++j;
				optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): " << basis_type(1).c_str() << "\n";
				break;
			}// 1x1S
                        case 2:
                        {
                                orb_idx[j] = 1;
				atom_idx[j] = i;
                                ++j;
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): " << basis_type(1).c_str() << "\n";
                                break;
                        }// 1x1S
                        case 3:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                        case 4:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                        case 5:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                        case 6: 
			{
				optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
				for (int k = 1; k <= 5; ++k)
				{
					orb_idx[j] = k;
					atom_idx[j] = i;
					optfile->ofile << basis_type(k).c_str() << " ";
					++j;
				}
				optfile->ofile << "\n";
				break;
			}// 1S,2S,2Px,2Py,2Pz
                        case 7:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                        case 8:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                        case 9:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                        case 10:
                        {
                                optfile->ofile << "Orbitals for atom " << i << " (" << atom_data[i].AtomLetter().c_str() << "): ";
                                for (int k = 1; k <= 5; ++k)
                                {
                                        orb_idx[j] = k;
                                        atom_idx[j] = i;
                                        optfile->ofile << basis_type(k).c_str() << " ";
                                        ++j;
                                }
                                optfile->ofile << "\n";
                                break;
                        }// 1S,2S,2Px,2Py,2Pz
                }
        }
	optfile->ofile << "|---------------------Orbitals-------------------|\n\n";
	
	num_elec = num_elec - ipt_parms.sys_charge;
	optfile->ofile << "\nNumber of Electrons: " << num_elec << "\n\n";
}

/*____________________________________________________________________________
                   ----Return Basis Type	    ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/15/2014     ----
                   ----Modified By:                 ----
Given an integer, returns the basis type corresponding. Used for output.
*/
string MemHandler::basis_type (int value)
{
	string rtnstr;
	switch (value)
        {
                case 1: {rtnstr = "1S";break;}// 1S
                case 2: {rtnstr = "2S";break;}// 2S
                case 3: {rtnstr = "2Px";break;}// 2Px
                case 4: {rtnstr = "2Py";break;}// 2Py
                case 5: {rtnstr = "2Pz";break;}// 2Pz
        }
	return rtnstr;
}

