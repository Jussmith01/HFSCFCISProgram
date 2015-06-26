#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../utils_cpp/lib_includes.h"
#include "classes.h"

// ********************************************************************* //
// ******************FOUR E- INT CLASS MEMBER FUNCTIONS****************** //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Allocates Necessary Memory ----
		   ----for input data.		  ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 11/11/2014   ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
int Four_Elec_Index::alloc_val_table(int orbitals)
{
	//Set Number of Degenerate Electron Repulsion integrals and the index
	int fset = orbitals * (orbitals + 1) / (int) (2);
	int degenerate = fset * (fset + 1) / (int) (2);

	/*Val_Table = new Val_Store [degenerate];
	for (int m = 0; m < degenerate; ++m)
		{Val_Table[m].key = 0;}
	*/
	//num_deg = degenerate;
	num_val = degenerate;
	return degenerate;
};

/*____________________________________________________________________________
                   ----Free Memory		  ----
                   ----for input data.            ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 11/11/2014  ----
                   ----Modified By:               ----
*/
void Four_Elec_Index::free_val_table()
{
	//delete [] Val_Table;
        //Val_Table = NULL;
	values.clear();
};

/*____________________________________________________________________________
                   ----Set Value at given Index   ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/16/2014  ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
void Four_Elec_Index::set_val(double val,int i,int j,int k,int l)
{
	std::string key = produce_unique_ident(i,j,k,l);

	/*int comp = 1;
	int m = 0;

	while (comp == 1)
	{
		int pos = (key + m) % num_val;
		//comp = Val_Table[pos].compare_key(key);
		if (Val_Table[pos].key == 0)
		{
			//cout << "pos: " << pos << "\n";
			Val_Table[pos].key = key;
			Val_Table[pos].value = val;
			break;
		}
		++m;
	}*/

	values[key.c_str()]=val;
};

/*____________________________________________________________________________
                   ----Set Value at given Index   ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/16/2014  ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
double Four_Elec_Index::get_val(int i,int j,int k,int l)
{
        /*int pos = get_mem_loc(i,j,k,l);

	return Val_Table[pos].value;*/

        std::string key = produce_unique_ident(i,j,k,l);
	double val=values[key.c_str()];

	return val;
};

/*____________________________________________________________________________
                   ----Get Memory Location        ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/16/2014  ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
/*int Four_Elec_Index::get_mem_loc(int i,int j,int k,int l)
{
        long int key = produce_unique_ident(i,j,k,l);

	int rtnval;
        int comp = 1;
        int m = 0;

        while (comp == 1)
        {
                int pos = (key + m) % num_val;
                //comp = Val_Table[pos].compare_key(key);
                if (Val_Table[pos].compare_key(key) == 0)
                {
			rtnval = pos;
                        break;
                }
                ++m;
        }

	//cout << "NUMSKIP: " << m << "\n";
	return rtnval;
};*/

/*____________________________________________________________________________
                   ----Return Ident               ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/16/2014  ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
long int Four_Elec_Index::produce_ident(int i,int j,int k,int l)
{
        long int key;

        string one = "1";
        string zero = "0";

        string svali = int_to_str(i);
        string svalj = int_to_str(j);
        string svalk = int_to_str(k);
        string svall = int_to_str(l);

        int leng = svali.length();
        for (int m = 0;m < 3 - leng;++m)
                {svali = zero + svali;}

        leng = svalj.length();
        for (int m = 0;m < 3 - leng;++m)
                {svalj = zero + svalj;}

        leng = svalk.length();
        for (int m = 0;m < 3 - leng;++m)
                {svalk = zero + svalk;}

        leng = svall.length();
        for (int m = 0;m < 3 - leng;++m)
                {svall = zero + svall;}

        string val_str = one + svali + svalj + svalk + svall;
        long int rtnval;

        istringstream(val_str) >> key;

	return key;
};

/*____________________________________________________________________________
                   ----Produce Unique Integral Ident ----
                   ----for input data.               ----
                   ----Author: Justin Smith          ----
                   ----Date Modified: 11/11/2014     ----
                   ----Modified By:                  ----

*/
std::string Four_Elec_Index::produce_unique_ident (int i,int j,int k,int l)
{
                int tmpi,tmpi2;
                int tmpj,tmpj2;
                int tmpk,tmpk2;
                int tmpl,tmpl2;

                int deg_val;


                //cout << "INDEX(" << count << "): " << i << j << k << l;

                if (j < i)
                {
                        tmpi = j;
                        tmpj = i;
                } else {
                        tmpi = i;
                        tmpj = j;
                }

                if (l < k)
                {
                        tmpk = l;
                        tmpl = k;
                } else {
                        tmpk = k;
                        tmpl = l;
                }

                if (tmpk < tmpi)
                {
                        tmpi2 = tmpk;
                        tmpj2 = tmpl;
                        tmpk2 = tmpi;
                        tmpl2 = tmpj;
                } else {
                        if (tmpk == tmpi && tmpl < tmpj)
                                {
                                        tmpi2 = tmpk;
                                        tmpj2 = tmpl;
                                        tmpk2 = tmpi;
                                        tmpl2 = tmpj;
                                } else {
                                        tmpi2 = tmpi;
                                        tmpj2 = tmpj;
                                        tmpk2 = tmpk;
                                        tmpl2 = tmpl;
                                }
                }

		long int rtnval = produce_ident(tmpi2,tmpj2,tmpk2,tmpl2);
		std::stringstream ss;
		ss << rtnval;

                //cout << "i: " << i << " j: " << j << " k: " << k << " l: " << l << " RETURN:" << rtnval << "\n";
                return ss.str();
}

