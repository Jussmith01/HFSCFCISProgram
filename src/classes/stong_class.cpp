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
                   ----Allocates Necessary Memory ----
		   ----for input data.		  ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/16/2014   ----
                   ----Modified By:               ----
This program allocates all necessary memory for the input data.
*/
void STOnG_handler::mem_alloc ()
{
	
	a = new double [N * norbs];//Size, number of contracted gaussians * number of orbitals
	d = new double [N * norbs];//Size, number of contracted gaussians * number of orbitals
	norm = 1;
}

/*____________________________________________________________________________
                   ----Frees All Memory Allocations ----
                   ----Author: Justin Smith         ----
                   ----Date Modified: 10/16/2014     ----
                   ----Modified By:                 ----
Frees any memory allocated by this class.
*/
void STOnG_handler::mem_free ()
{
        delete [] a;
        a = NULL;

        delete [] d;
        d = NULL;
}

/*____________________________________________________________________________
                   ----Reads Constants Input and Stores ----
                   ----Author: Justin Smith         	----
                   ----Date Modified: 4/21/2014     	----
                   ----Modified By: Justin Smith        ----
Reads the constants from the input and stores them.
*/
void STOnG_handler::read_input (int NCG,int atomic_number,dataOutput *optfile)
{
	N = NCG;

	char* Path = getenv("SCF_HOME");
	string scfpath;
	if (Path!=NULL) {
		scfpath = Path;
	} else {
		cout << "ERROR: SCF_HOME not set - cannot find STO-" << NCG << "G constants.\n";
		optfile->ofile << "ERROR: SCF_HOME not set - cannot find STO-" << NCG << "G constants.\n\n";
		optfile->close_output(1);
		exit(0);
	}

	string STO = "bs_const/STO";
	string pf = "G_constants.dat";

	string inputfile = scfpath + STO + int_to_str(N) + pf;
        string sestr = AtomLetter(atomic_number);

	mem_alloc();
	
	string line;
        ifstream iptfile (inputfile.c_str());
	int counter = 0;
	int start_read = 0;
	int skipper = 0;
	int i_data = 0;
        int i=0;

        //COUNTS NUMBER OF LINES IN THE INPUT FILE FOR DYNAMIC MEMORY ALLOCATION
	//cout << "FILENAME: " << inputfile.c_str() << " Search String: " << sestr << "\n";
        if (iptfile.is_open())
        {
                        while (!iptfile.eof())
                        {
				if (i == norbs * N) {break;}
				
                                getline(iptfile,line);
                                if (trim(line).compare(sestr) == 0) 
				{
					//cout << "LINE: " << trim(line) << " SESTR: " << sestr << " VAL: " << trim(line).compare(sestr) << "\n";
					start_read = counter;
				}

				if (start_read > 0 && counter >= start_read + 3 && counter < start_read + 3 + N)
				{
					a[i] = get_value(line,1);
					//cout << "1S LINE(" << counter << "): " << a[i];
					d[i] = get_value(line,2);
					//cout << " - " << d[i] << "\n";
					++i;
				}

                                if (start_read > 0 && counter >= start_read + 3 + N + 1 + skipper && counter < start_read + 3 + N + N + 1)
                                {
                                        a[i] = get_value(line,1);
					//cout << "2S LINE(" << counter << "): " << a[i];
                                        d[i] = get_value(line,2);
					//cout << " - " << d[i] << "\n";
                                        ++i;
                                }

                                if (start_read > 0 && counter >= start_read + 3 + N + N + 2 + skipper && counter < start_read + 3 + N + N + N + 2)
                                {
                                        a[i] = get_value(line,1);
					//cout << "2P LINE(" << counter << "): " << a[i];
                                        d[i] = get_value(line,2);
					//cout << " - " << d[i] << "\n";
                                        ++i;
                                }

				++counter;
                        }
        iptfile.close();

        } else {
		optfile->ofile << "***error -- unable to open STO constant file***\n";
		optfile->ofile << "FILENAME: " << inputfile.c_str() << "\n";
		optfile->close_output(1);
		exit(1);
	}
}

/*____________________________________________________________________________
                ----Atomic Letter Code Selector Program----
                ----Author: Justin Smith               ----
                ----Date Modified: 10/15/2014          ----
                ----Modified By:                       ----
This atoms class member function selects the proper letter code for the given 
atomic number defined. The main purpose of this is simply for output, but 
also, any unsupported atoms should not be included here. This is so that the 
unsupported atom error is thrown and the program is stopped at this point with 
the correct error message.
*/
string STOnG_handler::AtomLetter (int an)
{
        string Letter;
        switch (an)
        {
                case 1: {Letter = "h";norbs = 1;break;}
                case 2: {Letter = "he";norbs = 1;break;}
                case 3: {Letter = "li";norbs = 3;break;}
                case 4: {Letter = "be";norbs = 3;break;}
                case 5: {Letter = "b";norbs = 3;break;}
                case 6: {Letter = "c";norbs = 3;break;}
                case 7: {Letter = "n";norbs = 3;break;}
                case 8: {Letter = "o";norbs = 3;break;}
                case 9: {Letter = "f";norbs = 3;break;}
                case 10: {Letter = "ne";norbs = 3;break;}
                default:
                {
                Letter = "error";
                break;
                }
        }
	string STO = "-STO-";
	string G = "G";
	string rtnstr = Letter + STO + int_to_str(N) + G;

        return rtnstr;
}

/*____________________________________________________________________________
                ----Changes and integer to a string    ----
                ----Author: Justin Smith               ----
                ----Date Modified: 10/15/2014          ----
                ----Modified By:                       ----
Converts an integer to a string.
*/
string STOnG_handler::int_to_str (int a_value)
{
	const int n = 8;
        ostringstream out;

        if (a_value <= 1.0e-04 && a_value >= -1.0e-04)
        {
                int nval = 0;
                out << setprecision(n) << nval;
        } else {
                out << setprecision(n) << a_value;
        }

        return out.str();
}

/*____________________________________________________________________________
                ----Trims Whitespace from a string     ----
                ----Author: Justin Smith               ----
                ----Date Modified: 10/15/2014          ----
                ----Modified By:                       ----
Trims a string.
*/
string STOnG_handler::trim(string line)
{
        string trimmed_line;
        unsigned int first = line.find_first_not_of(" \t\r\n\x0b");
        unsigned int last = line.find_last_not_of(" \t\r\n\x0b") + 1;

        last = last - first;

        if (line.length() != 0)
        {
                trimmed_line = line.substr(first,last);
        }
        return trimmed_line;
}

/*____________________________________________________________________________
                ----Obtains Float Value from String    ----
                ----Author: Justin Smith               ----
                ----Date Modified: 10/15/2014          ----
                ----Modified By:                       ----
Returns double value from a line
*/
double STOnG_handler::get_value(string line,int i)
{
	double rtnval;
	string tmp1;

	int pos = 0;

	tmp1 = line;

	pos = tmp1.find_first_of("-.0123456789");
	tmp1 = tmp1.substr(pos);

	pos = tmp1.find_first_not_of("-.0123456789");

	if (i == 1)
	{
		tmp1 = tmp1.substr(0,pos + 1);
		istringstream(tmp1) >> rtnval;
	}

	if (i == 2) 
	{
		tmp1 = tmp1.substr(pos + 1);
		pos = tmp1.find_first_of("-.0123456789");
		tmp1 = tmp1.substr(pos);
		istringstream(tmp1) >> rtnval;
	}

	return rtnval;
}

/*____________________________________________________________________________
                ----Return Basis Index		       ----
                ----Author: Justin Smith               ----
                ----Date Modified: 10/15/2014          ----
                ----Modified By:                       ----
Returns the basis function index value so that correct basis function constants
can be selected easily.
*/
int STOnG_handler::rtn_bas_idx (int orb)
{
        int rtnval;
        switch (orb)
        {
                case 1: {rtnval = 1;break;}
                case 2: {rtnval = 2;break;}
                case 3: {rtnval = 3;break;}
                case 4: {rtnval = 3;break;}
                case 5: {rtnval = 3;break;}
                default:
                {
                cout << "error";
                break;
                }
        }
        return rtnval;
}

