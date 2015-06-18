#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <iostream>
#include "lib_includes.h"

using namespace std;

string trim (string line)
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

void parameters::define_parms (dataOutput *optfile)
{
	optfile->ofile << "\n|------------------------------------------------------|\n";
	optfile->ofile << "|--------------------Input Parameters------------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n";
	istringstream(parmchars[0]) >> nCGO;
	optfile->ofile << "Basis Selected: STO-" << nCGO << "G\n";
	istringstream(parmchars[1]) >> sys_charge;
	optfile->ofile << "System Charge: " << sys_charge << "\n";
	istringstream(parmchars[2]) >> j_conv;
	optfile->ofile << "Jacobi Algorithm Convergence: " << j_conv << "\n";
	istringstream(parmchars[3]) >> e_conv;
	optfile->ofile << "SCF Energy Convergence: " << e_conv << "\n";
	istringstream(parmchars[4]) >> opt_verb;
	optfile->ofile << "Output Verbosity: " << opt_verb << "\n";
	optfile->verbose = opt_verb;
	istringstream(parmchars[5]) >> cis_flag;
	optfile->ofile << "CIS Calculations Flag(1=y/0=no): " << cis_flag << "\n";
        istringstream(parmchars[6]) >> td_flag;
        optfile->ofile << "Transition Dipoles Flag(1=y/0=no): " << td_flag << "\n";
        dataout = trim(parmchars[7]);
        optfile->ofile << "Dataout Append Name: " << dataout << "\n";
	optfile->ofile << "|------------------------------------------------------|\n";
	optfile->ofile << "|-----------------End Input Parameters-----------------|\n";
	optfile->ofile << "|------------------------------------------------------|\n\n";
}

