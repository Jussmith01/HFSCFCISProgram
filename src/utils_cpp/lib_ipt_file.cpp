/*******************************************************************
************File Designed to Read in XYZ Trajectory File************
********************************************************************
DESIGNED BY: JUSTIN S. SMITH
CONTACT E-MAIL: JSMITH48@CHEM.UFL.EDU

SAMPLE INPUT FILE:

-----------------STARTS BELOW THIS LINE---------------------
$PARAM
3               !Number of contracted gaussians per orbital "n of STO-nG"
1               !Charge on the molecule
1.0E-14	        !SCF Jacobi Convergence
1.0E-6		!SCF Energy Convergence
1        	!Output verbosity (0 minimal, 1 all)
$ENDPARAM
$COORD
        1       -1.12798434        1.545547           2.2132132
        6       -0.1224696362       -3.1174832504      -1.5015549412
        8        2.1012819183       4.0397895568      -1.1329875333
        6        3.2785162937       3.3434825848      -0.7159896664
        6       14.1137391171      14.5455261160      -3.6515096380
$ENDCOORD
$END OF FILE$
-------------------ENDS ABOVE THIS LINE-----------------------
*********************************************************************/

//**C/C++ included libraries**
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>

//**Included custom headers**
#include "lib_includes.h"
#include "../classes/classes.h"

using namespace std;

//*****************************************************************
//***********Obtains line number of file line delimiters***********
//*****************************************************************
extern string int_to_str(int a_value)
{
        ostringstream out;

        if (a_value <= 1.0e-04 && a_value >= -1.0e-04)
        {
                int nval = 0;
                out << nval;
        } else {
                out << a_value;
        }

        return out.str();
}

//*****************************************************************
//***********Obtains line number of file line delimiters***********
//*****************************************************************
unsigned int get_line_num (const char *search_string, string *INPUT_DATA, int NUM_LINES)
{
	unsigned int LINE_NUM=0,i=0;
	while (i < NUM_LINES) 
	{
		if (strcmp(search_string,INPUT_DATA[i].c_str()) == 0) 
		{LINE_NUM = i;break;}
		i++;
	}
	return LINE_NUM;
}
//****************************************************************
//******Save coordinates and atomic number to atoms class*********
//****************************************************************
void save_coords_xyz (string *INPUT_DATA,MemHandler *data_mem,int LINE_BEG,int LINE_END)
{	
	int beg_sstr,end_sstr,first,last,sz_substr,i=0,error_chk;
	string comp_str = "";
	for (int j = LINE_BEG + 1; j < LINE_END; ++j)
        {
                error_chk = INPUT_DATA[j].find_first_not_of("-.0123456789 \t");
                if(error_chk >=  0 && error_chk <= INPUT_DATA[j].length())
                {
                        cout << "INPUT FILE DISCREPENCY DETECTED IN $COORD:\n";
                        cout << "POS OF DISCREPENCY - COL(" << error_chk << ") - LINE(" << j << ") - ATOM(" << i <<")\n";
                        cout << "***error in input file -- input line contains characters other than numeric characters -- EXIT***\n";
                        exit(1);
                }

	
		//SAVE THE ATOM TYPE BY ATOMIC NUMBER
		first = INPUT_DATA[j].find_first_of("-.1234567890");
		last = first + INPUT_DATA[j].substr(first).find_first_not_of("-.1234567890");
		error_chk = first + INPUT_DATA[j].substr(first,last).find_first_of("-.");
		if(error_chk >= first && error_chk <= last)
		{
			cout << "INPUT FILE DISCREPENCY DETECTED IN $COORD:\n";
			cout << "POS OF DISCREPENCY - COL(" << error_chk << ") - LINE(" << j << ") - ATOM(" << i <<")\n";
			cout << "***error in input file -- atomic number not an unsigned integer -- EXIT***\n";
			exit(1);
		}
		istringstream(INPUT_DATA[j].substr(first,last)) >> data_mem->atom_data[i].atomic_num;
		//cout << "ATOMIC NUMBER BEFORE SAVE:" << INPUT_DATA[j].substr(first,last) << "\n";
		//cout << "ATOMIC NUMBER AFTER SAVE:" << data_mem->atom_data[i].atomic_num << "\n";

		//SAVE THE XYZ POSITION COORDINATES
		for (int k = 0; k < 3; k++)
		{
			first = last + INPUT_DATA[j].substr(last).find_first_of(".-1234567890");
        	        last = first + INPUT_DATA[j].substr(first).find_first_not_of(".-1234567890");
			istringstream(INPUT_DATA[j].substr(first,last)) >> data_mem->atom_data[i].pos_xyz[k];
		}

		++i;
        }
}

//*******************************************************************
//**********************Read the input file**************************
//*******************************************************************
extern void read_input (char *inputfile,MemHandler *data_mem,dataOutput *optfile) 
{

	//ipt_data *DATA;
	string line;
	ifstream iptfile (inputfile);
	int i=0,NUM_LINES=0;

	//COUNTS NUMBER OF LINES IN THE INPUT FILE FOR DYNAMIC MEMORY ALLOCATION
	if (iptfile.is_open())
        	{
                	while (!iptfile.eof())
                	{
				getline(iptfile,line);
				if (line == "$END OF FILE$") {break;}
				NUM_LINES++;
	                }
        iptfile.close();

        } else {optfile->ofile << "***error -- unable to open file***\n";exit(1);}

	//ALLOCATION OF TEMP DATA STRING ARRAY
	string *DATA = NULL;
	DATA = new string [NUM_LINES];
	
	//READ INPUT LINES INTO STRING ARRAY: DATA
	ifstream ipt2file (inputfile);
	if (ipt2file.is_open())
	{
		while (!ipt2file.eof())
		{
			getline(ipt2file,line);
			if (line != "$END OF FILE$")
			{
				DATA[i] = line;
			} else {break;}
			i++;
		}
	ipt2file.close();

	} else {optfile->ofile << "***error -- unable to open file***\n";exit(1);}
	//PARSE DATA FUNCTIONS FOLLOW

	unsigned int LINE_BEG, LINE_END,char_search=0;
	string test_string_beg,test_string_end;
	i = 0;
	while (i <= 2) 
	{
		switch (i)
		{
			case 0: {test_string_beg = "$PARAM";test_string_end = "$ENDPARAM";break;}
			case 1: {test_string_beg = "$COORD";test_string_end = "$ENDCOORD";break;}
		}

                LINE_BEG = get_line_num(test_string_beg.c_str(),DATA,NUM_LINES);
                LINE_END = get_line_num(test_string_end.c_str(),DATA,NUM_LINES);

		switch (i)
		{
			case 0:
			{
				for (int j = LINE_BEG + 1; j < LINE_END; ++j)
				{
					//cout << "LINE(" << j << "): " << DATA[j] << "\n";
					char_search = DATA[j].find("!");
					data_mem->ipt_parms.parmchars[j-LINE_BEG-1] = DATA[j].substr(0,char_search);
					//cout << data_mem->ipt_parms.parmchars[j-LINE_BEG-1] << "\n";
				}
				break;
			}
                        case 1:
                        {
                                data_mem->ipt_parms.num_atoms = (LINE_END - LINE_BEG) - 1;
				data_mem->mem_alloc();
				save_coords_xyz(DATA,data_mem,LINE_BEG,LINE_END);
				break;
                        }
		}
		i++;
	}

	//***********CLEAR FUNCTION MEMORY ALLOCATIONS***********
	delete [] DATA;
	DATA = NULL; 
} 
