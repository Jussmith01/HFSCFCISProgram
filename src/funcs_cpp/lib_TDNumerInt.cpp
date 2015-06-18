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

//      *************************************************************     //
//                      	Integrand
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/23/2015  ----
                   ----Modified By:               ----
*/
double function(int l, double gamma,double P,double x)
{
        double eval;

        eval = exp(-gamma*(x-P)*(x-P))*x*pow(x-P,l);

        return eval;
};

//      *************************************************************     //
//           		Search the Function Space 
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/23/2015  ----
                   ----Modified By:               ----
*/
void SearchSpace(double &beg,double &end,int l,double gamma,double P)
{
	double nx = -100;
	bool chkfin = true;

	int cntr=0;

	while (chkfin)
	{
		double x = nx + cntr * 0.1;
		double tv = function(l,gamma,P,x);

		if (abs(tv) > 1.0E-5)
		{
			chkfin=false;
			beg = nx + (cntr-1) * 0.1;
		}

		++cntr;
	}

	nx = 100;
	chkfin = true;

	cntr = 0;

        while (chkfin)
        {
                double x = nx - cntr * 0.1;
                double tv = function(l,gamma,P,x);

                if (abs(tv) > 1.0E-5)
                {
                        chkfin=false;
                        end = nx - (cntr-1) * 0.1;
                }

                ++cntr;
        }
};

//      *************************************************************     //
//           Numerically Integrate the Transition Dipole Integral 
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 4/23/2015  ----
                   ----Modified By:               ----
*/
double TDNumerInt(int l,double gamma,double P)
{
	int N = 1500;
	double beg;
	double end;

	SearchSpace(beg,end,l,gamma,P);

	double ds = (end - beg)/N;
	double c = ds/3.0;

	double sum = 0.0;
	
	sum += function(l,gamma,P,beg);
	sum += function(l,gamma,P,beg + (N-1) * ds);

	for (int i=1;i<N-1;++i)
	{
        	double x = beg + i * ds;

        	//EXECUTE SIMPSONS RULE
        	sum += c*(2 * ((i % 2) + 1) * function(l,gamma,P,x));
	}

	//cout << "beg: " << beg << " end: " << end <<" l: " << l << " gamma: " << gamma << " P: " << P << " VAL: " << sum << endl;

        return sum;
};

