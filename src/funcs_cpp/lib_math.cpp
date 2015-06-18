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
//                 Calc Coefficients of the Binomial Expansion 
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/27/2014  ----
                   ----Modified By:               ----
*/
extern double fi(int in,int l1,int l2,double x1, double x2)
{
        double rtnval;

        if (l1 == 0 && l2 == 0) {

                rtnval = 1.0;

        } else if (l1 == 0 && l2 == 1) {
                switch (in)
                {
                        case 0: {rtnval = x2; break;}
                        case 1: {rtnval = 1.0; break;}
                }
        } else if (l1 == 1 && l2 == 0) {
               switch (in)
                {
                        case 0: {rtnval = x1; break;}
                        case 1: {rtnval = 1.0; break;}
                }
        } else if (l1 == 1 && l2 == 1) {
               switch (in)
                {
                        case 0: {rtnval = x1 * x2; break;}
                        case 1: {rtnval = x1 + x2; break;}
                        case 2: {rtnval = 1.0; break;}
                }
        } else {cout << "!!!ERROR: UNSUPPORTED ANGULAR VALUE!!!\n";exit(0);}

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//             Calculates the Coefficient of the NucAtt Factor
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/
extern double C_fact(double an,double bm)
{
        double rtnval;

        rtnval = (2 * M_PI) / (double)(an + bm);

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                      Calc Expansion Term for Fnu(x)
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/
double Diff(int nu,int n,double x)
{
	double rtnval = 1;
	if (n <= 5)
	{
		rtnval = pow(-x,n) / (double)(fact(n) * (2 * (nu + n) + 1));
	} else {
		for (int i = n + 1; i --> 1; )
		{
			rtnval = rtnval * (x / (double)i);
			//cout << "n: " << n << " i: " << i << "\n"; 
		}
		rtnval = rtnval * pow(-1,n) / (double)(2 * (nu + n) + 1);
	}
	return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                      Fnu(x) Function Approximation
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/
extern double Fnu(int nu,double x)
{
        double rtnval;
	double diff = 1;
	double eps = 1.0E-10;
	double xD = (double)x;
	int n = 0;

	//cout << "xD: " << xD << " x: " << x << "\n";

	if (xD < 1.0E-4)
	{
		//cout << "         ***Method X < 1.0E-5 -- Analytical***\n";
		rtnval = 1 / (double)(2 * nu + 1);
	} else if (xD <= 14 && x > 1.0E-4)
	{
		//cout << "         ***Method X < 10 and  > 1.0E-5 -- Numerical Expansion***\n";
		double pos = 0;
		double neg = 0;
		double pos_vals[50];
		double neg_vals[50];
		while (abs(diff) > eps || n % 2 != 0)
		{
			diff = Diff(nu,n,xD);

			if (rtnval != rtnval) {cout << "NaN Detected -- EXITING\n";exit(0);}			

			if (diff > 0)
			{
				pos_vals[n / 2] = diff;
			}
			else if (diff < 0)
			{
				neg_vals[(n - 1) / 2] = abs(diff);
			}
			//cout << setprecision(15) << "            (" << n << ")rtnval: " << rtnval << " diff: " << diff << " x: " << x << " nu: " << nu <<"\n";
			++n;
			if(n == 100) {break;}
		}

		//cout <<  "|--------------------------|\n";
		for (int i = n / 2; i --> 0;)
		{
			pos += pos_vals[i];
			neg += neg_vals[i];
			//cout << setprecision(15) << "            (" << i << ")pos: " << pos << " neg: " << neg << " x: " << x << " n: " << n <<"\n";
		}
		

		rtnval = pos - neg;
	} else if (xD > 14) {
                //cout << "         ***Method X > 10 -- Analytical Equation***\n";
		rtnval = 0.5 * sqrt(M_PI) * pow(x,-0.5);
                
		for (int i = 0; i < nu; ++i)
                {
                        rtnval = (1 / (double)(2 * xD)) * ((2 * i + 1) * rtnval - exp(-xD));
                }

	} else {cout << "**ERROR, RANGE NOT COVERED IN Fm(x)**";}
	
	double testrtn = (double)rtnval;
	//cout << "testrtn = " << testrtn << " rtnval: " << rtnval << " xD: " << xD << " nu: " << nu << "\n";
	//cout << "***FNISHED AFTER " << n << " CYCLES***\n\n";

        return testrtn;
}

//________________________________________________________________________//
//      *************************************************************     //
//                      Error Function Approximation
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/22/2014  ----
                   ----Modified By:               ----
*/
extern double erf(double x,int n)
{
        double rtnval = 0;
	int i = 0;
	double eps = 1.0E-10;
	double diff =100;

        rtnval = 2/(double)sqrt(M_PI);
        while (abs(diff) > eps)
        {
		diff = (pow(-1,i) * pow(x,2 * i + 1)) / (double)(fact(i) * (2 * n + 1));
		rtnval += diff;
		//cout << "rtnval: " << rtnval << "\n";
		++i;
		if(i == 100) {break;}
	}

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                      Normalization Constant
//      
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/20/2014  ----
                   ----Modified By:               ----
*/

extern double norm(int ax,int ay,int az,double an)
{
        double rtnval;

        rtnval = sqrt((pow(2,2 * (ax + ay + az) + 1.5) * pow(an,ax + ay + az + 1.5)) / (double)(d_fact(2 * ax - 1) * d_fact(2 * ay - 1) * d_fact(2 * az - 1) * pow(M_PI,1.5)));

        return rtnval;
}

/*extern double norm(int ax,int ay,int az,double an)
{
        double rtnval;

        rtnval = powf((2 * an)/(double)M_PI,0.75) * ((powf(4 * an,(ax + ay + az) / (double)2)) / (double)powf(d_fact(2 * ax - 1) * d_fact(2 * ay - 1) * d_fact(2 * az - 1),0.5));

        return rtnval;
}*/

extern double norm2(int ax,int ay,int az,double an)
{
        double rtnval;

        rtnval = pow((2 * an)/(double)M_PI,0.75) * sqrt((pow(8 * an, ax + ay + az) * fact(ax) * fact(ay) * fact(az)) / (double)(fact(2 * ax) * fact(2 * ay) * fact(2 * az)));
	//cout << "norm2: " << rtnval << "\n";
        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                   E Factor from The Gaussian Prod Rule     
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/15/2014  ----
                   ----Modified By:               ----
*/
extern double E_factor(double an, double bm, double AmB)
{
        double rtnval;

        rtnval = -((an * bm)/(double)(an + bm)) * AmB;
	//cout << "         RTNVAL: " << rtnval << " an: " << an << " bm: " << bm << " AmB: " << AmB;
	rtnval = exp((double)rtnval);
	//cout << " exprtn: " << rtnval << "\n";

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                     Combination Function (n choose k)
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/15/2014  ----
                   ----Modified By:               ----
*/
extern double nCk(int n, int k)
{
        double rtnval;
        
	rtnval = fact(n) / (fact(k) * fact(n - k));

	return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                           Factorial Function
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/15/2014  ----
                   ----Modified By:               ----
*/
extern double fact(int n)
{
        double rtnval = 1;

	unsigned long int diff;
	unsigned long int N = (unsigned long int)n;

        for (int i = 0; i < n; ++i)
        {
		diff = N - i;
        	rtnval = (unsigned long int)rtnval * diff;
        }

        return rtnval;
}

//________________________________________________________________________//
//      *************************************************************     //
//                       Double Factorial Function
//      *************************************************************     //
/*                 ----Author: Justin Smith       ----
                   ----Date Modified: 10/15/2014  ----
                   ----Modified By:               ----
*/
extern double d_fact(int n)
{
	double rtnval = 1;

	if (abs(n) % 2 == 0)
	{
		for (int i = 0; i < n; ++i)
		{
			if (n - 2 * i <= 1) {break;}
			rtnval = rtnval * (n - 2*i);
			//cout << "n-2i: " << n-2*i << " val: " << rtnval << "\n";
		}
	} else {
		for (int i = 0; i < n; ++i)
                {
			if (n - 2 * i <= 1) {break;}
                        rtnval = rtnval * (n - 2*i);
			//cout << "n-2i: " << n-2*i << " val: " << rtnval << "\n";
                }
	}
	return rtnval;
}

