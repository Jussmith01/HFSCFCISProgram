#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../utils_cpp/lib_includes.h"
#include "func_classes.h"

// ********************************************************************* //
// *********************2D ARRAY STRUCT FUNCTIONS*********************** //
// ********************************************************************* //

/*____________________________________________________________________________
                   ----Return 2D Array Value      ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/14/2014  ----
                   ----Modified By:               ----
Return the value from the row and column given of a 2d array stored in linear
memory. The array is assumed to be row major format.
*/
double array_idx_funcs::rtn_2dval (int row,int col,double *array)
{
	double rtnval = array[col + row * width];
	return rtnval;
}

/*____________________________________________________________________________
                   ----Replace 2D Array Value     ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/14/2014  ----
                   ----Modified By:               ----
Replace the value from the row and column given of a 2d array stored in linear
memory. The array is assumed to be row major format.
*/
void array_idx_funcs::rpl_2dval (int row,int col,double value,double *array)
{
        array[col + row * width] = value;
}

/*____________________________________________________________________________
                   ----Add 2D Array Value         ----
                   ----Author: Justin Smith       ----
                   ----Date Modified: 10/14/2014  ----
                   ----Modified By:               ----
Add the value from the row and column given of a 2d array stored in linear
memory. The array is assumed to be row major format.
*/
void array_idx_funcs::add_2dval (int row,int col,double value,double *array)
{
        array[col + row * width] += value;
}

